#!/usr/bin/env python3
"""Local instrument export server used by docker-compose for demos.

Serves deterministic paired-end FASTQ gzip files per ``run_id`` with a
manifest (paths, lengths, strong ETags, SHA-256), honours ``Range`` and
``If-Match``, and can simulate the three failure modes the acquisition
client must survive:

* **断流 (stream cut)**: a per-run number of responses are truncated
  mid-body and the connection is closed early.
* **换版 (version change)**: ``POST /admin/rotate`` bumps the global
  version, which changes every file's content, ETag, and SHA-256.
* **忽略 Range (ignore Range)**: a per-run flag makes file GETs always
  answer ``200`` with the full body, ignoring the ``Range`` header.

The run id ``demo`` starts with two truncated responses built in, so the
/docs walkthrough (import -> resume -> resume) works out of the box.  All
faults can also be set at runtime through the admin endpoints.  Standard
library only.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import random
import re
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse

PORT = int(os.environ.get("INSTRUMENT_PORT", "9000"))
SEQ_LEN = 150
BASE_RECORDS = 8000
# Extra records per version bump so a rotation also changes file lengths.
RECORDS_PER_VERSION = 137

_lock = threading.Lock()
_state = {
    "version": 1,
    "cache": {},  # (run_id, version, mate) -> gzip bytes
    # run_id -> {"ignore_range": bool, "cuts_remaining": int,
    #            "cut_fraction": float, "corrupt": bool}
    "behaviors": {},
}


def _default_behavior(run_id: str) -> dict:
    behavior = {
        "ignore_range": False,
        "cuts_remaining": 0,
        "cut_fraction": 0.4,
        "corrupt": False,
    }
    if run_id == "demo":
        # First two file responses are cut mid-stream, then it behaves.
        behavior["cuts_remaining"] = 2
    return behavior


def _behavior(run_id: str) -> dict:
    with _lock:
        behaviors = _state["behaviors"]
        if run_id not in behaviors:
            behaviors[run_id] = _default_behavior(run_id)
        return dict(behaviors[run_id])


def _take_cut(run_id: str) -> bool:
    """Atomically consume one of the run's remaining truncated responses."""

    with _lock:
        behavior = _state["behaviors"].setdefault(run_id, _default_behavior(run_id))
        if behavior["cuts_remaining"] > 0:
            behavior["cuts_remaining"] -= 1
            return True
        return False


def _record_count(run_id: str, version: int) -> int:
    return BASE_RECORDS + (version - 1) * RECORDS_PER_VERSION


def _pick_barcode(rng: random.Random) -> str:
    roll = rng.random()
    if roll < 0.45:
        return "ACGTACGT"
    if roll < 0.90:
        return "TTTTAAAA"
    if roll < 0.95:
        # One-mismatch variant of ACGTACGT (still assigned to it).
        pos = rng.randrange(8)
        replacement = rng.choice([b for b in "ACGT" if b != "ACGTACGT"[pos]])
        return "ACGTACGT"[:pos] + replacement + "ACGTACGT"[pos + 1 :]
    return "CCCCCCCC"  # unassigned


def _generate_mate(run_id: str, version: int, mate: int) -> bytes:
    seed = hashlib.sha256(f"{run_id}|{version}|{mate}".encode()).hexdigest()
    rng = random.Random(seed)
    parts = []
    for i in range(_record_count(run_id, version)):
        read_id = f"READ{i:07d}"
        seq = "".join(rng.choice("ACGT") for _ in range(SEQ_LEN))
        qual = "".join(rng.choice("IIIIJJJJFFFF####") for _ in range(SEQ_LEN))
        if mate == 1:
            header = f"@{read_id}/1 barcode={_pick_barcode(rng)}"
        else:
            header = f"@{read_id}/2"
        parts.append(f"{header}\n{seq}\n+\n{qual}\n")
    return gzip.compress("".join(parts).encode("ascii"), compresslevel=6, mtime=0)


def _mate_bytes(run_id: str, mate: int) -> bytes:
    with _lock:
        version = _state["version"]
        key = (run_id, version, mate)
        cached = _state["cache"].get(key)
    if cached is None:
        cached = _generate_mate(run_id, version, mate)
        with _lock:
            _state["cache"][key] = cached
    return cached


def _file_meta(run_id: str, mate: int) -> dict:
    data = _mate_bytes(run_id, mate)
    digest = hashlib.sha256(data).hexdigest()
    name = f"R{mate}"
    return {
        "name": name,
        "path": f"/runs/{run_id}/files/{name}.fastq.gz",
        "length": len(data),
        "etag": f'"{digest[:32]}"',
        "sha256": digest,
    }


def _manifest(run_id: str) -> dict:
    with _lock:
        version = _state["version"]
    return {
        "run_id": run_id,
        "version": version,
        "files": [_file_meta(run_id, 1), _file_meta(run_id, 2)],
    }


class Handler(BaseHTTPRequestHandler):
    protocol_version = "HTTP/1.1"
    server_version = "InstrumentSource/1.0"

    def log_message(self, fmt, *args):  # keep container logs useful
        print(f"{self.address_string()} - {fmt % args}", flush=True)

    # ---------------- helpers ----------------

    def _send_json(self, status: int, obj: dict, extra_headers: dict | None = None):
        body = json.dumps(obj).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        for key, value in (extra_headers or {}).items():
            self.send_header(key, value)
        self.end_headers()
        self.wfile.write(body)

    def _read_json_body(self) -> dict:
        length = int(self.headers.get("Content-Length") or 0)
        if length <= 0:
            return {}
        try:
            return json.loads(self.rfile.read(length))
        except (ValueError, UnicodeDecodeError):
            return {}

    # ---------------- routing ----------------

    def do_GET(self):
        path = urlparse(self.path).path
        if path == "/health":
            return self._send_json(200, {"status": "ok"})
        match = re.fullmatch(r"/runs/([^/]+)/manifest", path)
        if match:
            return self._send_json(200, _manifest(match.group(1)))
        match = re.fullmatch(r"/runs/([^/]+)/files/(R1|R2)\.fastq\.gz", path)
        if match:
            mate = 1 if match.group(2) == "R1" else 2
            return self._serve_file(match.group(1), mate)
        match = re.fullmatch(r"/admin/runs/([^/]+)/behavior", path)
        if match:
            return self._send_json(200, _behavior(match.group(1)))
        if path == "/admin/state":
            with _lock:
                snapshot = {
                    "version": _state["version"],
                    "behaviors": {k: dict(v) for k, v in _state["behaviors"].items()},
                }
            return self._send_json(200, snapshot)
        return self._send_json(404, {"detail": "not found"})

    def do_POST(self):
        path = urlparse(self.path).path
        if path == "/admin/rotate":
            with _lock:
                _state["version"] += 1
                _state["cache"].clear()
                version = _state["version"]
            return self._send_json(200, {"version": version})
        if path == "/admin/reset":
            with _lock:
                _state["version"] = 1
                _state["cache"].clear()
                _state["behaviors"].clear()
            return self._send_json(200, {"version": 1, "behaviors": {}})
        match = re.fullmatch(r"/admin/runs/([^/]+)/behavior", path)
        if match:
            run_id = match.group(1)
            updates = self._read_json_body()
            with _lock:
                behavior = _state["behaviors"].setdefault(
                    run_id, _default_behavior(run_id)
                )
                if "ignore_range" in updates:
                    behavior["ignore_range"] = bool(updates["ignore_range"])
                if "cuts" in updates:
                    behavior["cuts_remaining"] = max(0, int(updates["cuts"]))
                if "cut_fraction" in updates:
                    behavior["cut_fraction"] = min(
                        0.95, max(0.01, float(updates["cut_fraction"]))
                    )
                if "corrupt" in updates:
                    behavior["corrupt"] = bool(updates["corrupt"])
                snapshot = dict(behavior)
            return self._send_json(200, snapshot)
        return self._send_json(404, {"detail": "not found"})

    # ---------------- file serving ----------------

    def _serve_file(self, run_id: str, mate: int):
        data = _mate_bytes(run_id, mate)
        meta = _file_meta(run_id, mate)
        etag = meta["etag"]
        behavior = _behavior(run_id)

        body = data
        if behavior["corrupt"]:
            # Flip bytes mid-file but keep the manifest's ETag/SHA-256, so
            # the client detects the corruption at digest verification.
            mutated = bytearray(data)
            mutated[len(mutated) // 2] ^= 0xFF
            mutated[len(mutated) // 2 + 1] ^= 0x01
            body = bytes(mutated)

        if_match = self.headers.get("If-Match")
        if if_match is not None:
            candidates = [token.strip() for token in if_match.split(",")]
            if "*" not in candidates and etag not in candidates:
                return self._send_json(
                    412,
                    {"detail": "If-Match failed", "current_etag": etag},
                    extra_headers={"ETag": etag},
                )

        total = len(body)
        start, end = 0, total - 1
        status = 200
        range_header = self.headers.get("Range")
        if range_header and not behavior["ignore_range"]:
            match = re.fullmatch(r"bytes=(\d+)-(\d*)", range_header.strip())
            if not match:
                return self._send_json(
                    416,
                    {"detail": "unsupported Range unit or syntax"},
                    extra_headers={"Content-Range": f"bytes */{total}"},
                )
            start = int(match.group(1))
            if match.group(2):
                end = min(int(match.group(2)), total - 1)
            if start >= total or start > end:
                return self._send_json(
                    416,
                    {"detail": "range not satisfiable"},
                    extra_headers={"Content-Range": f"bytes */{total}"},
                )
            status = 206

        payload = body[start : end + 1]
        self.send_response(status)
        self.send_header("ETag", etag)
        self.send_header("Accept-Ranges", "bytes")
        self.send_header("Content-Length", str(len(payload)))
        if status == 206:
            self.send_header("Content-Range", f"bytes {start}-{end}/{total}")
        self.end_headers()

        if _take_cut(run_id):
            # 断流: promise the full length, deliver only a fraction, then
            # close the connection without completing the body.
            cut_at = max(1, int(len(payload) * behavior["cut_fraction"]))
            self.wfile.write(payload[:cut_at])
            self.wfile.flush()
            self.close_connection = True
            return
        self.wfile.write(payload)


def main() -> None:
    server = ThreadingHTTPServer(("0.0.0.0", PORT), Handler)
    server.daemon_threads = True
    print(f"instrument source listening on :{PORT}", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
