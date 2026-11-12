from __future__ import annotations

import gzip
import hashlib
import io
import json
import random
import re
import sqlite3
import threading
import zipfile
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import pytest
from fastapi.testclient import TestClient

from app import main as main_module

PAIRS = 2000
EXPECTED_COUNTS = {"ACGTACGT": 667, "TTTTAAAA": 667, "unassigned": 666}


def gz(data: bytes) -> bytes:
    return gzip.compress(data, mtime=0)


def rec(read_id: str, barcode: str | None, seq: str, mate: int) -> bytes:
    header = f"@{read_id}/{mate}"
    if barcode is not None:
        header += f" barcode={barcode}"
    return f"{header}\n{seq}\n+\n{'I' * len(seq)}\n".encode("ascii")


def make_run_data(pairs: int = PAIRS, seq_len: int = 150, seed: int = 1234):
    rng = random.Random(seed)
    barcodes = ["ACGTACGT", "TTTTAAAA", "CCCCCCCC"]
    r1_parts, r2_parts = [], []
    for i in range(pairs):
        read_id = f"READ{i:06d}"
        seq1 = "".join(rng.choice("ACGT") for _ in range(seq_len))
        seq2 = "".join(rng.choice("ACGT") for _ in range(seq_len))
        r1_parts.append(rec(read_id, barcodes[i % 3], seq1, 1))
        r2_parts.append(rec(read_id, None, seq2, 2))
    return gz(b"".join(r1_parts)), gz(b"".join(r2_parts))


class StubInstrument:
    """In-process instrument export server with controllable fault modes."""

    def __init__(self, r1: bytes, r2: bytes):
        self.lock = threading.Lock()
        self.bodies = {"R1": r1, "R2": r2}
        self.behavior = {"ignore_range": False, "cuts": 0, "cut_fraction": 0.5}
        self.manifest_mutator = None  # fn(files: list[dict]) applied before send
        self.response_etag_override: str | None = None
        self.lie_range_start = 0  # added to the Content-Range start
        self.lie_range_total: int | None = None  # replaces Content-Range total
        self.requests: list[dict] = []
        stub = self

        class Handler(BaseHTTPRequestHandler):
            protocol_version = "HTTP/1.1"

            def log_message(self, *args):
                pass

            def _json(self, status, obj, headers=None):
                body = json.dumps(obj).encode()
                self.send_response(status)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(body)))
                for key, value in (headers or {}).items():
                    self.send_header(key, value)
                self.end_headers()
                self.wfile.write(body)

            def do_GET(self):
                match = re.fullmatch(r"/runs/([^/]+)/manifest", self.path)
                if match:
                    return self._manifest(match.group(1))
                match = re.fullmatch(
                    r"/runs/([^/]+)/files/(R1|R2)\.fastq\.gz", self.path
                )
                if match:
                    return self._file(match.group(2))
                stub._record(self, 404)
                return self._json(404, {"detail": "not found"})

            def _manifest(self, run_id):
                with stub.lock:
                    files = [stub._meta(run_id, "R1"), stub._meta(run_id, "R2")]
                    mutator = stub.manifest_mutator
                if mutator is not None:
                    mutator(files)
                stub._record(self, 200)
                return self._json(200, {"run_id": run_id, "files": files})

            def _file(self, name):
                with stub.lock:
                    body = stub.bodies[name]
                    etag = stub._etag(name)
                    behavior = dict(stub.behavior)
                    etag_override = stub.response_etag_override
                    lie_start = stub.lie_range_start
                    lie_total = stub.lie_range_total
                sent_headers = {k.lower(): v for k, v in self.headers.items()}

                if_match = sent_headers.get("if-match")
                if if_match is not None:
                    candidates = [t.strip() for t in if_match.split(",")]
                    if "*" not in candidates and etag not in candidates:
                        stub._record(self, 412)
                        return self._json(
                            412, {"detail": "If-Match failed"}, {"ETag": etag}
                        )

                total = len(body)
                start, end, status = 0, total - 1, 200
                range_header = sent_headers.get("range")
                if range_header and not behavior["ignore_range"]:
                    match = re.fullmatch(r"bytes=(\d+)-(\d*)", range_header)
                    if match is None:
                        stub._record(self, 416)
                        return self._json(
                            416,
                            {"detail": "bad range"},
                            {"Content-Range": f"bytes */{total}"},
                        )
                    start = int(match.group(1))
                    if match.group(2):
                        end = min(int(match.group(2)), total - 1)
                    if start >= total or start > end:
                        stub._record(self, 416)
                        return self._json(
                            416,
                            {"detail": "unsatisfiable"},
                            {"Content-Range": f"bytes */{total}"},
                        )
                    status = 206

                payload = body[start : end + 1]
                self.send_response(status)
                self.send_header("ETag", etag_override or etag)
                self.send_header("Accept-Ranges", "bytes")
                self.send_header("Content-Length", str(len(payload)))
                if status == 206:
                    shown_start = start + lie_start
                    shown_total = lie_total if lie_total is not None else total
                    self.send_header(
                        "Content-Range",
                        f"bytes {shown_start}-{end}/{shown_total}",
                    )
                self.end_headers()
                stub._record(self, status)

                cut = False
                with stub.lock:
                    if stub.behavior["cuts"] > 0:
                        stub.behavior["cuts"] -= 1
                        cut = True
                if cut:
                    n = max(1, int(len(payload) * behavior["cut_fraction"]))
                    self.wfile.write(payload[:n])
                    self.wfile.flush()
                    self.close_connection = True
                    return
                self.wfile.write(payload)

        self._server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        self._server.daemon_threads = True
        self._thread = threading.Thread(target=self._server.serve_forever, daemon=True)
        self._thread.start()

    def _record(self, handler, status):
        with self.lock:
            self.requests.append(
                {
                    "path": handler.path,
                    "status": status,
                    "headers": {k.lower(): v for k, v in handler.headers.items()},
                }
            )

    def _etag(self, name):
        return f'"{hashlib.sha256(self.bodies[name]).hexdigest()[:32]}"'

    def _meta(self, run_id, name):
        digest = hashlib.sha256(self.bodies[name]).hexdigest()
        return {
            "name": name,
            "path": f"/runs/{run_id}/files/{name}.fastq.gz",
            "length": len(self.bodies[name]),
            "etag": f'"{digest[:32]}"',
            "sha256": digest,
        }

    @property
    def url(self):
        return f"http://127.0.0.1:{self._server.server_address[1]}"

    def file_requests(self, name: str | None = None):
        result = [r for r in self.requests if r["path"].endswith(".fastq.gz")]
        if name is not None:
            result = [r for r in result if f"/{name}." in r["path"]]
        return result

    def set_bodies(self, r1: bytes, r2: bytes):
        with self.lock:
            self.bodies = {"R1": r1, "R2": r2}

    def shutdown(self):
        self._server.shutdown()
        self._server.server_close()
        self._thread.join(timeout=5)


@pytest.fixture
def instrument():
    r1, r2 = make_run_data()
    stub = StubInstrument(r1, r2)
    yield stub
    stub.shutdown()


@pytest.fixture
def client(tmp_path, monkeypatch, instrument):
    monkeypatch.setenv("ACQ_DATA_DIR", str(tmp_path / "acq"))
    monkeypatch.setenv("INSTRUMENT_BASE_URL", instrument.url)
    monkeypatch.setattr(main_module, "TMP_ROOT", tmp_path / "demux-tmp")
    with TestClient(main_module.app) as test_client:
        yield test_client


def file_entry(status: dict, name: str) -> dict:
    return next(f for f in status["files"] if f["name"] == name)


def test_import_completes_then_status_and_demux(client, instrument, tmp_path):
    response = client.post("/api/runs/run1/import")
    assert response.status_code == 200
    status = response.json()
    assert status["state"] == "usable"
    for entry in status["files"]:
        assert entry["digest_ok"] is True
        assert entry["bytes_on_disk"] == entry["length"]

    # Query endpoints reflect the same state.
    assert client.get("/api/runs/run1").json()["state"] == "usable"
    runs = client.get("/api/runs").json()["runs"]
    assert [r["run_id"] for r in runs] == ["run1"]

    # Final files are on disk under the run directory.
    run_dir = tmp_path / "acq" / "runs" / "run1"
    assert (run_dir / "R1.fastq.gz").is_file()
    assert (run_dir / "R2.fastq.gz").is_file()
    assert not list(run_dir.glob("*.part"))

    # A usable run demultiplexes through the original logic.
    demux = client.post(
        "/api/runs/run1/demux",
        data={"barcodes": ["ACGTACGT", "TTTTAAAA"]},
    )
    assert demux.status_code == 200
    assert demux.headers["content-type"].startswith("application/zip")
    with zipfile.ZipFile(io.BytesIO(demux.content)) as archive:
        manifest = json.loads(archive.read("manifest.json"))
    counts = {g["name"]: g["pair_count"] for g in manifest["groups"]}
    assert counts == EXPECTED_COUNTS
    assert manifest["total_pairs"] == PAIRS


def test_cut_stream_keeps_prefix_and_resume_finishes(client, instrument, tmp_path):
    instrument.behavior["cuts"] = 1
    response = client.post("/api/runs/run1/import")
    assert response.status_code == 200
    status = response.json()
    assert status["state"] == "incomplete"
    r1 = file_entry(status, "R1")
    assert 0 < r1["bytes_on_disk"] < r1["length"]
    assert file_entry(status, "R2")["bytes_on_disk"] == 0
    assert status["last_error"]["code"] == "transfer_interrupted"

    # Recorded progress never exceeds the bytes actually on disk.
    part = tmp_path / "acq" / "runs" / "run1" / "R1.part"
    assert part.is_file()
    assert r1["bytes_on_disk"] <= part.stat().st_size

    response = client.post("/api/runs/run1/resume")
    assert response.status_code == 200
    assert response.json()["state"] == "usable"

    # The resume used Range + If-Match with the manifest's strong ETag.
    r1_requests = instrument.file_requests("R1")
    assert len(r1_requests) == 2
    resume_request = r1_requests[1]
    assert resume_request["status"] == 206
    assert resume_request["headers"].get("range") == f"bytes={r1['bytes_on_disk']}-"
    assert resume_request["headers"].get("if-match") == r1["etag"]


def test_restart_reconciles_progress_with_disk(tmp_path, monkeypatch, instrument):
    instrument.behavior["cuts"] = 1
    acq_dir = tmp_path / "acq"
    monkeypatch.setenv("ACQ_DATA_DIR", str(acq_dir))
    monkeypatch.setenv("INSTRUMENT_BASE_URL", instrument.url)
    with TestClient(main_module.app) as first:
        status = first.post("/api/runs/run1/import").json()
        assert status["state"] == "incomplete"

    # Simulate a crash that left the DB claiming more than was written.
    part = acq_dir / "runs" / "run1" / "R1.part"
    actual = part.stat().st_size
    db = sqlite3.connect(acq_dir / "acquisition.db")
    db.execute(
        "UPDATE run_files SET bytes_on_disk = ? WHERE run_id = 'run1' AND name = 'R1'",
        (actual + 5000,),
    )
    db.commit()
    db.close()

    # The next "process" clamps progress to what is really on disk.
    with TestClient(main_module.app) as second:
        status = second.get("/api/runs/run1").json()
        assert file_entry(status, "R1")["bytes_on_disk"] == actual
        response = second.post("/api/runs/run1/resume")
        assert response.json()["state"] == "usable"


def test_usable_run_with_deleted_files_falls_back_to_incomplete(
    tmp_path, monkeypatch, instrument
):
    acq_dir = tmp_path / "acq"
    monkeypatch.setenv("ACQ_DATA_DIR", str(acq_dir))
    monkeypatch.setenv("INSTRUMENT_BASE_URL", instrument.url)
    with TestClient(main_module.app) as first:
        assert first.post("/api/runs/run1/import").json()["state"] == "usable"

    for name in ("R1", "R2"):
        (acq_dir / "runs" / "run1" / f"{name}.fastq.gz").unlink()

    with TestClient(main_module.app) as second:
        status = second.get("/api/runs/run1").json()
        assert status["state"] == "incomplete"
        assert all(f["bytes_on_disk"] == 0 for f in status["files"])
        assert second.post("/api/runs/run1/resume").json()["state"] == "usable"


def test_ignored_range_restarts_file_from_zero(client, instrument):
    instrument.behavior["cuts"] = 1
    instrument.behavior["ignore_range"] = True
    status = client.post("/api/runs/run1/import").json()
    assert status["state"] == "incomplete"
    partial = file_entry(status, "R1")["bytes_on_disk"]
    assert partial > 0

    response = client.post("/api/runs/run1/resume")
    assert response.status_code == 200
    assert response.json()["state"] == "usable"

    # The source ignored Range and answered 200; the client must have
    # restarted the file (a silent append would break the SHA-256 digest).
    r1_requests = instrument.file_requests("R1")
    assert r1_requests[1]["headers"].get("range") == f"bytes={partial}-"
    assert r1_requests[1]["status"] == 200
    final = client.get("/api/runs/run1").json()
    assert file_entry(final, "R1")["digest_ok"] is True


def test_version_change_rejects_batch_and_clears_staging(
    client, instrument, tmp_path
):
    instrument.behavior["cuts"] = 1
    assert client.post("/api/runs/run1/import").json()["state"] == "incomplete"
    run_dir = tmp_path / "acq" / "runs" / "run1"
    assert (run_dir / "R1.part").is_file()

    # Source rotates to a new version: new content, new ETags.
    instrument.set_bodies(*make_run_data(seed=999))
    response = client.post("/api/runs/run1/resume")
    assert response.status_code == 409
    assert response.json()["detail"]["code"] == "source_version_changed"
    assert any(r["status"] == 412 for r in instrument.file_requests("R1"))

    status = client.get("/api/runs/run1").json()
    assert status["state"] == "failed"
    assert status["last_error"]["code"] == "source_version_changed"
    assert not run_dir.exists()  # staging cleared

    # A failed run cannot resume, but a fresh import picks up the new version.
    assert client.post("/api/runs/run1/resume").status_code == 409
    assert client.post("/api/runs/run1/import").json()["state"] == "usable"


def test_etag_mismatch_on_206_rejects_batch(client, instrument):
    instrument.behavior["cuts"] = 1
    assert client.post("/api/runs/run1/import").json()["state"] == "incomplete"
    instrument.response_etag_override = '"deadbeefdeadbeef"'
    response = client.post("/api/runs/run1/resume")
    assert response.status_code == 409
    assert response.json()["detail"]["code"] == "source_version_changed"


def test_etag_mismatch_on_fresh_download_rejects_batch(client, instrument):
    instrument.response_etag_override = '"deadbeefdeadbeef"'
    response = client.post("/api/runs/run1/import")
    assert response.status_code == 409
    assert response.json()["detail"]["code"] == "source_version_changed"


def test_bad_content_range_start_is_not_appended(client, instrument):
    instrument.behavior["cuts"] = 1
    status = client.post("/api/runs/run1/import").json()
    before = file_entry(status, "R1")["bytes_on_disk"]

    instrument.lie_range_start = 10
    response = client.post("/api/runs/run1/resume")
    assert response.status_code == 200
    status = response.json()
    assert status["state"] == "incomplete"
    assert status["last_error"]["code"] == "transfer_interrupted"
    # Nothing was appended at the wrong offset.
    assert file_entry(status, "R1")["bytes_on_disk"] == before

    instrument.lie_range_start = 0
    assert client.post("/api/runs/run1/resume").json()["state"] == "usable"


def test_bad_content_range_total_is_rejected(client, instrument):
    instrument.behavior["cuts"] = 1
    assert client.post("/api/runs/run1/import").json()["state"] == "incomplete"
    instrument.lie_range_total = 1  # smaller than any real total
    response = client.post("/api/runs/run1/resume")
    assert response.status_code == 200
    assert response.json()["state"] == "incomplete"
    assert response.json()["last_error"]["code"] == "transfer_interrupted"


def test_digest_mismatch_rejects_batch_and_clears_staging(
    client, instrument, tmp_path
):
    def corrupt_sha(files):
        for entry in files:
            entry["sha256"] = "0" * 64

    instrument.manifest_mutator = corrupt_sha
    response = client.post("/api/runs/run1/import")
    assert response.status_code == 409
    assert response.json()["detail"]["code"] == "digest_mismatch"
    assert client.get("/api/runs/run1").json()["state"] == "failed"
    assert not (tmp_path / "acq" / "runs" / "run1").exists()


def test_invalid_fastq_content_fails_validation(client, instrument, tmp_path):
    instrument.set_bodies(b"this is not gzip at all", b"neither is this")
    response = client.post("/api/runs/run1/import")
    assert response.status_code == 502
    assert response.json()["detail"]["code"] == "source_data_invalid"
    status = client.get("/api/runs/run1").json()
    assert status["state"] == "failed"
    assert status["last_error"]["code"] == "gzip_error"
    assert not (tmp_path / "acq" / "runs" / "run1").exists()


def test_manifest_limits_and_strong_etag_enforced(client, instrument):
    def too_big(files):
        for entry in files:
            entry["length"] = 64 * 1024 * 1024 + 1

    instrument.manifest_mutator = too_big
    response = client.post("/api/runs/run1/import")
    assert response.status_code == 502
    assert response.json()["detail"]["code"] == "invalid_manifest"

    def weak_etag(files):
        for entry in files:
            entry["etag"] = "W/" + entry["etag"]

    instrument.manifest_mutator = weak_etag
    response = client.post("/api/runs/run1/import")
    assert response.status_code == 502
    assert response.json()["detail"]["code"] == "invalid_manifest"

    def foreign_path(files):
        for entry in files:
            entry["path"] = "https://evil.example.com/x"

    instrument.manifest_mutator = foreign_path
    response = client.post("/api/runs/run1/import")
    assert response.status_code == 502
    assert response.json()["detail"]["code"] == "invalid_manifest"

    instrument.manifest_mutator = None
    assert client.post("/api/runs/run1/import").json()["state"] == "usable"


def test_concurrent_acquisition_is_rejected(client):
    acquisition = client.app.state.acquisition
    with acquisition.acquisition_lock("run1"):
        response = client.post("/api/runs/run1/import")
        assert response.status_code == 409
        assert response.json()["detail"]["code"] == "acquisition_in_progress"
        response = client.post("/api/runs/run1/resume")
        assert response.status_code == 409


def test_demux_requires_usable_run(client, instrument):
    assert client.post("/api/runs/never/demux", data={"barcodes": ["ACGTACGT"]}).status_code == 404
    instrument.behavior["cuts"] = 1
    assert client.post("/api/runs/run1/import").json()["state"] == "incomplete"
    response = client.post("/api/runs/run1/demux", data={"barcodes": ["ACGTACGT"]})
    assert response.status_code == 409
    assert response.json()["detail"]["code"] == "run_not_usable"


def test_unknown_runs_and_run_id_validation(client):
    assert client.get("/api/runs/nope").status_code == 404
    assert client.post("/api/runs/nope/resume").status_code == 404
    assert client.post("/api/runs/bad id/import").status_code == 422
    assert client.post("/api/runs/-leading-dash/import").status_code == 422
    assert client.post("/api/runs/bad%24dollar/import").status_code == 422
    detail = client.post("/api/runs/bad id/import").json()["detail"]
    assert detail["code"] == "invalid_run_id"


def test_reimport_of_usable_run_is_idempotent(client, instrument):
    assert client.post("/api/runs/run1/import").json()["state"] == "usable"
    assert client.post("/api/runs/run1/import").json()["state"] == "usable"
    # The second import did not download anything again.
    assert len(instrument.file_requests()) == 2


def test_source_unavailable_gives_502(tmp_path, monkeypatch):
    monkeypatch.setenv("ACQ_DATA_DIR", str(tmp_path / "acq"))
    monkeypatch.setenv("INSTRUMENT_BASE_URL", "http://127.0.0.1:9")
    with TestClient(main_module.app) as test_client:
        response = test_client.post("/api/runs/run1/import")
        assert response.status_code == 502
        assert response.json()["detail"]["code"] == "source_unavailable"
