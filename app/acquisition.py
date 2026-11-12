"""Resumable run acquisition from the fixed instrument export server.

A sequencing-center operator imports a run by ``run_id``.  The source base
URL is fixed by configuration (``INSTRUMENT_BASE_URL``); the source's
manifest pins the two mate files (R1/R2) with their paths, lengths, strong
ETags, and SHA-256 digests (each file at most 64 MiB).  The manifest is
frozen into SQLite at import time together with per-file progress, so a
stream break or a process restart never loses the already-downloaded
prefix: progress recorded in SQLite never exceeds the bytes actually on
disk, and on startup (and before every transfer) the recorded progress is
clamped to the real file sizes.

Resume uses ``Range`` + ``If-Match``.  A ``206`` response is accepted only
if its ETag matches the frozen manifest and its ``Content-Range`` starts at
exactly the current offset with the manifest's total length.  A ``200``
response to a resume request means the source ignored ``Range``: the file
is then truncated and rewritten from byte zero, never appended.  If the
source version changed (``412``, or a different ETag) or a completed file
fails its SHA-256 check, the whole batch is rejected and the staging area
is cleared.

Only when *both* files match their digests *and* pass the existing paired
FASTQ + gzip validation is the run marked ``usable``; before that no
demultiplexing is possible.  A single API process runs at most one
acquisition per run at a time (guarded by an in-process lock).
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import sqlite3
import threading
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator

import httpx

from .demux import DemuxError, validate_paired_fastq

MAX_FILE_BYTES = 64 * 1024 * 1024  # 64 MiB per mate file
MATE_NAMES = ("R1", "R2")

STATE_INCOMPLETE = "incomplete"
STATE_USABLE = "usable"
STATE_FAILED = "failed"

_RUN_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$")
_SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")
_CONTENT_RANGE_RE = re.compile(r"^bytes (\d+)-(\d+)/(\d+|\*)$")

_SCHEMA = """
CREATE TABLE IF NOT EXISTS runs (
    run_id TEXT PRIMARY KEY,
    state TEXT NOT NULL,
    last_error_code TEXT,
    last_error_message TEXT,
    created_at TEXT NOT NULL,
    updated_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS run_files (
    run_id TEXT NOT NULL,
    name TEXT NOT NULL,
    path TEXT NOT NULL,
    length INTEGER NOT NULL,
    etag TEXT NOT NULL,
    sha256 TEXT NOT NULL,
    bytes_on_disk INTEGER NOT NULL DEFAULT 0,
    digest_ok INTEGER NOT NULL DEFAULT 0,
    PRIMARY KEY (run_id, name)
);
"""


def _utcnow() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _is_strong_etag(etag: object) -> bool:
    return isinstance(etag, str) and len(etag) >= 2 and not etag.startswith("W/")


def _etags_equal(left: str, right: str) -> bool:
    """Strong comparison: both tags must be strong and byte-identical."""

    return _is_strong_etag(left) and _is_strong_etag(right) and left == right


class AcquisitionError(Exception):
    """Terminal failure of an import/resume; maps to an HTTP error response."""

    def __init__(
        self,
        code: str,
        message: str,
        status_code: int = 502,
        extra: dict | None = None,
    ):
        super().__init__(message)
        self.code = code
        self.message = message
        self.status_code = status_code
        self.extra = extra or {}

    def as_detail(self) -> dict:
        detail = {"code": self.code, "message": self.message}
        detail.update(self.extra)
        return detail


class TransferInterrupted(Exception):
    """The stream broke or the source misbehaved transiently.

    The downloaded prefix stays on disk and the run stays ``incomplete`` so
    a later resume picks up exactly where this attempt stopped.
    """


class Acquisition:
    """Owns the SQLite manifest/progress store and the transfer engine."""

    def __init__(
        self,
        data_dir: str | Path,
        base_url: str,
        *,
        client: httpx.Client | None = None,
    ):
        self.data_dir = Path(data_dir)
        self.base_url = base_url.rstrip("/")
        self.runs_dir = self.data_dir / "runs"
        self.runs_dir.mkdir(parents=True, exist_ok=True)
        self._db_path = self.data_dir / "acquisition.db"
        self._db = sqlite3.connect(self._db_path, check_same_thread=False)
        self._db.row_factory = sqlite3.Row
        self._db_lock = threading.Lock()
        with self._db_lock, self._db:
            self._db.execute("PRAGMA journal_mode=WAL")
            self._db.executescript(_SCHEMA)
        self._client = client or httpx.Client(
            base_url=self.base_url,
            timeout=httpx.Timeout(connect=5.0, read=15.0, write=15.0, pool=5.0),
            follow_redirects=False,
        )
        self._locks: dict[str, threading.Lock] = {}
        self._locks_guard = threading.Lock()
        self.reconcile()

    def close(self) -> None:
        with self._db_lock:
            self._db.close()
        self._client.close()

    # ------------------------------------------------------------------
    # paths and low-level store helpers
    # ------------------------------------------------------------------

    def run_dir(self, run_id: str) -> Path:
        return self.runs_dir / run_id

    def part_path(self, run_id: str, name: str) -> Path:
        return self.run_dir(run_id) / f"{name}.part"

    def final_path(self, run_id: str, name: str) -> Path:
        return self.run_dir(run_id) / f"{name}.fastq.gz"

    def _query_one(self, sql: str, args: tuple) -> sqlite3.Row | None:
        with self._db_lock:
            return self._db.execute(sql, args).fetchone()

    def _query_all(self, sql: str, args: tuple = ()) -> list[sqlite3.Row]:
        with self._db_lock:
            return self._db.execute(sql, args).fetchall()

    def _execute(self, sql: str, args: tuple) -> None:
        with self._db_lock, self._db:
            self._db.execute(sql, args)

    @staticmethod
    def check_run_id(run_id: str) -> None:
        if not _RUN_ID_RE.fullmatch(run_id):
            raise AcquisitionError(
                "invalid_run_id",
                "run_id must be 1-64 characters of [A-Za-z0-9._-] and start "
                "with a letter or digit.",
                status_code=422,
            )

    # ------------------------------------------------------------------
    # status queries
    # ------------------------------------------------------------------

    def get_run(self, run_id: str) -> dict | None:
        row = self._query_one("SELECT * FROM runs WHERE run_id = ?", (run_id,))
        if row is None:
            return None
        files = self._query_all(
            "SELECT * FROM run_files WHERE run_id = ? ORDER BY name", (run_id,)
        )
        return self._status_from_rows(row, files)

    def list_runs(self) -> list[dict]:
        rows = self._query_all("SELECT * FROM runs ORDER BY created_at, run_id")
        result = []
        for row in rows:
            files = self._query_all(
                "SELECT * FROM run_files WHERE run_id = ? ORDER BY name",
                (row["run_id"],),
            )
            result.append(self._status_from_rows(row, files))
        return result

    @staticmethod
    def _status_from_rows(run: sqlite3.Row, files: list[sqlite3.Row]) -> dict:
        last_error = None
        if run["last_error_code"] is not None:
            last_error = {
                "code": run["last_error_code"],
                "message": run["last_error_message"],
            }
        return {
            "run_id": run["run_id"],
            "state": run["state"],
            "files": [
                {
                    "name": f["name"],
                    "path": f["path"],
                    "length": f["length"],
                    "etag": f["etag"],
                    "sha256": f["sha256"],
                    "bytes_on_disk": f["bytes_on_disk"],
                    "digest_ok": bool(f["digest_ok"]),
                }
                for f in files
            ],
            "last_error": last_error,
            "created_at": run["created_at"],
            "updated_at": run["updated_at"],
        }

    # ------------------------------------------------------------------
    # crash/restart reconciliation: progress never exceeds bytes on disk
    # ------------------------------------------------------------------

    def reconcile(self) -> None:
        """Clamp every recorded progress to what is actually on disk.

        Runs marked ``usable`` whose final files vanished (data directory
        partially wiped) fall back to ``incomplete`` with zeroed progress so
        they can simply be resumed.
        """

        for run in self._query_all("SELECT * FROM runs"):
            run_id = run["run_id"]
            files = self._query_all(
                "SELECT * FROM run_files WHERE run_id = ?", (run_id,)
            )
            if run["state"] == STATE_USABLE:
                finals_exist = all(
                    self.final_path(run_id, f["name"]).is_file() for f in files
                )
                if not finals_exist:
                    with self._db_lock, self._db:
                        self._db.execute(
                            "UPDATE runs SET state = ?, updated_at = ? WHERE run_id = ?",
                            (STATE_INCOMPLETE, _utcnow(), run_id),
                        )
                        self._db.execute(
                            "UPDATE run_files SET bytes_on_disk = 0, digest_ok = 0 "
                            "WHERE run_id = ?",
                            (run_id,),
                        )
                continue
            if run["state"] != STATE_INCOMPLETE:
                continue
            for f in files:
                self._sync_offset(run_id, f["name"], f["length"])

    def _sync_offset(self, run_id: str, name: str, length: int) -> int:
        """Reconcile one file's recorded progress with its ``.part`` file.

        Returns the trustworthy offset: never more than the bytes actually
        on disk, never more than the manifest length.  Any on-disk tail
        beyond the recorded progress (a crash between write and progress
        update) is truncated away so appends always join a known prefix.
        """

        row = self._query_one(
            "SELECT bytes_on_disk, digest_ok FROM run_files WHERE run_id = ? AND name = ?",
            (run_id, name),
        )
        assert row is not None
        part = self.part_path(run_id, name)
        actual = part.stat().st_size if part.exists() else 0
        offset = min(row["bytes_on_disk"], actual, length)
        if actual > offset:
            with open(part, "r+b") as handle:
                handle.truncate(offset)
        digest_ok = row["digest_ok"] and offset == length
        if offset != row["bytes_on_disk"] or digest_ok != row["digest_ok"]:
            self._execute(
                "UPDATE run_files SET bytes_on_disk = ?, digest_ok = ? "
                "WHERE run_id = ? AND name = ?",
                (offset, int(digest_ok), run_id, name),
            )
        return offset

    # ------------------------------------------------------------------
    # run lifecycle mutations
    # ------------------------------------------------------------------

    def _reset_run(self, run_id: str, files: list[dict]) -> None:
        """Freeze a fresh manifest for the run and clear all staging."""

        shutil.rmtree(self.run_dir(run_id), ignore_errors=True)
        self.run_dir(run_id).mkdir(parents=True, exist_ok=True)
        now = _utcnow()
        with self._db_lock, self._db:
            self._db.execute(
                "INSERT INTO runs (run_id, state, created_at, updated_at) "
                "VALUES (?, ?, ?, ?) "
                "ON CONFLICT(run_id) DO UPDATE SET state = excluded.state, "
                "last_error_code = NULL, last_error_message = NULL, "
                "updated_at = excluded.updated_at",
                (run_id, STATE_INCOMPLETE, now, now),
            )
            self._db.execute("DELETE FROM run_files WHERE run_id = ?", (run_id,))
            for f in files:
                self._db.execute(
                    "INSERT INTO run_files "
                    "(run_id, name, path, length, etag, sha256) "
                    "VALUES (?, ?, ?, ?, ?, ?)",
                    (run_id, f["name"], f["path"], f["length"], f["etag"], f["sha256"]),
                )

    def _set_progress(self, run_id: str, name: str, offset: int, length: int) -> None:
        # digest_ok can only stay set while the file is complete on disk.
        self._execute(
            "UPDATE run_files SET bytes_on_disk = ?, "
            "digest_ok = CASE WHEN ? >= ? THEN digest_ok ELSE 0 END "
            "WHERE run_id = ? AND name = ?",
            (offset, offset, length, run_id, name),
        )

    def _set_digest_ok(self, run_id: str, name: str) -> None:
        self._execute(
            "UPDATE run_files SET digest_ok = 1 WHERE run_id = ? AND name = ?",
            (run_id, name),
        )

    def _set_last_error(self, run_id: str, code: str, message: str) -> None:
        self._execute(
            "UPDATE runs SET last_error_code = ?, last_error_message = ?, "
            "updated_at = ? WHERE run_id = ?",
            (code, message, _utcnow(), run_id),
        )

    def _fail_run(self, run_id: str, code: str, message: str) -> None:
        """Reject the whole batch: mark failed and clear the staging area."""

        self._execute(
            "UPDATE runs SET state = ?, last_error_code = ?, "
            "last_error_message = ?, updated_at = ? WHERE run_id = ?",
            (STATE_FAILED, code, message, _utcnow(), run_id),
        )
        shutil.rmtree(self.run_dir(run_id), ignore_errors=True)

    def _mark_usable(self, run_id: str) -> None:
        self._execute(
            "UPDATE runs SET state = ?, last_error_code = NULL, "
            "last_error_message = NULL, updated_at = ? WHERE run_id = ?",
            (STATE_USABLE, _utcnow(), run_id),
        )

    # ------------------------------------------------------------------
    # per-run mutual exclusion (single process, one acquisition per run)
    # ------------------------------------------------------------------

    @contextmanager
    def acquisition_lock(self, run_id: str) -> Iterator[None]:
        with self._locks_guard:
            lock = self._locks.setdefault(run_id, threading.Lock())
        if not lock.acquire(blocking=False):
            raise AcquisitionError(
                "acquisition_in_progress",
                f"An acquisition for run {run_id!r} is already in progress; "
                "query its status and retry when it finishes.",
                status_code=409,
            )
        try:
            yield
        finally:
            lock.release()

    # ------------------------------------------------------------------
    # import / resume entry points
    # ------------------------------------------------------------------

    def import_run(self, run_id: str) -> dict:
        """Fetch and freeze the manifest, then transfer until done/interrupted."""

        self.check_run_id(run_id)
        manifest = self._fetch_manifest(run_id)
        files = self._validate_manifest(run_id, manifest)

        existing = self.get_run(run_id)
        if existing is not None and existing["state"] == STATE_USABLE:
            # The run already holds fully validated data; importing again is
            # an idempotent no-op and never discards good data.
            return existing
        if (
            existing is not None
            and existing["state"] == STATE_INCOMPLETE
            and self._manifest_matches(existing, files)
        ):
            # Same manifest as the frozen one: keep the downloaded prefix
            # and continue (equivalent to a resume).
            pass
        else:
            # First import, a rejected batch, or the source now serves a
            # different version: drop whatever was staged and start over.
            self._reset_run(run_id, files)
        return self._transfer_and_finalize(run_id)

    def resume_run(self, run_id: str) -> dict:
        self.check_run_id(run_id)
        run = self.get_run(run_id)
        if run is None:
            raise AcquisitionError(
                "run_not_found",
                f"Run {run_id!r} was never imported.",
                status_code=404,
            )
        if run["state"] == STATE_USABLE:
            return run
        if run["state"] == STATE_FAILED:
            raise AcquisitionError(
                "run_failed",
                f"Run {run_id!r} was rejected "
                f"({run['last_error']['code']}); re-import it to start over.",
                status_code=409,
            )
        return self._transfer_and_finalize(run_id)

    @staticmethod
    def _manifest_matches(run: dict, files: list[dict]) -> bool:
        frozen = {f["name"]: f for f in run["files"]}
        if set(frozen) != {f["name"] for f in files}:
            return False
        return all(
            frozen[f["name"]]["path"] == f["path"]
            and frozen[f["name"]]["length"] == f["length"]
            and frozen[f["name"]]["etag"] == f["etag"]
            and frozen[f["name"]]["sha256"] == f["sha256"]
            for f in files
        )

    # ------------------------------------------------------------------
    # manifest handling
    # ------------------------------------------------------------------

    def _fetch_manifest(self, run_id: str) -> dict:
        try:
            response = self._client.get(f"/runs/{run_id}/manifest")
        except httpx.HTTPError as exc:
            raise AcquisitionError(
                "source_unavailable",
                f"Could not reach the instrument source: {exc}",
                status_code=502,
            ) from exc
        if response.status_code != 200:
            raise AcquisitionError(
                "source_unavailable",
                f"Manifest request for run {run_id!r} returned HTTP "
                f"{response.status_code}.",
                status_code=502,
            )
        try:
            doc = response.json()
        except json.JSONDecodeError as exc:
            raise AcquisitionError(
                "invalid_manifest",
                "The source manifest is not valid JSON.",
                status_code=502,
            ) from exc
        return doc

    @staticmethod
    def _validate_manifest(run_id: str, doc: object) -> list[dict]:
        def fail(reason: str) -> AcquisitionError:
            return AcquisitionError(
                "invalid_manifest",
                f"Invalid manifest for run {run_id!r}: {reason}",
                status_code=502,
            )

        if not isinstance(doc, dict):
            raise fail("top level must be an object.")
        if doc.get("run_id") != run_id:
            raise fail(f"manifest run_id {doc.get('run_id')!r} does not match.")
        files = doc.get("files")
        if not isinstance(files, list) or len(files) != len(MATE_NAMES):
            raise fail("'files' must list exactly the R1 and R2 entries.")

        validated: list[dict] = []
        seen: set[str] = set()
        for entry in files:
            if not isinstance(entry, dict):
                raise fail("each file entry must be an object.")
            name = entry.get("name")
            if name not in MATE_NAMES or name in seen:
                raise fail("file entries must be exactly one 'R1' and one 'R2'.")
            seen.add(name)

            path = entry.get("path")
            if (
                not isinstance(path, str)
                or not path.startswith("/")
                or path.startswith("//")
                or "://" in path
                or ".." in path.split("/")
            ):
                raise fail(f"{name}: 'path' must be an absolute path on the "
                           "configured source.")

            length = entry.get("length")
            if isinstance(length, bool) or not isinstance(length, int):
                raise fail(f"{name}: 'length' must be an integer.")
            if length < 0 or length > MAX_FILE_BYTES:
                raise fail(
                    f"{name}: length {length} exceeds the 64 MiB per-file limit."
                    if length > MAX_FILE_BYTES
                    else f"{name}: length must not be negative."
                )

            etag = entry.get("etag")
            if not _is_strong_etag(etag):
                raise fail(f"{name}: 'etag' must be a strong entity tag.")

            sha256 = entry.get("sha256")
            if not isinstance(sha256, str) or not _SHA256_RE.fullmatch(sha256):
                raise fail(f"{name}: 'sha256' must be 64 hex characters.")

            validated.append(
                {
                    "name": name,
                    "path": path,
                    "length": length,
                    "etag": etag,
                    "sha256": sha256.lower(),
                }
            )
        validated.sort(key=lambda f: MATE_NAMES.index(f["name"]))
        return validated

    # ------------------------------------------------------------------
    # transfer engine
    # ------------------------------------------------------------------

    def _transfer_and_finalize(self, run_id: str) -> dict:
        run = self.get_run(run_id)
        assert run is not None

        try:
            for file_entry in run["files"]:
                self._download_file(run_id, file_entry)
        except TransferInterrupted as exc:
            self._set_last_error(run_id, "transfer_interrupted", str(exc))
            return self.get_run(run_id)

        # Both files are fully on disk: verify the pinned SHA-256 digests.
        for file_entry in run["files"]:
            if file_entry["digest_ok"]:
                continue
            digest = self._sha256_file(self.part_path(run_id, file_entry["name"]))
            if digest != file_entry["sha256"]:
                message = (
                    f"{file_entry['name']} SHA-256 is {digest}, but the frozen "
                    f"manifest pins {file_entry['sha256']}; the batch is rejected."
                )
                self._fail_run(run_id, "digest_mismatch", message)
                raise AcquisitionError("digest_mismatch", message, status_code=409)
            self._set_digest_ok(run_id, file_entry["name"])

        # Digests match: run the existing paired FASTQ + gzip validation.
        try:
            with open(self.part_path(run_id, "R1"), "rb") as r1_handle, open(
                self.part_path(run_id, "R2"), "rb"
            ) as r2_handle:
                validate_paired_fastq(r1_handle, r2_handle)
        except DemuxError as exc:
            message = (
                f"Downloaded files failed FASTQ validation: {exc.message}"
            )
            self._fail_run(run_id, exc.code, message)
            raise AcquisitionError(
                "source_data_invalid",
                message,
                status_code=502,
                extra={"validation": exc.as_detail()},
            ) from exc

        # Everything passed: promote both parts to final files and flip the
        # run to usable in one step, so it never becomes usable half-way.
        for name in MATE_NAMES:
            os.replace(self.part_path(run_id, name), self.final_path(run_id, name))
        self._mark_usable(run_id)
        return self.get_run(run_id)

    def _download_file(self, run_id: str, file_entry: dict) -> None:
        name = file_entry["name"]
        length = file_entry["length"]
        etag = file_entry["etag"]
        self.run_dir(run_id).mkdir(parents=True, exist_ok=True)

        while True:
            offset = self._sync_offset(run_id, name, length)
            if offset == length:
                return

            headers = {}
            if offset > 0:
                headers["Range"] = f"bytes={offset}-"
                headers["If-Match"] = etag

            try:
                with self._client.stream(
                    "GET", file_entry["path"], headers=headers
                ) as response:
                    if response.status_code == 412:
                        self._reject_version_changed(
                            run_id,
                            f"The source rejected If-Match {etag} for {name} "
                            "(HTTP 412): the served version changed.",
                        )
                    if response.status_code == 200:
                        # Fresh download, or the source ignored our Range
                        # header: either way the body starts at byte 0 and
                        # the file is rewritten from scratch, never appended.
                        self._check_response_etag(run_id, name, response, etag)
                        write_offset = 0
                    elif response.status_code == 206:
                        self._check_response_etag(run_id, name, response, etag)
                        write_offset = self._check_content_range(
                            response, offset, length
                        )
                    else:
                        raise TransferInterrupted(
                            f"GET {file_entry['path']} returned unexpected "
                            f"HTTP {response.status_code}."
                        )
                    self._consume(run_id, name, response, write_offset, length)
            except (TransferInterrupted, AcquisitionError):
                raise
            except httpx.HTTPError as exc:
                raise TransferInterrupted(
                    f"Connection to the source broke while fetching {name}: {exc}"
                ) from exc
            except OSError as exc:
                raise TransferInterrupted(
                    f"Local disk error while storing {name}: {exc}"
                ) from exc

    def _check_response_etag(
        self, run_id: str, name: str, response: httpx.Response, etag: str
    ) -> None:
        served = response.headers.get("etag")
        if served is None:
            raise TransferInterrupted(
                f"The source sent no ETag for {name}; cannot pin the version."
            )
        if not _etags_equal(served, etag):
            self._reject_version_changed(
                run_id,
                f"The source now serves ETag {served} for {name}, but the "
                f"frozen manifest pins {etag}: the served version changed.",
            )

    def _reject_version_changed(self, run_id: str, message: str) -> None:
        self._fail_run(run_id, "source_version_changed", message)
        raise AcquisitionError(
            "source_version_changed", message, status_code=409
        )

    @staticmethod
    def _check_content_range(
        response: httpx.Response, offset: int, length: int
    ) -> int:
        header = response.headers.get("content-range")
        if header is None:
            raise TransferInterrupted("206 response carries no Content-Range.")
        match = _CONTENT_RANGE_RE.match(header.strip())
        if match is None:
            raise TransferInterrupted(
                f"Cannot parse Content-Range {header!r}."
            )
        start = int(match.group(1))
        end = int(match.group(2))
        total = match.group(3)
        if total == "*" or int(total) != length:
            raise TransferInterrupted(
                f"Content-Range total {total!r} differs from the manifest "
                f"length {length}."
            )
        if start != offset:
            raise TransferInterrupted(
                f"Content-Range starts at byte {start}, but the local prefix "
                f"is {offset} bytes; refusing to append at the wrong offset."
            )
        if end < start or end >= length:
            raise TransferInterrupted(
                f"Content-Range end {end} is outside the manifest length."
            )
        return start

    def _consume(
        self,
        run_id: str,
        name: str,
        response: httpx.Response,
        offset: int,
        length: int,
    ) -> None:
        """Write the response body at ``offset``, checkpointing per chunk.

        Progress is recorded in SQLite only after the bytes are flushed to
        the file, so the stored progress never exceeds what is on disk.
        """

        part = self.part_path(run_id, name)
        written = offset
        try:
            with open(part, "r+b" if part.exists() else "wb") as handle:
                handle.seek(offset)
                handle.truncate()
                # No chunking of our own: transport chunks are written as
                # they arrive so a broken stream leaves the longest
                # possible prefix on disk.
                for chunk in response.iter_bytes():
                    if written + len(chunk) > length:
                        raise TransferInterrupted(
                            f"The source sent more than the manifest length "
                            f"({length} bytes) for {name}."
                        )
                    handle.write(chunk)
                    written += len(chunk)
                    handle.flush()
                    self._set_progress(run_id, name, written, length)
                handle.flush()
                os.fsync(handle.fileno())
        except httpx.HTTPError as exc:
            raise TransferInterrupted(
                f"Stream broke after {written} bytes of {name}: {exc}"
            ) from exc

    @staticmethod
    def _sha256_file(path: Path) -> str:
        digest = hashlib.sha256()
        with open(path, "rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        return digest.hexdigest()
