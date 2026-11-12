from __future__ import annotations

import logging
import os
from contextlib import asynccontextmanager
from datetime import datetime, timezone
from pathlib import Path

import anyio.to_thread
from fastapi import FastAPI, File, Form, Request, UploadFile
from fastapi.responses import JSONResponse, RedirectResponse, StreamingResponse
from pydantic import BaseModel
from starlette.background import BackgroundTask

from .acquisition import (
    STATE_USABLE,
    Acquisition,
    AcquisitionError,
)
from .demux import DemuxError, DemuxPackage, normalize_barcodes, prepare_package

logger = logging.getLogger("demux")
logging.basicConfig(level=os.getenv("LOG_LEVEL", "INFO"))

TMP_ROOT = Path(os.getenv("DEMUX_TMP_ROOT", "/tmp/demux"))
STREAM_CHUNK = 128 * 1024

DESCRIPTION = """
按 8 位 A/C/G/T 条码把一个双端 FASTQ 混样拆分为样本包。

### 从仪器源采集 run（断网续传）

测序中心数据员可脱离仪器网络，按 `run_id` 从配置固定的仪器导出服务器
（`INSTRUMENT_BASE_URL`）采集 R1/R2：

1. `POST /api/runs/{run_id}/import`：拉取源清单（两份文件的路径、长度、
   强 ETag、SHA-256，每份 ≤ 64 MiB），冻结进 SQLite 后开始分块落盘。
2. 断流或进程重启后，`POST /api/runs/{run_id}/resume` 用
   `Range` + `If-Match` 续传；进度只记录实际落盘字节。源若忽略
   `Range` 返回 200，则该文件从头覆盖重下；源换版（412 或 ETag 变化）
   或 SHA-256 不符会拒绝整批并清除暂存。
3. `GET /api/runs/{run_id}` 随时查询状态和每份文件的落盘进度。
4. 两份文件都通过摘要及 FASTQ 配对、gzip 校验后，run 才整体标为
   `usable`；此前不能拆样。`POST /api/runs/{run_id}/demux`
   对可用 run 执行与上传接口相同的拆样逻辑并下载样本 ZIP。

服务以单 API 进程运行，同一 run 同一时间只允许一个采集；重复导入/续传
请求会收到 409。

### 直接上传拆样

* `r1`、`r2`：两份 **gzip** 压缩的 FASTQ，必须都是严格的四行一条记录。
* `barcodes`：1 到 16 个唯一的 8 位 A/C/G/T 样本条码。Swagger UI 中可点
  “Add item” 添加多个字段，也可以在一个字段里用逗号/换行分隔。
* R1 头必须为 `@ID/1 barcode=<观测条码>`，R2 头必须为 `@ID/2`；两条记录的
  ID 必须相同。
* 读段长度不超过 1000，序列与质量字符串长度必须相同。

### 条码判定

观测条码与每个样本条码逐位比较。**仅当恰好一个样本条码的 Hamming
距离 ≤ 1 时**，该读段对归入这个样本；距离为 0 的精确匹配也遵守这一
规则，所以如果精确匹配的同时又落入另一个条码的 1 位范围内，也会进入
`unassigned`。没有候选或候选多于一个也进入 `unassigned`。一对 R1/R2 永远
一起移动，各组内保持原始顺序和原始四行文本。

### 输出

完整读取两路、验证 gzip 尾部 CRC/长度后，才返回 ZIP。ZIP 包含
`manifest.json`、`pair_counts.tsv`，以及每个组（所有提交条码和
`unassigned`）的 R1/R2 gzip。空组也有合法的空 gzip 文件，便于下游
发现样本。

### 错误

坏记录、R1/R2 ID 错配、条数不相等、压缩损坏均返回 400 JSON，并在
`detail.file`、`detail.record`、`detail.line` 中给出位置；请求不会留下
半成品 ZIP。临时文件在失败、下载完成或客户端断开时删除。

### 示例文件

仓库的 `examples/` 目录包含 `normal`、`ambiguous`、`missing_pair` 三组
可直接上传的 gzip 文件；正常样例条码为 `ACGTACGT` 和 `TTTTAAAA`。
"""


@asynccontextmanager
async def lifespan(app: FastAPI):
    # Single-process deployment: the acquisition store is opened once and
    # reconciled with the files on disk before any request is served.
    acquisition = Acquisition(
        data_dir=Path(os.getenv("ACQ_DATA_DIR", "data")),
        base_url=os.getenv("INSTRUMENT_BASE_URL", "http://localhost:9000"),
    )
    app.state.acquisition = acquisition
    yield
    acquisition.close()


app = FastAPI(
    title="Paired-end FASTQ barcode demultiplexer",
    version="1.1.0",
    description=DESCRIPTION,
    docs_url="/docs",
    redoc_url=None,
    lifespan=lifespan,
)


class ErrorDetail(BaseModel):
    code: str
    message: str
    file: str | None = None
    record: int | None = None
    line: int | None = None
    locations: dict[str, dict[str, int]] | None = None


class ErrorResponse(BaseModel):
    detail: ErrorDetail


class RunError(BaseModel):
    code: str
    message: str


class AcquiredFile(BaseModel):
    name: str
    path: str
    length: int
    etag: str
    sha256: str
    bytes_on_disk: int
    digest_ok: bool


class RunStatus(BaseModel):
    run_id: str
    state: str
    files: list[AcquiredFile]
    last_error: RunError | None
    created_at: str
    updated_at: str


class RunList(BaseModel):
    runs: list[RunStatus]


@app.exception_handler(DemuxError)
async def demux_error_handler(_request, exc: DemuxError) -> JSONResponse:
    return JSONResponse(status_code=exc.status_code, content={"detail": exc.as_detail()})


@app.exception_handler(AcquisitionError)
async def acquisition_error_handler(_request, exc: AcquisitionError) -> JSONResponse:
    return JSONResponse(status_code=exc.status_code, content={"detail": exc.as_detail()})


@app.get("/", include_in_schema=False)
def root() -> RedirectResponse:
    return RedirectResponse(url="/docs")


@app.get("/health", summary="Health check")
def health() -> dict[str, str]:
    return {"status": "ok"}


def _acquisition(request: Request) -> Acquisition:
    return request.app.state.acquisition


async def _stream_zip(package: DemuxPackage):
    """Yield the already-complete ZIP, deleting request files in all exits."""

    try:
        with open(package.zip_path, "rb") as handle:
            while True:
                chunk = await anyio.to_thread.run_sync(handle.read, STREAM_CHUNK)
                if not chunk:
                    break
                yield chunk
    finally:
        # Runs after a normal response, a write error, or client disconnect
        # (Starlette closes the async iterator in those cases).
        package.cleanup()


def _zip_response(package: DemuxPackage) -> StreamingResponse:
    filename = "demux-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + ".zip"
    headers = {
        "Content-Disposition": f'attachment; filename="{filename}"',
        "Content-Length": str(package.zip_path.stat().st_size),
        "Cache-Control": "no-store",
        "X-Content-Type-Options": "nosniff",
        "X-Demux-Total-Pairs": str(package.manifest["total_pairs"]),
    }
    response = StreamingResponse(
        _stream_zip(package),
        media_type="application/zip",
        headers=headers,
    )
    # Belt and braces: the generator's finally is the primary cleanup path,
    # and this background task makes cleanup happen even for unusual ASGI
    # response paths. cleanup() is idempotent.
    response.background = BackgroundTask(package.cleanup)
    return response


def _build_package(holder: dict[str, DemuxPackage]) -> DemuxPackage:
    package = holder.get("package")
    if package is None:  # defensive; prepare_package either returns or raises
        raise DemuxError("internal_error", "Package was not created.", status_code=500)
    return package


@app.post(
    "/api/demux",
    summary="Upload R1/R2 and download a demultiplexed ZIP",
    description=(
        "Upload two gzip FASTQ files and sample barcodes. The response is only "
        "produced after both streams have been fully read and validated. On "
        "success the body is `application/zip`; on a record, pairing, or gzip "
        "error the body is a JSON error with a file/record/line location."
    ),
    response_class=StreamingResponse,
    responses={
        200: {
            "description": "Complete ZIP package",
            "content": {
                "application/zip": {
                    "schema": {"type": "string", "format": "binary"}
                }
            },
        },
        400: {"model": ErrorResponse, "description": "FASTQ or gzip error"},
        422: {"model": ErrorResponse, "description": "Invalid request"},
        500: {"model": ErrorResponse, "description": "Internal error"},
    },
)
async def upload_and_demux(
    r1: UploadFile = File(..., description="gzip-compressed R1 FASTQ"),
    r2: UploadFile = File(..., description="gzip-compressed R2 FASTQ"),
    barcodes: list[str] = Form(
        ...,
        description="1-16 unique 8-base A/C/G/T barcodes (repeat field or separate with commas/newlines)",
    ),
) -> StreamingResponse:
    holder: dict[str, DemuxPackage] = {}

    def process(barcodes_for_run: list[str]) -> None:
        # UploadFile.file is a spooled temporary file; seek before reading.
        r1.file.seek(0)
        r2.file.seek(0)
        holder["package"] = prepare_package(
            r1.file,
            r2.file,
            barcodes_for_run,
            tmp_root=TMP_ROOT,
        )

    try:
        normalized = normalize_barcodes(barcodes)
        await anyio.to_thread.run_sync(process, normalized)
    except DemuxError:
        raise
    except Exception as exc:  # pragma: no cover - exercised by disk failures
        logger.exception("Unexpected demultiplexing failure")
        raise DemuxError(
            "internal_error",
            "An internal error occurred while building the package.",
            status_code=500,
        ) from exc
    except BaseException:
        # Request cancellation can arrive while the worker finishes. Reclaim
        # the package even though its future response will never be sent.
        package = holder.get("package")
        if package is not None:
            package.cleanup()
        raise
    finally:
        # Both source spool files are no longer needed: the ZIP is complete.
        try:
            r1.file.close()
            r2.file.close()
        except Exception:
            logger.warning("Could not close upload spool files", exc_info=True)

    return _zip_response(_build_package(holder))


# ----------------------------------------------------------------------
# Run acquisition from the configured instrument source
# ----------------------------------------------------------------------

RUN_RESPONSES = {
    200: {"model": RunStatus, "description": "Run status after the operation"},
    404: {"model": ErrorResponse, "description": "Unknown run"},
    409: {
        "model": ErrorResponse,
        "description": "Acquisition already running, batch rejected "
        "(version change / digest mismatch), or run not usable",
    },
    422: {"model": ErrorResponse, "description": "Invalid run_id"},
    502: {"model": ErrorResponse, "description": "Source or source-data error"},
}


@app.post(
    "/api/runs/{run_id}/import",
    summary="Import a run from the instrument source and start/resume the transfer",
    description=(
        "Fetch the source manifest (two files with path, length, strong ETag "
        "and SHA-256, each at most 64 MiB), freeze it into SQLite, and "
        "download both files in chunks. If the stream breaks, the downloaded "
        "prefix is kept and the run stays `incomplete`; call `resume` to "
        "continue. If the frozen manifest still matches the source, an "
        "interrupted run keeps its progress; a changed source version or a "
        "digest mismatch rejects the whole batch and clears the staging area."
    ),
    response_model=RunStatus,
    responses=RUN_RESPONSES,
)
async def import_run(run_id: str, request: Request) -> dict:
    acquisition = _acquisition(request)

    def op() -> dict:
        with acquisition.acquisition_lock(run_id):
            return acquisition.import_run(run_id)

    return await anyio.to_thread.run_sync(op)


@app.post(
    "/api/runs/{run_id}/resume",
    summary="Resume an interrupted run transfer (Range + If-Match)",
    description=(
        "Continue downloading from the bytes already on disk. Resume sends "
        "`Range` and `If-Match`; a 206 is accepted only if its ETag and "
        "Content-Range match the frozen manifest, while a 200 (source "
        "ignored Range) restarts that file from byte zero. When both files "
        "are complete their SHA-256 digests and the paired FASTQ/gzip "
        "validation must pass before the run is marked `usable`."
    ),
    response_model=RunStatus,
    responses=RUN_RESPONSES,
)
async def resume_run(run_id: str, request: Request) -> dict:
    acquisition = _acquisition(request)

    def op() -> dict:
        with acquisition.acquisition_lock(run_id):
            return acquisition.resume_run(run_id)

    return await anyio.to_thread.run_sync(op)


@app.get(
    "/api/runs/{run_id}",
    summary="Query run acquisition status and per-file progress",
    response_model=RunStatus,
    responses={404: {"model": ErrorResponse, "description": "Unknown run"}},
)
async def get_run(run_id: str, request: Request) -> dict:
    run = await anyio.to_thread.run_sync(_acquisition(request).get_run, run_id)
    if run is None:
        raise AcquisitionError(
            "run_not_found", f"Run {run_id!r} was never imported.", status_code=404
        )
    return run


@app.get(
    "/api/runs",
    summary="List all known runs",
    response_model=RunList,
)
async def list_runs(request: Request) -> dict:
    runs = await anyio.to_thread.run_sync(_acquisition(request).list_runs)
    return {"runs": runs}


@app.post(
    "/api/runs/{run_id}/demux",
    summary="Demultiplex a usable run and download the sample ZIP",
    description=(
        "Only runs marked `usable` (both files passed SHA-256 and paired "
        "FASTQ/gzip validation) can be demultiplexed. Runs the same "
        "splitting logic as the upload endpoint and returns the same ZIP."
    ),
    response_class=StreamingResponse,
    responses={
        200: {
            "description": "Complete ZIP package",
            "content": {
                "application/zip": {
                    "schema": {"type": "string", "format": "binary"}
                }
            },
        },
        404: {"model": ErrorResponse, "description": "Unknown run"},
        409: {"model": ErrorResponse, "description": "Run is not usable yet"},
        400: {"model": ErrorResponse, "description": "FASTQ or gzip error"},
        422: {"model": ErrorResponse, "description": "Invalid request"},
    },
)
async def demux_acquired_run(
    run_id: str,
    request: Request,
    barcodes: list[str] = Form(
        ...,
        description="1-16 unique 8-base A/C/G/T barcodes (repeat field or separate with commas/newlines)",
    ),
) -> StreamingResponse:
    acquisition = _acquisition(request)
    run = await anyio.to_thread.run_sync(acquisition.get_run, run_id)
    if run is None:
        raise AcquisitionError(
            "run_not_found", f"Run {run_id!r} was never imported.", status_code=404
        )
    if run["state"] != STATE_USABLE:
        raise AcquisitionError(
            "run_not_usable",
            f"Run {run_id!r} is {run['state']}; both files must pass digest "
            "and FASTQ validation before demultiplexing.",
            status_code=409,
        )

    r1_path = acquisition.final_path(run_id, "R1")
    r2_path = acquisition.final_path(run_id, "R2")
    holder: dict[str, DemuxPackage] = {}

    def process(barcodes_for_run: list[str]) -> None:
        with open(r1_path, "rb") as r1_handle, open(r2_path, "rb") as r2_handle:
            holder["package"] = prepare_package(
                r1_handle,
                r2_handle,
                barcodes_for_run,
                tmp_root=TMP_ROOT,
            )

    try:
        normalized = normalize_barcodes(barcodes)
        await anyio.to_thread.run_sync(process, normalized)
    except DemuxError:
        raise
    except Exception as exc:  # pragma: no cover - exercised by disk failures
        logger.exception("Unexpected demultiplexing failure")
        raise DemuxError(
            "internal_error",
            "An internal error occurred while building the package.",
            status_code=500,
        ) from exc
    except BaseException:
        package = holder.get("package")
        if package is not None:
            package.cleanup()
        raise

    return _zip_response(_build_package(holder))
