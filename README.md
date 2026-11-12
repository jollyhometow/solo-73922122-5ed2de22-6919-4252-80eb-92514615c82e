# 双端 FASTQ 条码拆样服务

一个用 FastAPI 实现的服务，支持两条输入路径：

1. **从仪器源采集 run**：测序中心数据员按 `run_id` 从配置固定的仪器导出
   服务器拉取 R1/R2，断网、断流、进程重启后都能续传；两份文件都通过
   SHA-256 摘要和 FASTQ 配对、gzip 校验后才标为可用，可用后即可拆样
   下载样本 ZIP。
2. **直接上传拆样**：接收 R1、R2 两份 gzip FASTQ 和最多 16 个唯一的
   8 位 A/C/G/T 样本条码，流式拆分为样本包 ZIP。

## 快速启动

```bash
docker compose up --build -d
```

这会启动两个服务：

| 服务 | 端口 | 说明 |
|---|---|---|
| `demux` | 8000 | 拆样 API，文档在 <http://localhost:8000/docs> |
| `instrument` | 9000 | 本地仪器导出源，支持断流、换版、忽略 Range 三种故障注入 |

`demux` 的采集数据（SQLite 清单/进度 + 分块暂存文件）保存在命名卷
`demux-data`（容器内 `/data`），容器重启后进度保留。`/health` 用于健康
检查。停止并清理：

```bash
docker compose down        # 保留数据卷
docker compose down -v     # 连同采集数据一起删除
```

也可以用本地 Python 启动（Python 3.11+）：

```bash
python3 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
python3 instrument/server.py &                 # 仪器源，监听 9000
ACQ_DATA_DIR=./data INSTRUMENT_BASE_URL=http://localhost:9000 \
  uvicorn app.main:app --host 0.0.0.0 --port 8000
```

## 在 /docs 演示采集、续传和拆样

仪器源内置了一个 `demo` run：前两次文件响应会在 40% 处断流，之后恢复
正常，因此无需任何额外操作就能演示完整流程。

1. **导入（被断流）**：打开 <http://localhost:8000/docs>，展开
   `POST /api/runs/{run_id}/import`，Try it out，`run_id` 填 `demo`，
   执行。响应 200，`state` 为 `incomplete`，`files[0].bytes_on_disk`
   约为 `length` 的 40%，`last_error` 是 `transfer_interrupted`。
2. **查询**：用 `GET /api/runs/{run_id}` 查看进度；也可以此时
   `docker compose restart demux` 再查询——进度不变（重启后从落盘字节
   恢复，SQLite 里的进度绝不超过实际落盘字节）。
3. **续传（再次被断流）**：`POST /api/runs/{run_id}/resume`，
   `run_id=demo`。R1 进度前进到约 64%，仍为 `incomplete`。续传请求
   带 `Range: bytes=<已落盘>-` 和 `If-Match: <清单ETag>`。
4. **续传（完成）**：再次 `resume`。两份文件下完后依次校验 SHA-256、
   跑 FASTQ 配对 + gzip 校验，全部通过，`state` 变为 `usable`，
   两份 `.part` 晋级为 `R1.fastq.gz` / `R2.fastq.gz`。
5. **拆样**：`POST /api/runs/{run_id}/demux`，`run_id=demo`，
   `barcodes` 填 `ACGTACGT` 和 `TTTTAAAA`，下载样本 ZIP。

### 演示换版（源端版本变化）

```bash
curl -X POST -H 'Content-Type: application/json' \
  -d '{"cuts": 1}' http://localhost:9000/admin/runs/lab/behavior
# /docs 中 import run_id=lab -> incomplete（断流）
curl -X POST http://localhost:9000/admin/rotate     # 仪器源换版
# /docs 中 resume run_id=lab -> 409 source_version_changed
```

`If-Match` 不再匹配，源返回 412，API 拒绝整批、清除暂存，run 标为
`failed`；再次 `import` 会按新清单重新采集。

### 演示忽略 Range

```bash
curl -X POST -H 'Content-Type: application/json' \
  -d '{"cuts": 1, "ignore_range": true}' \
  http://localhost:9000/admin/runs/ignr/behavior
# /docs 中 import run_id=ignr -> incomplete
# /docs 中 resume run_id=ignr -> 源无视 Range 返回 200，
#   API 将该文件截断从头覆盖（不追加），最终 usable
```

### 演示摘要不符

```bash
curl -X POST -H 'Content-Type: application/json' \
  -d '{"corrupt": true}' http://localhost:9000/admin/runs/evil/behavior
# /docs 中 import run_id=evil -> 409 digest_mismatch，暂存被清除
```

仪器源管理接口（直接打 9000 端口）：

| 接口 | 说明 |
|---|---|
| `GET /runs/{run_id}/manifest` | 清单：两文件路径、长度、强 ETag、SHA-256 |
| `GET /runs/{run_id}/files/{R1\|R2}.fastq.gz` | 文件下载，支持 `Range`/`If-Match` |
| `POST /admin/rotate` | 全局版本 +1，所有文件的 ETag/SHA-256/长度改变 |
| `POST /admin/runs/{run_id}/behavior` | 设置故障：`{"cuts": n, "cut_fraction": f, "ignore_range": b, "corrupt": b}` |
| `GET /admin/runs/{run_id}/behavior` | 查看该 run 的故障配置 |
| `POST /admin/reset` | 版本归零、清空全部故障配置 |
| `GET /admin/state` | 查看当前版本与全部故障配置 |

## 采集 API

所有采集接口都在 `INSTRUMENT_BASE_URL` 固定的源上操作；服务以**单 API
进程**运行，同一 run 同一时间只允许一个采集，并发请求收到
`409 acquisition_in_progress`。

### `POST /api/runs/{run_id}/import`

拉取源清单并冻结进 SQLite（两文件路径、长度、强 ETag、SHA-256，每份
≤ 64 MiB，缺一项即 `502 invalid_manifest`），然后分块落盘到
`runs/{run_id}/{R1,R2}.part`。结果：

* 两份文件下完且校验通过 → `state=usable`；
* 中途断流/源异常 → `state=incomplete`，前缀保留，可续传；
* 源换版或摘要不符 → `409`，整批拒绝并清除暂存，`state=failed`。

对已经 `usable` 的 run 重复 import 是幂等空操作；清单未变的
`incomplete` run 重复 import 等于续传；清单变了（或之前失败）则清空
暂存按新清单重采。

### `POST /api/runs/{run_id}/resume`

从落盘字节处续传：请求带 `Range: bytes=<offset>-` 和
`If-Match: <冻结ETag>`。

* `206`：核对响应 ETag 与冻结清单一致、`Content-Range` 起点等于本地
  偏移、总长度等于清单长度，才追加写盘；
* `200`（源忽略 Range）：ETag 核对一致后**从头覆盖**该文件，绝不追加；
* `412` 或 ETag 变化：源换版 → 拒绝整批、清除暂存、`409`；
* 连接中断：保留已落盘前缀，`state=incomplete`，可再次续传。

两份文件都齐后：先逐份核对 SHA-256（不符 → 拒绝整批），再跑与上传
接口相同的 FASTQ 配对 + gzip 校验（不符 → `502 source_data_invalid`，
整批拒绝），全部通过才一起标为 `usable`。

### `GET /api/runs/{run_id}` / `GET /api/runs`

查询单个 run 或全部 run 的状态：每份文件的清单长度、已落盘字节、
摘要是否通过、最近错误和时间戳。

### `POST /api/runs/{run_id}/demux`

仅 `usable` 的 run 可用（否则 `409 run_not_usable`）。表单字段
`barcodes` 与上传接口相同，对落盘的两份文件执行同一套拆样逻辑，
返回同样的样本 ZIP。

### 状态机与错误码

```text
(无记录) --import--> incomplete --传输完成+校验通过--> usable
                       |  ^-- resume/重启后续传            |
                       |                                 (文件丢失 -> incomplete)
                       +-- 换版/摘要不符/FASTQ无效 --> failed --import--> incomplete
```

采集相关 `detail.code`：`invalid_run_id`(422)、`run_not_found`(404)、
`acquisition_in_progress`(409)、`run_failed`(409)、`run_not_usable`(409)、
`source_version_changed`(409)、`digest_mismatch`(409)、
`source_unavailable`(502)、`invalid_manifest`(502)、
`source_data_invalid`(502)。断流不算错误响应：返回 200 且
`state=incomplete`、`last_error.code=transfer_interrupted`。

### 可靠性设计

* **进度不超过落盘字节**：每写一块、flush 之后才更新 SQLite 进度；
  启动时（及每次传输前）把进度钳制到 `.part` 实际大小，崩溃留下的
  多余尾巴先截断再续传。
* **暂存与晋级**：下载中的文件是 `runs/{run_id}/*.part`；两份都通过
  摘要和 FASTQ/gzip 校验后，才原子改名为 `R1.fastq.gz`/`R2.fastq.gz`
  并在同一事务里把 run 标为 `usable`，不会出现"半个可用 run"。
* **整批拒绝**：换版、摘要不符、FASTQ 无效都会删除整个 run 的暂存
  目录并把 run 标为 `failed`，只能重新 import。
* **配置**：`INSTRUMENT_BASE_URL`（源地址，默认
  `http://localhost:9000`）、`ACQ_DATA_DIR`（SQLite + 暂存目录，默认
  `./data`，compose 中为 `/data` 卷）、`DEMUX_TMP_ROOT`（拆样临时目录）。

## 在 /docs 演示上传和下载

1. 打开 `/docs`，展开 `POST /api/demux`，点 **Try it out**。
2. `r1` 选择 `examples/normal/R1.fastq.gz`，`r2` 选择
   `examples/normal/R2.fastq.gz`。
3. `barcodes` 点 **Add item** 两次，分别填 `ACGTACGT` 和 `TTTTAAAA`
   （也可以在一个字段里填 `ACGTACGT, TTTTAAAA`）。
4. 执行后响应是 ZIP 文件，浏览器会下载；用解压工具查看
   `manifest.json`、`pair_counts.tsv` 和每个组的 R1/R2 gzip。

等价的命令行调用：

```bash
curl -sS -D headers.txt \
  -F 'r1=@examples/normal/R1.fastq.gz' \
  -F 'r2=@examples/normal/R2.fastq.gz' \
  -F 'barcodes=ACGTACGT' \
  -F 'barcodes=TTTTAAAA' \
  -o normal.zip \
  http://localhost:8000/api/demux
```

`normal.zip` 的期望结果：

| 组 | 配对数 | 原因 |
|---|---:|---|
| `ACGTACGT` | 3 | 两个精确匹配，一个 1 位错配且只接近该样本 |
| `TTTTAAAA` | 2 | 一个精确匹配，一个 1 位错配且只接近该样本 |
| `unassigned` | 1 | `CCCCCCCC` 与两个样本都相差超过 1 位 |

### 随附样例

| 目录 | 条码 | 期望 |
|---|---|---|
| `examples/normal` | `ACGTACGT`, `TTTTAAAA` | 正常拆分，包含精确匹配、1 位错配和未分配 |
| `examples/ambiguous` | `ACGTACGT`, `ACGTACGA` | 三个观测条码都同时落入两个样本的 1 位范围，全部进入 `unassigned` |
| `examples/missing_pair` | `ACGTACGT`, `TTTTAAAA` | R1 有 3 条、R2 只有 2 条，返回 `count_mismatch`，位置为 R1/record 3/line 9 |

每个目录都同时提供 `.fastq` 原文和可直接上传的 `.fastq.gz`，`barcodes.txt`
是建议的条码。重新生成样例：

```bash
python3 examples/generate_examples.py
```

## 拆样规则

### `POST /api/demux`

`multipart/form-data`：

| 字段 | 类型 | 说明 |
|---|---|---|
| `r1` | file | gzip 压缩的 R1 FASTQ |
| `r2` | file | gzip 压缩的 R2 FASTQ |
| `barcodes` | string，可重复 | 1–16 个唯一条码，每个恰好 8 位 `A/C/G/T`。重复字段或以逗号/换行分隔均可 |

成功时返回 `application/zip`。ZIP 内容：

```text
manifest.json
pair_counts.tsv
<barcode>/<barcode>_R1.fastq.gz
<barcode>/<barcode>_R2.fastq.gz
unassigned/unassigned_R1.fastq.gz
unassigned/unassigned_R2.fastq.gz
```

所有提交的样本条码（包括 0 条的组）和 `unassigned` 都会有各自的
R1/R2 gzip；0 条组是合法的空 gzip，便于下游脚本统一处理。
`manifest.json` 里有总配对数、每个组的 `pair_count` 和文件路径；
`pair_counts.tsv` 提供同一信息的表格版本。

输入要求：

* 严格四行一条的 FASTQ，最后一条也应有换行。
* R1 头：`@ID/1 barcode=<8 位观测条码>`；R2 头：`@ID/2`。
* 两个文件按顺序记录 ID 相同、条数相等。
* 每个读段长度 ≤ 1000，序列和 quality 长度严格相等。
* R1/R2 都是完整、未损坏的 gzip（读取结束时会校验 CRC32 和 ISIZE 尾部）。

### 条码判定

观测条码与每个样本条码逐位（Hamming）比较：

```text
候选 = {样本条码 | 与观测条码的 Hamming 距离 <= 1}
候选恰好 1 个 -> 归入该样本
否则（0 个或多于 1 个） -> unassigned
```

精确匹配也参与这个规则：如果观测条码正好等于样本 A，但它同时与样本 B
只差 1 位，结果为 `unassigned`（歧义），不会强行归入 A。R1/R2 作为一整对
读段处理，不会拆开；每组内保持输入出现顺序，并逐字节保留原始四行文本。

### 错误响应

坏记录、ID 错配、条数不等或 gzip 损坏返回 HTTP 400；条码等请求参数问题
返回 HTTP 422。响应格式：

```json
{
  "detail": {
    "code": "count_mismatch",
    "message": "R1 has an extra record at pair 3, ...",
    "file": "R1",
    "record": 3,
    "line": 9,
    "locations": {
      "R1": {"record": 3, "line": 9},
      "R2": {"record": 3, "line": 9}
    }
  }
}
```

`code` 包括 `bad_record`、`truncated_record`、`id_mismatch`、
`count_mismatch`、`length_mismatch`、`read_too_long`、`gzip_error` 和
`invalid_barcodes`。`record` 是从 1 开始的读段序号，`line` 是从 1 开始的
物理行号（四行一条时第 N 条头行是 `4N-3`）。错误发生在构建阶段，响应
不会包含半成品 ZIP。

## 资源与内存

* 读取：每次只在内存中保留一对 R1/R2 记录（单行设上限），不建读段列表。
* 分组：最多 16 个样本加 `unassigned`，使用固定数量的 gzip 临时文件，
  逐个读段对追加；内存与读段总数无关。
* 打包：ZIP 先落盘，再以 128 KiB 块流式发送；`ZIP_STORED` 不重复压缩
  已压缩的组文件，内存同样不随读段数增长。
* 清理：构建失败时立即删除本请求临时目录；下载正常结束、写响应失败或
  客户端断开时，响应生成器的 `finally` 和幂等后台清理都会删除临时目录；
  上传 spool 文件也随即关闭。`DemuxPackage` 还有析构兜底。
* 采集：按传输块落盘，内存只有当前块；进度和清单在 SQLite，暂存文件
  在 `ACQ_DATA_DIR`，进程重启后原样续传。

临时目录默认 `/tmp/demux`，可用环境变量 `DEMUX_TMP_ROOT` 修改。

## 测试

```bash
pip install -r requirements-dev.txt
pytest -q
```

核心测试只用标准库；API 测试需要 FastAPI/httpx。采集测试
（`tests/test_acquisition.py`）自带一个可控的仪器源桩，覆盖断流续传、
换版拒绝、忽略 Range 覆盖重下、摘要不符、进度与落盘一致性等场景。
当前仓库附带的 `tests/test_demux.py` 也可不装依赖直接运行：

```bash
python3 -m unittest tests.test_demux -v
```
