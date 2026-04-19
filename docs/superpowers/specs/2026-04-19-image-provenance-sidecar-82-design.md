# #82 图像 sidecar provenance 设计（通用 Python 模块方案）

## 背景与目标

Issue #82 目标是在受控流程内，为图像产物建立可验证 sidecar 溯源。

本次确认的试点链路：`scripts/plot_scan_csv.py`。

本次采用 B 方案：抽取通用 Python 模块，支持后续多脚本复用，而不是仅在单脚本内内聚实现。

## 非目标

- 不在图像文件本体嵌入元数据。
- 不在本 issue 覆盖所有历史绘图脚本。
- 不引入 C2PA/区块链/水印等开放传播溯源机制。

## 方案总览

### 新增通用模块

新增：`scripts/common/provenance_image.py`

职责：

1. 计算文件哈希（SHA256）。
2. 构造 v1 图像 sidecar 元数据对象。
3. 写入 `<image>.provenance.json`。
4. 对图像 + sidecar 做强绑定校验并返回 pass/fail 与原因。

### 模块接口（v1）

建议对外函数：

- `compute_sha256(path: Path) -> str`
- `build_image_provenance(...) -> dict`
- `write_image_provenance_sidecar(image_path: Path, metadata: dict) -> Path`
- `verify_image_provenance(image_path: Path) -> tuple[bool, str]`

## v1 字段规范

sidecar 文件命名：`<image_basename>.provenance.json`

必填字段（与 #82 对齐）：

- `schema_version`
- `generated_at_utc`
- `image_path`
- `image_sha256`
- `script_path`
- `command`
- `git_commit`
- `julia_version`
- `input_data_hashes`（数组，元素含 `path`、`sha256`）

可选字段：

- `config_path`
- `config_hash`
- `seed`
- `artifact_paths`
- `notes`

## 在 `plot_scan_csv.py` 的接入方式

1. 在所有 `savefig` 落点后统一调用 sidecar 写入。
2. 默认开启 sidecar 输出（不改现有图像命名与目录逻辑）。
3. 收集并传入上下文：
   - 当前命令行（command）
   - `git_commit`
   - 解释器/运行时版本映射到 `julia_version` 字段（按 issue 字段要求保留该键）
   - 输入 CSV 的哈希清单
4. sidecar 与图像同目录落盘。

## 校验入口

在 `plot_scan_csv.py` 增加 CLI 校验模式（如 `--verify-provenance <image>`）：

- 读取同名 sidecar。
- 复算图像 `sha256`。
- 返回：
  - `PASS`：哈希一致且必填字段齐全。
  - `FAIL`：sidecar 缺失 / 字段缺失 / 哈希不匹配 / JSON 非法。

## 测试策略（TDD）

新增 Python 单测（放在仓库现有 Python 测试约定路径）：

1. 生成临时图像 + 写 sidecar 成功。
2. 校验通过（pass）。
3. 篡改图像后校验失败（fail）。
4. sidecar 缺字段时校验失败并给出明确原因。

## 文档更新

在脚本说明或相关 README 增加：

- sidecar 的边界说明：仅保证受控环境复现，不保证开放传播链路。
- 归档要求：图像与 `*.provenance.json` 必须同存档。
- 校验命令示例。

## 验证命令

按 issue #82 给定：

```bash
julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'
julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'
```

并补充 Python 单测命令（以仓库实际测试入口为准）。

## 风险与缓解

1. 图像与 sidecar 分离造成断链
   - 缓解：校验入口强制 `image_sha256` 绑定；文档明确同存档要求。

2. 字段命名与现有 provenance 风格不一致
   - 缓解：对齐 `scripts/relaxtime/provenance_metadata.jl` 的组织风格（版本、哈希、运行上下文）。

3. 仅覆盖单条链路导致代表性不足
   - 缓解：优先选择 `scripts/plot_scan_csv.py` 这种通用绘图入口作为首条落地。

## DoD 对齐

- 至少 1 条图像链路（`plot_scan_csv.py`）落地 sidecar。
- `*.provenance.json` 覆盖 v1 必填字段。
- 校验入口能识别 pass/fail（含篡改检测）。
- 文档中明确 sidecar 方案边界。
