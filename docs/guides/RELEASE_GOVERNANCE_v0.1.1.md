# v0.1.1 Patch 发布准入清单（防漂移）

本清单用于约束 `v0.1.1` 作为 **patch** 发布：默认仅允许修复与文档收敛，不引入公共能力语义漂移。

## 1. 版本定位

- 版本类型：`patch`
- 允许内容：
  - bug fix（行为修复）
  - 文档修正（路径/契约/示例）
  - 非破坏性治理增强（检查脚本、测试补强）
- 不允许内容：
  - 新的公共 API 面
  - 既有公共契约字段删改
  - 默认行为的未声明改变

## 2. 必过门禁（与 v0.1.0 对齐）

- [ ] `julia --project=. -e 'ENV["UNIT_FILES"]="simulation/test_fullserver_compute_handlers.jl;simulation/test_fullserver_pnjl_handlers.jl;simulation/test_pnjl_scan_jobs.jl"; include("tests/unit/runtests.jl")'`
- [ ] `julia --project=. -e 'ENV["UNIT_FILES"]="relaxtime/test_transport_coefficients.jl;models/test_transport_workflow.jl"; include("tests/unit/runtests.jl")'`
- [ ] `julia --project=. scripts/dev/check_docs_consistency.jl`
- [ ] `julia --project=. scripts/dev/archive_docs.jl -c`

## 3. 契约不漂移检查（必须逐条确认）

### 3.1 FullServer 错误响应契约

- [ ] 错误最小字段仍为：`status/error_code/error/message_id`
- [ ] `pnjl_scan_jobs` 对外不暴露内部堆栈文本

### 3.2 scan jobs 治理能力

- [ ] TTL/keep_max 仍可通过环境变量配置
- [ ] `governance` 指标字段仍可观测且语义不变

### 3.3 run_脚本能力合同（Mott workflow）

- [ ] `run_mott_phase_scan.jl` / `run_mott_phase_derived_csv.jl` / `run_mott_phase_plot_modes.jl` 仍可按文档执行
- [ ] `effective_config.json` 与 `run_manifest.json` 仍产出

## 4. Release Notes 最小结构（v0.1.1）

发布说明建议固定为以下四段：

1. `Fixes`
2. `Contract Stability`
3. `Verification`
4. `Upgrade Notes`

## 5. 触发 minor/major 的升级条件

若出现以下情况，不应继续使用 patch 版本号：

- 增加新的稳定入口脚本或公开 API -> 至少 `minor`
- 改变既有错误响应字段语义或删除字段 -> 至少 `major`（或显式迁移期）
- 改变默认工作流行为且影响复现结果 -> 至少 `minor`
