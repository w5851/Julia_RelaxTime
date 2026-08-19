# Analysis Governance

本目录收纳 `docs/analysis` 分析产物的治理、清理和 provenance 记录，不承载新的科学结果。

## Figure asset registry

[`figure_asset_registry_v1/`](figure_asset_registry_v1/) 是历史 PNG/PDF/SVG 资产 inventory、作者审核和 allowlist 执行快照。其 registry、preflight、retirement 和 relocation JSON 保留生成时的路径、hash 和状态，不因本次 namespace 迁移重写。

治理快照不授权删除、移动或覆盖正式图像；任何后续资产动作仍需独立的作者批准和 allowlist 审阅。
