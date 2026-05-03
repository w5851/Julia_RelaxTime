# 系统状态总结

本页用于给出“当前可用能力 + 已知风险”，状态口径与 `README.md` 保持一致。

## 1. 总体状态

- **Web/API 演示链路**：可用，适合联调与能力展示
- **PNJL 相图与统一扫描链路**：可用，稳定脚本入口已统一到白名单
- **截面/弛豫时间链路**：可用，建议持续执行关键点回归与跨实现对比

## 2. 关键入口

- 稳定脚本白名单：`docs/guides/scripts/README.md`
- 相图主产线：`scripts/pnjl/calculate_phase_structure.jl`
- 统一扫描入口：`scripts/models/run_unified_scan.jl`
- 服务端：`scripts/server/server_full.jl`
- 前端页面：`web/index.html`

如需查看更深层实现说明：

- 相图 API：`docs/api/models/phase/README.md`
- 扫描 API：`docs/api/models/scans/README.md`
- workflow API：`docs/api/models/workflows/README.md`
- transport API：`docs/api/relaxtime/transport/README.md`

## 3. 推荐验证命令

```powershell
# 相图最小产线
julia --project=. scripts/pnjl/calculate_phase_structure.jl --preset=smoke --output_dir=data/outputs/results/phase_smoke

# 单元 smoke
julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'

# 文档与入口治理
julia --project=. scripts/dev/check_docs_consistency.jl
julia --project=. scripts/dev/check_script_entrypoints.jl
```

## 4. 已知注意事项

- 稳定用户入口以 `docs/guides/scripts/README.md` 为准；并非所有 `run_*.jl` 都自动属于白名单。
- API 主题主导航以 `docs/api/models/*` 与 `docs/api/relaxtime/*` 为准；历史 `docs/api/pnjl/` 页面只保留兼容层说明。
- 不建议将“前端可用”解读为“全部参数区间与全部场景均已覆盖验收”。
- 对研究结论请优先参考 `docs/reference/` 与 `docs/dev/archived/` 的比对记录。

## 5. 输出目录口径

- 默认运行产物目录：`data/outputs/`
- 根目录 `outputs/`：仅历史兼容，不作为默认结果落盘目录
