# 系统状态总结

本页用于给出“当前可用能力 + 已知风险”，状态口径与 `README.md` 保持一致。

Latest Release: [v0.1.0](https://github.com/w5851/Julia_RelaxTime/releases/tag/v0.1.0)

## 1. 总体状态

- **Web/API 演示链路**：可用（`/health`、`/compute`）
- **PNJL 求解与扫描链路**：可用（建议结合对比报告验证数值一致性）
- **截面/弛豫时间链路**：已验证可用（建议持续执行关键点回归与跨实现对比）

## 2. 关键入口

- 服务端：`scripts/server/server_full.jl`
- 一键启动：`scripts/server/start.bat`
- 前端页面：`web/index.html`
- 单元测试入口：`tests/unit/runtests.jl`

## 3. 推荐验证命令

```powershell
# 启动服务
julia --project=. scripts/server/server_full.jl

# 运行单元 smoke
$env:UNIT_PROFILE='smoke'
julia --project=. tests/unit/runtests.jl
```

## 4. 已知注意事项

- 不建议将“前端可用”解读为“全部参数区间与全部场景均已覆盖验收”
- 对研究结论请优先参考 `docs/reference/` 与 `docs/dev/archived/` 的比对记录
- 发生路径疑问时，统一以仓库根目录当前结构为准（`scripts/server/`、`tests/unit/`、`docs/reference/`）

## 4.1 近期变更（2026-02-19）

- PNJL/rPNJL 配置链路已补“关键参数校验 + 可控日志开关”，用于异常回溯：
	- `PNJL_CONFIG_LOG=1`：输出 PNJL 配置解析来源
	- `RPNJL_CONFIG_LOG=1`：输出 rPNJL 配置解析来源
- `docs/dev/active/2026-02-19_PNJL集成方向提炼待办.md` 已完成并归档到 `docs/dev/archived/`。

## 5. 输出目录口径

- 默认运行产物目录：`data/outputs/`
- 根目录 `outputs/`：仅历史兼容，不作为默认结果落盘目录
