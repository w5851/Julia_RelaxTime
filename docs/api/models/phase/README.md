# 相图主题 API

本目录是相图主题在新架构下的主入口，面向 `src/models/phase/` 实现域与 `Models` 聚合入口之间的 API 文档治理。

## 阅读顺序

- 面向用户入口：见 [Overview.md](docs/api/models/phase/Overview.md)
- 算法核心：见 [Algorithms.md](docs/api/models/phase/Algorithms.md)
- 导出 API 全集：见 [generated/Exports.md](docs/api/models/phase/generated/Exports.md)

## 适用范围

本主题覆盖以下能力：

- `Models.run_phase_pipeline` 的通用 / research 主工作流调用
- `Models.run_production_phase_pipeline` 的 production / baseline 主工作流调用
- `Models.find_cep` 的独立 CEP 诊断
- `Models.build_phase_artifacts`、`Models.resolve_phase_output_target`、`Models.promote_phase_artifacts` 的工件治理
- `Models.analyze_pm_branch_competition` 的 compare-only `P-mu` 诊断入口
- `src/models/phase/` 下的关键算法关系，例如 S-shape、Maxwell、crossover 与自适应 `rho` 加密
- `Models.RhoSupportConfig` 与 opt-in `rho_support_cascade` / `rho_support_hybrid` 的分层几何证书和缓存诊断
- `Models.RhoHybridVerificationConfig` 的离散极值外侧 guard、端点局部几何证书与 Stage-C 局部验证合同

## 目录职责

- [Overview.md](docs/api/models/phase/Overview.md)：给首次使用者的入口页，先讲如何跑通完整流程
- [Algorithms.md](docs/api/models/phase/Algorithms.md)：给维护者或需要理解判据关系的读者
- [PhaseTransition.md](docs/api/models/phase/PhaseTransition.md)：S-shape、Maxwell 与温度切片整理
- [Crossover.md](docs/api/models/phase/Crossover.md)：crossover 检测与扫描
- [AdaptiveRhoRefinement.md](docs/api/models/phase/AdaptiveRhoRefinement.md)：自适应 `rho` 加密辅助层
- [PMPhaseDiagnostic.md](docs/api/models/phase/PMPhaseDiagnostic.md)：compare-only `P-mu / Omega-mu` 双分支竞争诊断
- [generated/Exports.md](docs/api/models/phase/generated/Exports.md)：自动生成的公开导出索引，作为完整性基线

## 双入口口径

- `Models.run_phase_pipeline`
  - 面向通用相图工作流
  - 默认适合研究态比较、策略试验和 `find_cep` / crossover 组合分析
- `Models.run_production_phase_pipeline`
  - 面向 production/baseline 同步
  - 显式承载温度扫掠、unknown budget、fallback 与高精度 CEP 收口逻辑

两者都返回 `PhasePipelineResult`，但 production 入口会额外依赖 `FirstOrderSweepResult` 与 `ProductionPipelineConfig` 这两个稳定类型来表达扫掠合同。

`ProductionPipelineConfig.rho_refinement_policy` 默认是 `:uniform_nested`。只有显式
选择 `:rho_support_cascade` 或 `:rho_support_hybrid` 才启用 PNJL rho-support 路由；它们不改变默认 uniform
production，也不允许绕过 geometry gate 或将单层 no-S-shape 提升为 monotone 证书。
`rho_support_hybrid` 是五级链的显式 shadow/候选策略：Stage A cascade、Stage B
memoized dense，冲突时才在离散极值外侧 guard 内执行 Stage C `Δrho=0.003125`
局部验证。guard 是两个 μ 极值外侧的首个严格 Stage-B 采样点；不插值、不二分、
不使用固定 padding。Stage-C 始终保留完整 Stage-B 全域曲线，并只追加 guard 内按
Stage-B 特征排序选出的点。

当 Stage-B 的公共 Maxwell 结果是唯一三交点且 endpoint-dependent 时，可显式设置
`RhoHybridVerificationConfig(endpoint_policy=:three_crossing_endpoint_local_v2)`。该 v2
策略用两个实际 Stage-B 外支点 bracket 右 Maxwell crossing，不把右支越过 μ 极值作为
额外条件；完整 Stage-B 曲线始终保留，Stage-C 只在左 crossing 的 active bracket 内加入
midpoint，并把 anchor 与最多 12 个 refinement 点分别写入诊断。左 bracket 下界保持为零并
达到预算时证书类型为 `endpoint_limited_first_order`（对外 `rho_hadron=0.0`，另存上界）；
下界变为正值且连续两级 geometry 通过时为
`endpoint_local_geometry_first_order`（保留正的 `rho_hadron`）。两者都只属于内部诊断，
物理输出仍使用三态合同；候选不唯一、右 crossing 无实际 bracket、solver/geometry 失败或
预算耗尽均保持 `ambiguous_near_critical`。旧的
`endpoint_policy=:bounded_zero_density_v1` 继续保留以复现历史产物，默认 uniform/cascade
语义不变。

production 的 Maxwell 容差由统一合同派生：内部二分 tolerance 为所有 active area
acceptance gates 最小值的 `0.1` 倍；rho/temperature geometry gate 仍独立记录 coarse/fine
离散收敛误差。有效 solver tolerance 会写入 config snapshot/hash 和 diagnostics，默认
`maxwell_construction` 调用口径保持兼容。

## 入口约束

- 相图主题的唯一主入口为本页与 [Overview.md](docs/api/models/phase/Overview.md)
- 不再保留旧 `pnjl` 相图页作为并列导航入口
