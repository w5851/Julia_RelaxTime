# 相图主题 API

本目录是相图主题在新架构下的主入口，面向 `src/models/phase/` 实现域与 `Models` 聚合入口之间的 API 文档治理。

## 阅读顺序

- 面向用户入口：见 [Overview.md](docs/api/models/phase/Overview.md)
- 算法核心：见 [Algorithms.md](docs/api/models/phase/Algorithms.md)
- 导出 API 全集：见 [generated/Exports.md](docs/api/models/phase/generated/Exports.md)

## 适用范围

本主题覆盖以下能力：

- `Models.run_phase_pipeline` 的主工作流调用
- `Models.find_cep` 的独立 CEP 诊断
- `Models.build_phase_artifacts`、`Models.resolve_phase_output_target`、`Models.promote_phase_artifacts` 的工件治理
- `Models.analyze_pm_branch_competition` 的 compare-only `P-mu` 诊断入口
- `src/models/phase/` 下的关键算法关系，例如 S-shape、Maxwell、crossover 与自适应 `rho` 加密

## 目录职责

- [Overview.md](docs/api/models/phase/Overview.md)：给首次使用者的入口页，先讲如何跑通完整流程
- [Algorithms.md](docs/api/models/phase/Algorithms.md)：给维护者或需要理解判据关系的读者
- [PhaseTransition.md](docs/api/models/phase/PhaseTransition.md)：S-shape、Maxwell 与温度切片整理
- [Crossover.md](docs/api/models/phase/Crossover.md)：crossover 检测与扫描
- [AdaptiveRhoRefinement.md](docs/api/models/phase/AdaptiveRhoRefinement.md)：自适应 `rho` 加密辅助层
- [PMPhaseDiagnostic.md](docs/api/models/phase/PMPhaseDiagnostic.md)：compare-only `P-mu / Omega-mu` 双分支竞争诊断
- [generated/Exports.md](docs/api/models/phase/generated/Exports.md)：自动生成的公开导出索引，作为完整性基线

## 入口约束

- 相图主题的唯一主入口为本页与 [Overview.md](docs/api/models/phase/Overview.md)
- 不再保留旧 `pnjl` 相图页作为并列导航入口
