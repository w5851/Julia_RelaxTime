# PNJL Physics Simulation Legacy Reuse Audit

## Scope

本记录审计已归档的旧目录 `Julia/PNJL_Physics_Simulation`，判断其是否仍有内容可以直接进入 `Julia_RelaxTime` 主线。

- 来源仓库：`https://github.com/w5851/Julia.git`
- 来源仓库 HEAD：`4de0865`
- safety export：`D:\Desktop\_cleanup_quarantine\2026-05-27\_safety_exports\Julia\`
- Git bundle SHA-256：`A8A7730B3AD57338ED19288638C31D3D995176AABB6648351D96FD2C94513C4F`
- 审计日期：`2026-08-16`

审计期间可复核的归档内容曾是 Git bundle、tracked working-tree diff 和未跟踪路径清单；
本项已按用户授权清理 `safety export`，上述路径当前不再存在。`_safety_exports\full_archives\`
在清理前即为空，因此不存在可复核的完整 working-tree tarball。
bundle 可恢复完整 tracked Git 历史（审计 HEAD=`4de0865`），但未跟踪文件正文不在
safety export 中，不能把它们当成该 dirty 工作区的完整快照。桌面上发现的同名文件属于
另一个 `w5851/PNJL_Simulation` 仓库且哈希不同，只能作为相关项目线索，不能替代本归档。

该目录是一个早期 include/package 混合式 PNJL 重构项目，包含 `PNJL`、`PNJL_aniso`、`Rotation` 和 `GasLiquid` 四类旧模型实现，以及一套通用积分、方程工厂和 phase scan 接口。

## 主要内容与当前对应物

### 通用积分接口

旧项目的 `src/core/integration_interface.jl` 提供了：

- `IntegrationMethod`、`IntegrationGrid` 抽象类型；
- `MomentumGrid`、`AngleGrid`、`ProductGrid`；
- `integrate`、`integrate_2d`、`discrete_sum` 和 `angular_momentum_sum`；
- Gauss-Legendre 与一个自适应积分方法的统一外壳。

当前主项目已经在不同边界上实现了更明确的职责划分：

- `src/integration/GaussLegendre.jl`；
- `src/integration/IntervalQuadratureStrategies.jl`；
- `src/models/variants/*/core/` 中的模型专用积分；
- `Models` 入口与各扫描 workflow。

旧接口依赖 include 顺序、宽泛的 `Function` 参数和早期配置类型；其中自适应积分也不能绕过当前可选数值 oracle 的隔离规则。因此不能直接复制到 `src/`。

### 方程工厂与统一公共接口

旧项目的 `equation_factory.jl` 和 `unified_physics_public_interface.jl` 尝试把方程绑定、平衡求解、物理性质提取和 phase scan 放在一个公共层。其主要概念已经被当前主线更严格地拆开：

- `src/models/Models.jl` 与 `src/models/entrypoints.jl`：统一公开入口；
- `src/models/phase/`：phase pipeline、CEP、Maxwell/branch governance；
- `src/models/scans/ScanCommon.jl`：续扫键、候选 seed、失败语义与 pressure 选择；
- `tests/unit/models/`、`tests/integration/models/` 和 `tests/regression/models/`：分层验证。

旧项目的 `PhaseScanInterface.scan_phase_space` 没有当前主线所需的 provenance、失败点治理、候选选择契约和 production/research 区分，不能直接移植。

### 旧模型实现与报告

旧目录中的 Gas-Liquid、PNJL-aniso 和 Rotation 文件是早期 include 驱动实现，当前主线已有对应的 `Models` 变体和更近的模型核心。归档中的 modernization/performance 报告给出“测试通过率”“性能提升”“机器精度”等结论，但没有当前主项目的 commit、配置指纹、脚本 hash 或可复现 manifest；这些数字只能作为历史说明，不能进入 baseline、回归测试或论文证据。

## 决定

- **直接迁移代码：不执行。** 旧接口与当前 `Models`/scan/数值治理边界不兼容。
- **迁移旧输出：不执行。** 归档中的 CSV、DAT 和 phase scan 输出没有当前 provenance，不进入 `data/outputs/` 或 regression target。
- **立即可做的主项目贡献：完成本审计文档。** 它记录旧架构与当前对应物，避免清理后重复评估或误迁移。
- **未来可能的独立任务：** 若确实需要通用多维积分或 Gas-Liquid 参数反演，应基于当前稳定入口重新设计 contract、单位、失败语义和验证层，而不是恢复旧模块。

## 清理条件

本归档没有未完成的主线迁移依赖。该项的 review 解包目录已清理；
由于完整 working-tree tarball 实际不存在，也不应声称未跟踪正文已被完整保全。
主项目当前已有的用户 dirty 文件不在清理范围内。
