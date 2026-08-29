# 完整 KMT 介子相互作用核 Phase 0/1

更新日期：2026-08-29

当前状态：Phase 0、Phase 1、Phase 1.5 与 Phase 2（中性 RPA 代数后端）已完成；尚未接入现有夸克泡数值工作流。本任务只建立来源对齐的公式合同和并行后端，不修改上游 PNJL 平衡求解或既有介子数密度生产链。

## 1. 目标与边界

本任务把完整三味 KMT 四费米有效相互作用核作为一个与旧接口并行的后端，首先覆盖中性
`(0,3,8)` 基底的 `K_{03}=K_{30}`、`K_{38}=K_{83}`，以及味道分辨的带电
`K_{12}`、`K_{45}`、`K_{67}`。输入只来自平均场背景中的
`(phi_u, phi_d, phi_s, G, K)`。

本阶段明确不做：

- 修改 `PNJLCore` 的平均场方程、非对角夸克自能或凝聚变量；
- 修改旧 `EffectiveCouplings` 的 `NamedTuple` 合同；
- 修改旧 `MesonPropagator` 的 `0/8` 二维接口；
- 接入现有夸克泡数值工作流、介子相移或 `Omega_M` 热力学反馈；
- 将新后端的数值结果标为 production。

## 2. Phase 0：公式与约定冻结

- [x] 记录从 KMT 六费米项到二次 RPA 相互作用核的层次边界。
- [x] 锁定 `lambda_0=sqrt(2/3) I`、中性基底 `(0,3,8)`、`phi_f` 为物理夸克凝聚以及 P/S 通道符号约定。
- [x] 写出同位旋对称极限 `phi_u=phi_d` 下 `K_{03}=K_{38}=0` 与 `K_{45}=K_{67}` 的验收关系。
- [x] 将公式、来源和“约定依赖、数值接入仍需审计”的状态写入公式参考文档。

公式合同见：
`docs/reference/formula/relaxtime/couplings/KMT_MFA_to_RPA_QuadraticKernel.md`。

## 3. Phase 1：并行纯后端

- [x] 新增 `MesonInteractionKernel` 模块和不可变 `FullKMTInteraction` 数据结构。
- [x] 实现从 `SVector`、长度为 3 的向量/元组以及 equilibrium `x_state[1:3]` 构造完整核。
- [x] 暴露中性矩阵、带电味道通道和 `K_ab` 查询函数；保持旧模块 API 不变。
- [x] 将模块按依赖顺序接入 `RelaxTime`，补充 API 文档。
- [x] 增加对称极限、KMT 关闭、P/S 通道互补、旧对角耦合兼容和非法输入测试。

## 3.5 Phase 1.5：外部来源对齐

- [x] 以 Tian 等，Phys. Rev. D 114, 034012 (2026)，Eq. (3) 为当前实现的明确代数来源。
- [x] 按 `phi_f = sigma_f` 和 P/S 约定修正 `K_{03}`、`K_{38}` 的整体符号。
- [x] 更新公式合同、`KernelConvention` 元数据和不对称凝聚单元测试。

## 3.6 Phase 2：中性 RPA 代数后端

- [x] 新增 `MesonRPA`，将 `(Pi_u, Pi_d, Pi_s)` 组装为 `(0, 3, 8)` 极化矩阵。
- [x] 实现并导出 `I - 2 K Pi`、`2 K * inv(I - 2 K Pi)` 和其行列式。
- [x] 锁定矩阵乘法顺序、同位旋对称解耦、KMT 关闭但极化仍可非对角等边界行为。
- [x] 补充 API 文档与 31 项纯代数单元测试。

## 4. 验证门禁与后续入口

Phase 1.5/2 只要求来源对齐、纯代数和接口测试通过。当前 `MesonRPA` 接受调用方已计算好的三味夸克泡，不负责把 `PolarizationAniso` 的某一数值约定自动解释成论文 Eq. (26) 的 `Pi_f`。后续 Phase 3 才评估完整核对上游凝聚/质量和极化归一化的影响，Phase 4 才做带电 K/π 的 A/B 数值对照。

允许继续保留旧 `EffectiveCouplings` 与旧 `MesonPropagator` 的生产/回归语义。新后端的任何数值扫描必须使用独立 case slug，并标注 `diagnostic`。

## 5. 本阶段验证记录

- `julia --project=. tests/unit/relaxtime/test_meson_interaction_kernel.jl`：54/54 通过。
- `julia --project=. tests/unit/relaxtime/test_meson_rpa.jl`：31/31 通过。
- `julia --project=. -e 'include("src/relaxtime/RelaxTime.jl"); ...'`：`RelaxTime` include/re-export smoke 通过。
- `git diff --check`：通过。
- 单元测试覆盖：不对称 `K_{03}/K_{38}`、同位旋对称退化、KMT 关闭、P/S 互补、旧 0/8 耦合映射、equilibrium `x_state` 适配、来源对齐的 `(0,3,8)` 极化矩阵、RPA 矩阵顺序、KMT 关闭下的极化非对角项、非法输入与通道查询。

上述验证只证明来源对齐的纯代数后端接口和内部公式合同，不证明 `PolarizationAniso` 数值极化与外部论文逐项等价，也不证明完整 RPA/介子密度结果的物理正确性。

## 6. 证据和引用

- 项目背景调研：`docs/analysis/relaxtime/meson_isospin_workflow_literature_survey_2026-08-28.md`。
- 外部结构来源：Tian 等，Phys. Rev. D 114, 034012 (2026)，DOI：[10.1103/d7nm-y2vp](https://doi.org/10.1103/d7nm-y2vp)。该文的 NJL+磁场/Ritus/Pauli-Villars 和极点质量流程不作为当前 PNJL/BU 数值基线。
- 公式与实现约定：`docs/reference/formula/relaxtime/couplings/KMT_MFA_to_RPA_QuadraticKernel.md`。
- API 合同：`docs/api/relaxtime/propagator/MesonInteractionKernel.md`。
- 单元测试：`tests/unit/relaxtime/test_meson_interaction_kernel.jl`。
- 中性 RPA 单元测试：`tests/unit/relaxtime/test_meson_rpa.jl`。
