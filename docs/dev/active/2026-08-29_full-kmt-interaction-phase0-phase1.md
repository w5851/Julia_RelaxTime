# 完整 KMT 介子相互作用核 Phase 0/1

更新日期：2026-08-29

当前状态：Phase 0 与 Phase 1 已完成；Phase 2（完整中性 RPA 极化接入）尚未开始。本任务只建立公式合同和纯后端，不修改上游 PNJL 平衡求解或既有介子数密度生产链。

## 1. 目标与边界

本任务把完整三味 KMT 四费米有效相互作用核作为一个与旧接口并行的后端，首先覆盖中性
`(0,3,8)` 基底的 `K_{03}=K_{30}`、`K_{38}=K_{83}`，以及味道分辨的带电
`K_{12}`、`K_{45}`、`K_{67}`。输入只来自平均场背景中的
`(phi_u, phi_d, phi_s, G, K)`。

本阶段明确不做：

- 修改 `PNJLCore` 的平均场方程、非对角夸克自能或凝聚变量；
- 修改旧 `EffectiveCouplings` 的 `NamedTuple` 合同；
- 修改旧 `MesonPropagator` 的 `0/8` 二维接口；
- 接入完整 `0/3/8` RPA 极化矩阵、介子相移或 `Omega_M` 热力学反馈；
- 将新后端的数值结果标为 production。

## 2. Phase 0：公式与约定冻结

- [x] 记录从 KMT 六费米项到二次 RPA 相互作用核的层次边界。
- [x] 锁定 `lambda_0=sqrt(2/3) I`、中性基底 `(0,3,8)`、`phi_f` 为物理夸克凝聚以及 P/S 通道符号约定。
- [x] 写出同位旋对称极限 `phi_u=phi_d` 下 `K_{03}=K_{38}=0` 与 `K_{45}=K_{67}` 的验收关系。
- [x] 将候选公式、来源和“符号依赖约定、尚需独立文献审计”的状态写入公式参考文档。

公式合同见：
`docs/reference/formula/relaxtime/couplings/KMT_MFA_to_RPA_QuadraticKernel.md`。

## 3. Phase 1：并行纯后端

- [x] 新增 `MesonInteractionKernel` 模块和不可变 `FullKMTInteraction` 数据结构。
- [x] 实现从 `SVector`、长度为 3 的向量/元组以及 equilibrium `x_state[1:3]` 构造完整核。
- [x] 暴露中性矩阵、带电味道通道和 `K_ab` 查询函数；保持旧模块 API 不变。
- [x] 将模块按依赖顺序接入 `RelaxTime`，补充 API 文档。
- [x] 增加对称极限、KMT 关闭、P/S 通道互补、旧对角耦合兼容和非法输入测试。

## 4. 验证门禁与后续入口

Phase 1 只要求纯代数和接口测试通过。后续 Phase 2 才能把该后端接到完整中性 RPA 极化矩阵；接入前必须再次核对 `K_{03}`、`K_{38}` 的整体符号、凝聚定义和 KMT 耦合 `K` 的配置符号。Phase 3 才评估完整核对上游凝聚/质量的影响，Phase 4 才做带电 K/π 的 A/B 数值对照。

允许继续保留旧 `EffectiveCouplings` 与旧 `MesonPropagator` 的生产/回归语义。新后端的任何数值扫描必须使用独立 case slug，并标注 `diagnostic`。

## 5. 本阶段验证记录

- `julia --project=. tests/unit/relaxtime/test_meson_interaction_kernel.jl`：52/52 通过。
- `julia --project=. -e 'include("src/relaxtime/RelaxTime.jl"); ...'`：`RelaxTime` include/re-export smoke 通过。
- `git diff --check`：通过。
- 单元测试覆盖：不对称 `K_{03}/K_{38}`、同位旋对称退化、KMT 关闭、P/S 互补、旧 0/8 耦合映射、equilibrium `x_state` 适配、非法输入与通道查询。

上述验证只证明纯代数后端的接口和内部公式合同，不证明候选符号已经完成外部论文逐项复核，也不证明完整 RPA/介子密度结果的物理正确性。

## 6. 证据和引用

- 项目背景调研：`docs/analysis/relaxtime/meson_isospin_workflow_literature_survey_2026-08-28.md`。
- 公式与实现约定：`docs/reference/formula/relaxtime/couplings/KMT_MFA_to_RPA_QuadraticKernel.md`。
- API 合同：`docs/api/relaxtime/propagator/MesonInteractionKernel.md`。
- 单元测试：`tests/unit/relaxtime/test_meson_interaction_kernel.jl`。
