# Models 公式索引（主线内化版）

本目录用于沉淀 `Models` 相关公式文档，强调三点：

- 与当前仓库的单位/命名口径一致（内部自然单位，外部输入显式 `*_MeV`）。
- 公式不做“原样搬运”，而是提炼成可直接映射到本仓库实现与测试的表达。
- 对尚未落地的模型能力，先给出“理论定义 + 迁移接口约束”，避免与现有主线冲突。

## 当前主题

- `njl/NJL_core.md`：NJL 核心公式与能隙方程。
- `pnjl/Omega_各向同性.md`：PNJL 各向同性巨热力学势。
- `pnjl/Omega_RS各向异性.md`：PNJL 动量各向异性热项改造。
- `rpnjl/rPNJL_core.md`：含矢量相互作用的 rPNJL 核心关系。
- `pnjl_magnetic/PNJL_magnetic_core.md`：外磁场 PNJL 核心关系。

## 本批新增（PR1）

- `gas_liquid/GasLiquid_RMF_CoreEquations.md`：Gas-Liquid (RMF/Walecka) 核心方程与约束求解口径。
- `rotation/Rotation_PNJL_CoreEquations.md`：Rotation-PNJL 在 `T-mu-omega` 语义下的最小公式集。
- `shared/Lagrangian_to_GrandPotential_MFA_Workflow.md`：从拉氏量到巨热力学势的统一平均场推导流程。
- `shared/OmegaTotal_并入介子压强后的统一AD热力学流程.md`：固定 meson pressure 并入后，派生量仍走统一 AD 热力学主线。

## 使用建议

- 若先做文档后做代码：优先使用 `shared/` 文档统一符号，再按模型主题实现。
- 若做 API 文档：将稳定入口同步到 `docs/api/models/*`，本目录保留理论与公式层。
- 若发现公式与代码不一致：以代码为准修正文档，并记录偏差原因与修订日期。

## 文献对齐记录

- `rotation/Rotation_PNJL_CoreEquations.md` 已按 arXiv:2307.14402v1 做了一轮公式口径对齐（含 `\epsilon_n`、Bessel 权重、Polyakov 势参数化与驻点方程）。
- `gas_liquid/GasLiquid_RMF_CoreEquations.md` 已按 `thesis.pdf` 第 2.2 节做逐式对照（Eq. 2.40-2.70）与符号映射；对提取不确定项已显式标注“待复核”。
- `gas_liquid/GasLiquid_RMF_CoreEquations.md` 已补“基于当前实现反推的不确定项”，用于指导后续参数配置收敛（区分占位实现与论文定式）。
- `gas_liquid/GasLiquid_RMF_CoreEquations.md` 已吸收 `tmp/公式文档temp.md` 的可验证信息，并把耦合比口径统一为原文符号 `f_sigma/f_omega/f_rho`。
