---
title: models Ω 公式一致性核查草表（本地草案）
archived: true
original: docs/dev/active/2026-04-20_models_omega_formula_audit_draft.md
archived_date: 2026-04-20
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# models Ω 公式一致性核查草表（本地草案）

> 日期：2026-04-20
> 目的：先给 issue #92 的“Formula Fidelity Audit + Ω 复用矩阵”提供可分派底稿。
> 说明：这里是草案；本文先固化 #93 首轮（PNJL / rPNJL / PNJLMagnetic）核查结果，后续再扩展到 NJL/NJL2 与变体模型。
> 状态：核心结论已被实现覆盖，见 commits `d7adcec` / `d96fd64`。

---

## 1) 审计口径

- 公式来源目录：`docs/reference/formula/models/**`
- 代码来源目录：`src/models/**`
- 关注对象：`Ω` 组装与组分项（`chi/vac/therm/poly`）及 `P=-Ω` 派生链路

状态定义：

- `match`：与公式文档主表达基本一致
- `partial`：主结构一致，但有实现级近似/简化
- `mismatch`：关键定义不一致
- `placeholder`：明显占位/工程近似，未达公式复现

---

## 2) Ω-term 复用矩阵（草案）

| 模型 | 公式文档 | Ω 主通道 | chi | vac | therm | poly | 初判 |
|---|---|---|---|---|---|---|---|
| NJL | `docs/reference/formula/models/njl/NJL_core.md` | 通用 `omega.jl` | 模型实现 | 模型实现 | 模型实现 | 固定为 0 | match |
| NJL2 | `docs/reference/formula/models/njl/NJL_core.md`（两味补充） | 通用 `omega.jl` | 模型实现 | 模型实现 | 模型实现 | 固定为 0 | match |
| PNJL | `docs/reference/formula/models/pnjl/Omega_各向同性.md` | 通用 `omega.jl` | 模型实现 | 模型实现 | 模型实现 | 模型实现 | match |
| rPNJL | `docs/reference/formula/models/rpnjl/rPNJL_core.md` | 通用 `omega.jl` | 模型实现（含八夸克） | 模型实现 | 复用 PNJL 热项路径 | 模型实现（含 Vandermonde） | match/partial |
| PNJLMagnetic | `docs/reference/formula/models/pnjl_magnetic/PNJL_magnetic_core.md` | 模型特化 `omega_components` | 专用模块（含 `G(B)`） | 专用模块（Landau） | 专用模块（Landau） | 复用 PNJL Polyakov | partial |
| Rotation | `docs/reference/formula/models/rotation/Rotation_PNJL_CoreEquations.md` | 模型特化 `omega_components` | 简化实现 | 简化实现 | 特化核 | 特化核 | placeholder/partial |
| GasLiquid | `docs/reference/formula/models/gas_liquid/GasLiquid_RMF_CoreEquations.md` | 模型特化 `omega_components` | 近似/映射 | 0 占位 | 近似反推 | 0 占位 | placeholder |

---

## 3) 关键证据索引（代码）

- 通用 Ω 组装：`src/models/omega.jl:21`
- 通用热力学派生：`src/models/thermo_kernel.jl:24`
- 抽象 term 契约：`src/models/abstract_model.jl:90`

- PNJL term 实现：`src/models/pnjl_physics/PNJLModel.jl:47`, `src/models/pnjl_physics/PNJLModel.jl:55`, `src/models/pnjl_physics/PNJLModel.jl:59`, `src/models/pnjl_physics/PNJLModel.jl:70`
- rPNJL term 实现：`src/models/rpnjl/RPNJLModel.jl:202`, `src/models/rpnjl/RPNJLModel.jl:229`, `src/models/rpnjl/RPNJLModel.jl:249`, `src/models/rpnjl/RPNJLModel.jl:294`
- PNJLMagnetic 特化通道：`src/models/pnjl_physics/PNJLMagneticModel.jl:61`, `src/models/pnjl_physics/core/MagneticThermodynamics.jl:124`
- NJL term 实现：`src/models/njl/NJLModel.jl:97`, `src/models/njl/NJLModel.jl:101`, `src/models/njl/NJLModel.jl:109`, `src/models/njl/NJLModel.jl:123`
- NJL2 term 实现：`src/models/njl/NJL2Model.jl:82`, `src/models/njl/NJL2Model.jl:87`, `src/models/njl/NJL2Model.jl:97`, `src/models/njl/NJL2Model.jl:111`

- Rotation 特化与简化项：`src/models/variants/rotation/RotationModel.jl:64`, `src/models/variants/rotation/RotationModel.jl:71`, `src/models/variants/rotation/RotationModel.jl:96`
- GasLiquid 特化与占位项：`src/models/variants/gas_liquid/GasLiquidModel.jl:64`, `src/models/variants/gas_liquid/GasLiquidModel.jl:69`, `src/models/variants/gas_liquid/GasLiquidModel.jl:86`

---

## 4) 关键证据索引（公式文档）

- PNJL Ω：`docs/reference/formula/models/pnjl/Omega_各向同性.md:13`
- rPNJL Ω：`docs/reference/formula/models/rpnjl/rPNJL_core.md:69`
- Rotation Ω：`docs/reference/formula/models/rotation/Rotation_PNJL_CoreEquations.md:31`
- GasLiquid Ω：`docs/reference/formula/models/gas_liquid/GasLiquid_RMF_CoreEquations.md:96`
- GasLiquid 文档中已记录占位风险：`docs/reference/formula/models/gas_liquid/GasLiquid_RMF_CoreEquations.md:243`
- NJL Ω：`docs/reference/formula/models/njl/NJL_core.md:270`

---

## 5) “只改 poly、复用其余项”可行性（草案）

可行条件：

1. 模型走通用 `omega_components(model, ...)` 主线（`src/models/omega.jl:21`）
2. `chi/vac/therm` 不覆盖（或覆盖为通用实现）
3. 仅覆写 `polyakov_potential(model, Φ, Φbar, T; ...)`

结论（当前代码态）：

- PNJL/rPNJL 路径可直接支持这种复用模式。
- Rotation/GasLiquid 当前因为覆盖了整套 `omega_components`，暂不满足“仅改 poly”最小覆写路径。

---

## 6) #93 首轮逐项核查（PNJL / rPNJL / PNJLMagnetic）

### 6.1 PNJL（首轮结论：`match`）

- 公式口径：
  - `Ω = χ + U + Ω_vac + Ω_therm`（`docs/reference/formula/models/pnjl/Omega_各向同性.md:13`）
  - `P = -Ω`（`docs/reference/formula/models/pnjl/Omega_各向同性.md:159`）
- 代码映射：
  - 统一组装 `χ + poly + vac + therm`：`src/models/omega.jl:21`
  - `χ`：`calculate_chiral` -> `PNJLCore.chiral_potential`：`src/models/pnjl_physics/PNJLModel.jl:51`
  - `U`：`polyakov_potential`：`src/models/pnjl_physics/PNJLModel.jl:55`
  - `Ω_vac`：三味真空积分加和：`src/models/pnjl_physics/PNJLModel.jl:59`
  - `Ω_therm`：PNJL 热项积分核：`src/models/pnjl_physics/PNJLModel.jl:70`
- 判定：主结构与项级映射一致，首轮按 `match` 处理。

### 6.2 rPNJL（首轮结论：`match/partial`）

- 公式口径：
  - 在 PNJL 主体上加入八夸克项（`g1,g2`）与 Vandermonde 项（`kappa*log J`）：`docs/reference/formula/models/rpnjl/rPNJL_core.md:69`
- 代码映射：
  - 质量方程中加入 `g1/g2` 修正：`src/models/rpnjl/RPNJLModel.jl:202`
  - 手征势中加入八夸克项：`src/models/rpnjl/RPNJLModel.jl:229`
  - Polyakov 势加入 `-kappa*log(jac)`：`src/models/rpnjl/RPNJLModel.jl:249`
  - 热项复用 PNJL 路径：`src/models/rpnjl/RPNJLModel.jl:294`
- 判定：结构匹配；`partial` 原因是默认配置 `g1/g2/kappa` 可退化到 0（退化回 PNJL 极限），需在回归集固定“rPNJL 非零参数”样本点验证。

### 6.3 PNJLMagnetic（首轮结论：`partial`）

- 公式口径：
  - Landau 量子化真空项/热项与 `E_{f,n}`：`docs/reference/formula/models/pnjl_magnetic/PNJL_magnetic_core.md:37`, `docs/reference/formula/models/pnjl_magnetic/PNJL_magnetic_core.md:47`, `docs/reference/formula/models/pnjl_magnetic/PNJL_magnetic_core.md:66`, `docs/reference/formula/models/pnjl_magnetic/PNJL_magnetic_core.md:84`
  - 磁场依赖耦合 `G(B)`：`docs/reference/formula/models/pnjl_magnetic/PNJL_magnetic_core.md:126`
- 代码映射：
  - 专用 `omega_components` 入口：`src/models/pnjl_physics/PNJLMagneticModel.jl:61`
  - Landau 求和真空/热项：`src/models/pnjl_physics/core/MagneticThermodynamics.jl:157`
  - `G(B)` 参数化：`src/models/pnjl_physics/core/MagneticThermodynamics.jl:88`
- 判定：`eB != 0` 路径与文档主结构一致；`partial` 原因是 `eB -> 0` 分支当前将 `vac/therm` 归并为总和输出（组件字段置零），项级可追溯性较弱，但总 Ω 口径保持一致。

### 6.4 NJL / NJL2（增补结论：`match`）

- 公式口径：
  - NJL `Ω = chiral + vacuum + thermal`，并满足 `P=-Ω`（`docs/reference/formula/models/njl/NJL_core.md:270`）。
  - NJL2 采用两味退化口径（同分项结构，无 Polyakov 项）。
- 代码映射：
  - 通用组装 `chi + poly + vac + therm`：`src/models/omega.jl:21`
  - NJL：
    - `chi`：`src/models/njl/NJLModel.jl:101`
    - `poly=0`：`src/models/njl/NJLModel.jl:105`
    - `vac`：`src/models/njl/NJLModel.jl:109`
    - `therm`：`src/models/njl/NJLModel.jl:123`
  - NJL2：
    - `chi`：`src/models/njl/NJL2Model.jl:87`
    - `poly=0`：`src/models/njl/NJL2Model.jl:92`
    - `vac`：`src/models/njl/NJL2Model.jl:97`
    - `therm`：`src/models/njl/NJL2Model.jl:111`
- 判定：两者均可直接映射到通用四项分解，首轮按 `match` 处理。

## 7) #93 下一步（可直接分派）

- [ ] 任务 A：将 PNJL/rPNJL/PNJLMagnetic 首轮结论回填到 issue #93（含证据链接）
- [ ] 任务 B：在回归层新增/确认 rPNJL 非零 `g1/g2/kappa` 固定点样本，锁定 `match/partial` 中的 `partial` 风险项
- [ ] 任务 C：为 `PNJLMagnetic` 增加 `eB -> 0` 项级可追溯注记（避免仅总和可追溯）
- [x] 任务 D：继续扩展到 NJL/NJL2，并形成全模型同口径矩阵（首轮完成，后续补固定点验证）

## 8) 全模型可机读矩阵（首轮）

| model | formula_doc | omega_path | chi | vac | therm | poly | status | evidence_code | evidence_formula |
|---|---|---|---|---|---|---|---|---|---|
| NJL | `docs/reference/formula/models/njl/NJL_core.md` | generic | implemented | implemented | implemented | zero | match | `src/models/njl/NJLModel.jl:101`; `src/models/njl/NJLModel.jl:109`; `src/models/njl/NJLModel.jl:123`; `src/models/njl/NJLModel.jl:105` | `docs/reference/formula/models/njl/NJL_core.md:270` |
| NJL2 | `docs/reference/formula/models/njl/NJL_core.md` | generic | implemented | implemented | implemented | zero | match | `src/models/njl/NJL2Model.jl:87`; `src/models/njl/NJL2Model.jl:97`; `src/models/njl/NJL2Model.jl:111`; `src/models/njl/NJL2Model.jl:92` | `docs/reference/formula/models/njl/NJL_core.md:407` |
| PNJL | `docs/reference/formula/models/pnjl/Omega_各向同性.md` | generic | implemented | implemented | implemented | implemented | match | `src/models/pnjl_physics/PNJLModel.jl:51`; `src/models/pnjl_physics/PNJLModel.jl:59`; `src/models/pnjl_physics/PNJLModel.jl:70`; `src/models/pnjl_physics/PNJLModel.jl:55` | `docs/reference/formula/models/pnjl/Omega_各向同性.md:13` |
| rPNJL | `docs/reference/formula/models/rpnjl/rPNJL_core.md` | generic | implemented_ext | implemented | implemented_reuse | implemented_ext | match/partial | `src/models/rpnjl/RPNJLModel.jl:229`; `src/models/rpnjl/RPNJLModel.jl:282`; `src/models/rpnjl/RPNJLModel.jl:294`; `src/models/rpnjl/RPNJLModel.jl:249` | `docs/reference/formula/models/rpnjl/rPNJL_core.md:69` |
| PNJLMagnetic | `docs/reference/formula/models/pnjl_magnetic/PNJL_magnetic_core.md` | specialized | implemented_ext | implemented_landau | implemented_landau | implemented_reuse | partial | `src/models/pnjl_physics/PNJLMagneticModel.jl:61`; `src/models/pnjl_physics/core/MagneticThermodynamics.jl:157`; `src/models/pnjl_physics/core/MagneticThermodynamics.jl:167`; `src/models/pnjl_physics/core/MagneticThermodynamics.jl:182` | `docs/reference/formula/models/pnjl_magnetic/PNJL_magnetic_core.md:37` |
| Rotation | `docs/reference/formula/models/rotation/Rotation_PNJL_CoreEquations.md` | specialized | implemented/placeholder | placeholder | implemented_core | implemented | placeholder/partial | `src/models/variants/rotation/RotationModel.jl:96`; `src/models/variants/rotation/core/RotationThermo.jl:154`; `src/models/variants/rotation/RotationModel.jl:71`; `src/models/variants/rotation/core/RotationThermo.jl:156` | `docs/reference/formula/models/rotation/Rotation_PNJL_CoreEquations.md:36` |
| GasLiquid | `docs/reference/formula/models/gas_liquid/GasLiquid_RMF_CoreEquations.md` | specialized | implemented_core/placeholder_top | implemented_core/placeholder_top | implemented_core/placeholder_top | zero | placeholder/partial | `src/models/variants/gas_liquid/core/Thermodynamics.jl:56`; `src/models/variants/gas_liquid/core/Thermodynamics.jl:57`; `src/models/variants/gas_liquid/core/Thermodynamics.jl:58`; `src/models/variants/gas_liquid/GasLiquidModel.jl:54` | `docs/reference/formula/models/gas_liquid/GasLiquid_RMF_CoreEquations.md:101` |

字段说明：

- `omega_path`：`generic` 表示走 `src/models/omega.jl` 通用组装；`specialized` 表示模型自定义 `omega_components`。
- `implemented_ext`：在基础项上包含扩展物理项（如 rPNJL 八夸克、Vandermonde）。
- `placeholder_top`：模型顶层接口仍含占位实现；`implemented_core` 表示核心热力学模块已有分项实现。
