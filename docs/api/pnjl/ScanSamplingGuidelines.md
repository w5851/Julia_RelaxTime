# PNJL 扫描采样模板与禁用区间

本文档给出 `run_tmu_scan` / `run_trho_scan` 的推荐采样模板与禁用区间，用于降低不收敛、分支跳变与相图伪信号风险。

## 1. 适用范围

- T-μ 扫描：`PNJL.run_tmu_scan`
- T-ρ 扫描：`PNJL.run_trho_scan`
- 相结构流程：`scripts/pnjl/calculate_phase_structure.jl`

---

## 2. 推荐采样模板

### 2.1 快速冒烟模板（smoke）

用途：流程验证、接口联调、CI 级快速检查。

- T-μ：
  - `T = 100:20:140 MeV`
  - `μ = 0:60:240 MeV`
  - `p_num=12`, `t_num=4`
  - `use_phase_aware=true`
- T-ρ：
  - `T = 80:20:160 MeV`
  - `ρ = 0.0:0.2:2.0`
  - `reverse_rho=true`

与模板脚本 `scripts/pnjl/run_aniso_phase_template.jl --profile=smoke` 一致。

### 2.2 常规研究模板（standard）

用途：常规参数扫描与趋势分析。

- T-μ：
  - `T = 50:10:200 MeV`
  - `μ = 0:10:400 MeV`
  - `p_num=24`, `t_num=8`
  - `use_phase_aware=true`
- T-ρ：
  - 温度同上
  - `ρ = DEFAULT_RHO_VALUES`（分段多分辨率）
  - `reverse_rho=true`

其中 `DEFAULT_RHO_VALUES` 来自 `build_default_rho_grid()`：
- 粗网格：`Δρ=0.05`（全区间）
- 中网格：`Δρ=0.02`（`ρ<=1.0`）
- 细网格：`Δρ=0.01`（`ρ<=0.3`）
- 超细网格：`Δρ=0.005`（`ρ<=0.15`）

### 2.3 相图产出模板（phase-diagram）

用途：边界线/CEP/spinodal/crossover 产出。

- 建议先运行 T-ρ 主扫描（`reverse_rho=true`），再做 Maxwell 与 CEP。
- `calculate_phase_structure.jl` 推荐：
  - `T=30:10:350 MeV`
  - `ρ=0.0:0.05:4.0`（需要时在低密度端加密）
- 若 Maxwell 失败点集中在某温度区间，优先局部降 `Δρ`（而不是盲目扩大全局网格）。

---

## 3. 禁用区间与禁用组合

以下组合默认视为“禁用”或“仅诊断可用，不可作为主线结果”：

1. **T-ρ 扫描含 `ρ=0` 且 `reverse_rho=false`**
   - 原因：低密度端奇异/不稳定更容易破坏连续性跟踪。
   - 处理：使用 `reverse_rho=true`（默认）。

2. **一阶相变敏感区关闭 `phase_aware`**
   - 表现：分支跳变、同一网格路径结果不稳定。
   - 处理：保持 `use_phase_aware=true`，首点保留 multiseed 自举。

3. **相图主线使用过粗 ρ 采样（`Δρ > 0.1`）**
   - 风险：S 形丢失、Maxwell 无法闭合、CEP 误判。
   - 处理：至少使用 `DEFAULT_RHO_VALUES` 或同等分段加密。

4. **同时降低积分分辨率与网格分辨率做物理结论**
   - 例如：`p_num<=8` 且 `t_num<=4` 且粗步长网格。
   - 处理：此配置仅用于 smoke，不用于边界线结论。

---

## 4. 最小验收检查

每次扫描建议至少检查：

1. 输出中 `converged=false` 比例是否异常升高；
2. 同一 `xi` 下边界线是否连续（无明显锯齿跳点）；
3. Maxwell 失败原因是否集中在 `no_s_shape/no_sign_change`；
4. 低密度端（`ρ<0.3`）是否仍有足够采样点支撑插值。

---

## 5. 相关文档

- `docs/api/pnjl/TmuScan.md`
- `docs/api/pnjl/TrhoScan.md`
- `docs/guides/examples/pnjl_aniso_phase_repro_template.md`
- `docs/reference/domain-knowledge/pnjl/trho_seed_chain.md`
- `docs/reference/domain-knowledge/pnjl/maxwell_equal_area_options.md`
