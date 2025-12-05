# 计划：PNJL 相图扫描与分析功能实现

基于现有的 `AnisoGapSolver` 求解器和种子系统，实现 T-μ / T-ρ 网格扫描，并通过曲线分析（S 形检测、Maxwell 构造）定位临界点和一阶相变线，最终绘制完整相图。

## 步骤总览

| 步骤 | 内容 | 输出 |
|------|------|------|
| 1 | T-μ 网格扫描模块 | `src/pnjl/scans/TmuScan.jl` + CSV |
| 2 | T-ρ 网格扫描模块 | `src/pnjl/scans/TrhoScan.jl` + CSV |
| 3 | 绘图工具（ρ-μ / P-μ 曲线） | `scripts/pnjl/plot_curves.py` 或 `.jl` |
| 4 | S 形检测与稳态判断 | `src/pnjl/analysis/PhaseTransition.jl` |
| 5 | 临界点检测 | `src/pnjl/analysis/CEPFinder.jl` |
| 6 | 一阶相变线判断（Maxwell 构造） | 扩展 `PhaseTransition.jl` |
| 7 | 完整相图绘制 | `scripts/pnjl/plot_phase_diagram.py` |

---

## 步骤 1：T-μ 网格扫描模块

**目标**：遍历 (T, μ) 网格，调用 `solve_fixed_mu`，输出热力学量到 CSV。

**文件**：`src/pnjl/scans/TmuScan.jl`

**接口设计**：
```julia
function run_tmu_scan(;
    T_range = 50.0:10.0:200.0,   # MeV
    mu_range = 0.0:10.0:400.0,   # MeV
    xi = 0.0,
    output_path = "data/outputs/results/pnjl/tmu_scan.csv",
    seed_path = DEFAULT_SEED_PATH,
    p_num = 24, t_num = 8
)
```

**CSV 列**：`T_MeV, mu_MeV, xi, pressure, rho, entropy, energy, phi_u, phi_d, phi_s, Phi1, Phi2, converged, iterations`

**要点**：
- 使用 `SeedCache.find_initial_seed` 获取暖启动
- 支持断点续扫（检查已存在的行，跳过已计算点）
- 失败点记录但不中断

---

## 步骤 2：T-ρ 网格扫描模块

**目标**：遍历 (T, ρ) 网格，调用 `solve_fixed_rho`，输出到 CSV。

**文件**：`src/pnjl/scans/TrhoScan.jl`

**接口设计**：
```julia
function run_trho_scan(;
    T_range = 50.0:10.0:200.0,   # MeV
    rho_range = 0.0:0.05:3.0,    # ρ/ρ₀
    xi = 0.0,
    output_path = "data/outputs/results/pnjl/trho_scan.csv",
    linesearch = LineSearches.BackTracking(),  # 现已作为默认设置，可覆盖
    ...
)
```

**CSV 列**：`T_MeV, rho, xi, mu_u_MeV, mu_d_MeV, mu_s_MeV, mu_avg_MeV, pressure, entropy, energy, phi_u, phi_d, phi_s, Phi1, Phi2, converged`

**要点**：
- 使用 BackTracking 线搜索提升 fixed-ρ 收敛性
- 对于多解区域（S 形），从两端分别扫描以捕获不同分支

---

## 步骤 3：绘图工具（ρ-μ / P-μ 曲线）

**目标**：读取扫描 CSV，绘制固定温度下的 ρ(μ) 和 P(μ) 曲线。

**文件**：`scripts/pnjl/plot_curves.py`（推荐 Python + matplotlib，与 PNJL_Simulation 保持一致）

**功能**：
- `plot_rho_mu(csv_path, T_list)`: 绘制多条 ρ-μ 曲线（不同温度）
- `plot_P_mu(csv_path, T_list)`: 绘制多条 P-μ 曲线
- 自动标注 S 形区域和相变点

---

## 步骤 4：S 形检测与稳态判断

**目标**：分析 T-ρ 扫描结果，检测 ρ(μ) 曲线的 S 形（非单调性），判断稳态/亚稳态/不稳定区。

**文件**：`src/pnjl/analysis/PhaseTransition.jl`

**算法**：
1. 对固定 T 的 ρ-μ 数据，计算 dρ/dμ
2. 若存在 dρ/dμ < 0 区域 → S 形存在 → 一阶相变
3. S 形的两个极值点对应稳态与亚稳态的边界（spinodal 线）
4. 中间区域为热力学不稳定区

**接口**：
```julia
struct SShapeResult
    has_s_shape::Bool
    mu_spinodal_low::Float64   # 下 spinodal 点
    mu_spinodal_high::Float64  # 上 spinodal 点
    rho_at_spinodals::Tuple{Float64, Float64}
end

function detect_s_shape(T_mev, trho_data) -> SShapeResult
```

---

## 步骤 5：临界点检测

**目标**：扫描温度，找到 S 形消失的临界点 (T_CEP, μ_CEP)。

**文件**：`src/pnjl/analysis/CEPFinder.jl`

**算法**：
1. 从低温开始，逐步升温
2. 对每个 T 调用 `detect_s_shape`
3. S 形存在 → 继续升温；S 形消失 → 二分法细化 T
4. CEP 处：S 形刚好消失，dρ/dμ 有拐点但无极值

**接口**：
```julia
struct CEPResult
    T_cep_MeV::Float64
    mu_cep_MeV::Float64
    rho_cep::Float64
    confidence::Symbol  # :high, :medium, :low
end

function find_cep(trho_data; T_range, tol=1.0) -> CEPResult
```

---

## 步骤 6：一阶相变线判断（Maxwell 构造）

**目标**：从 P-μ 曲线确定一阶相变的化学势（两相共存点）。

**扩展**：`src/pnjl/analysis/PhaseTransition.jl`

**算法（Maxwell 等面积法则）**：
1. 对固定 T 的 P(μ) 曲线，若存在回环（S 形对应的多值区）
2. 找到 μ_coex 使得：∫[μ_low→μ_coex] (P - P_coex) dμ = ∫[μ_coex→μ_high] (P_coex - P) dμ
3. 或等价地：找 P_coex 使上下两块面积相等

**接口**：
```julia
struct MaxwellResult
    mu_coex_MeV::Float64      # 相变化学势
    P_coex::Float64           # 共存压力
    rho_gas::Float64          # 气相密度
    rho_liquid::Float64       # 液相密度
end

function maxwell_construction(T_mev, tmu_data) -> MaxwellResult
```

---

## 步骤 7：完整相图绘制

**目标**：在 T-μ 平面绘制：
- 一阶相变线（T < T_CEP 时的 μ_coex(T)）
- CEP 标记
- Crossover 区域（T > T_CEP）
- 可选：spinodal 线（亚稳态边界）

**文件**：`scripts/pnjl/plot_phase_diagram.py`

**输入**：
- 步骤 5 的 CEP 结果
- 步骤 6 对多个温度的 Maxwell 结果（形成一阶相变线）

**输出**：`data/outputs/figures/pnjl/pnjl_phase_diagram.png`

---

## 依赖关系图

```
步骤 1 (TmuScan) ──────┬──→ 步骤 3 (绘图 P-μ) ──→ 步骤 6 (Maxwell) ──┐
                       │                                              │
                       └──────────────────────────────────────────────┼──→ 步骤 7 (相图)
                                                                       │
步骤 2 (TrhoScan) ──→ 步骤 3 (绘图 ρ-μ) ──→ 步骤 4 (S形检测) ──→ 步骤 5 (CEP) ──┘
```

---

## 进一步考虑

1. **各向异性扩展**：所有模块支持 `xi` 参数，可生成不同 ξ 下的相图族
2. **并行化**：扫描模块可用 `@threads` 或 `Distributed` 加速
3. **服务器集成**：扫描可封装为长任务 API（`/scan/tmu`, `/scan/trho`）
4. **与 PNJL_Simulation 对齐**：复用其 `calculate_t_rho` 和 `find_cep_T` 逻辑
