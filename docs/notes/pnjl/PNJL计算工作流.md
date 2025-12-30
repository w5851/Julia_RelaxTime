# PNJL 模型计算工作流

> 文档版本：2025-12-30
> 状态：已实现

---

## 一、工作流概述

PNJL 模型的完整计算流程需要按特定顺序执行，因为后续步骤依赖前序步骤的输出数据。

```
┌─────────────────────────────────────────────────────────────────┐
│                      PNJL 计算工作流                             │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  Step 1: 相结构计算 ✅                                           │
│  ├── T-ρ 扫描                                                   │
│  ├── CEP 搜索 (临界终点，二分法)                                  │
│  ├── Maxwell 等面积构造                                          │
│  ├── Spinodal 检测                                              │
│  └── 输出: cep.csv, boundary.csv, spinodals.csv                 │
│           ↓                                                     │
│  Step 2: 相图绘制 ✅                                             │
│  ├── T-μ 相图（含 spinodal 线）                                  │
│  ├── T-ρ 相图（含 spinodal 线）                                  │
│  └── 输出: phase_diagram_*.png                                  │
│           ↓                                                     │
│  Step 3: 参数空间扫描                                            │
│  ├── T-μ 扫描 (使用 PhaseAwareContinuitySeed)                   │
│  ├── T-ρ 扫描                                                   │
│  └── 输出: tmu_scan.csv, trho_scan.csv                          │
│           ↓                                                     │
│  Step 4: 热力学导数计算                                          │
│  ├── 质量导数 (∂M/∂T, ∂M/∂μ)                                    │
│  ├── 热力学响应函数                                              │
│  ├── 体粘滞系数                                                  │
│  └── 输出: derivatives.csv, bulk_viscosity.csv                  │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

---

## 二、已实现的脚本

### 2.1 相结构计算

**脚本**: `scripts/pnjl/calculate_phase_structure.jl`

**功能**:
1. T-ρ 扫描（反向扫描避免 ρ=0 奇异点）
2. CEP 搜索（二分法，精度 0.01 MeV）
3. Maxwell 等面积构造
4. Spinodal 检测（二次插值细化）

**用法**:
```bash
julia scripts/pnjl/calculate_phase_structure.jl --xi=0.0 --T_min=50 --T_max=150 --T_step=5
julia scripts/pnjl/calculate_phase_structure.jl --xi=0.2 --T_min=50 --T_max=140 --T_step=5
julia scripts/pnjl/calculate_phase_structure.jl --xi=0.4 --T_min=50 --T_max=130 --T_step=5
```

**输出文件** (`data/reference/pnjl/`):

| 文件 | 内容 | 列 |
|------|------|-----|
| `cep.csv` | 临界终点 | xi, T_CEP_MeV, mu_CEP_MeV |
| `boundary.csv` | 一阶相变线 | xi, T_MeV, mu_transition_MeV, rho_hadron, rho_quark |
| `spinodals.csv` | Spinodal 线 | xi, T_MeV, mu_spinodal_hadron_MeV, mu_spinodal_quark_MeV, rho_spinodal_hadron, rho_spinodal_quark |

### 2.2 相图绘制

**脚本**: `scripts/pnjl/plot_phase_diagram.py`

**功能**:
- T-μ 相图（一阶相变线 + CEP + spinodal 虚线）
- T-ρ 相图（共存区 + spinodal 虚线）
- 组合相图

**用法**:
```bash
# 绘制所有 ξ 值
python scripts/pnjl/plot_phase_diagram.py --no-show --type all

# 绘制单个 ξ 值
python scripts/pnjl/plot_phase_diagram.py --no-show --type T-mu --xi 0.0

# 不绘制 spinodal 线
python scripts/pnjl/plot_phase_diagram.py --no-show --no-spinodal
```

**输出文件** (`data/outputs/figures/pnjl/`):
- `phase_diagram_T_mu.png`
- `phase_diagram_T_rho.png`
- `phase_diagram_combined.png`

### 2.3 数据验证

**脚本**: `scripts/pnjl/validate_phase_data.py`

**功能**:
- 检测数据跳变（光滑性检查）
- 检测单调性异常
- 输出统计信息和问题报告

**用法**:
```bash
python scripts/pnjl/validate_phase_data.py
```

---

## 三、核心模块

### 3.1 相变分析模块

**位置**: `src/pnjl/analysis/PhaseTransition.jl`

**导出函数**:
- `detect_s_shape(mu_vals, rho_vals)`: 检测 μ(ρ) 曲线的 S 形特征
- `maxwell_construction(mu_vals, rho_vals)`: Maxwell 等面积构造
- `group_curves_by_temperature(rows)`: 按温度分组数据

**数据结构**:
```julia
struct SShapeResult
    has_s_shape::Bool
    mu_spinodal_hadron::Union{Nothing, Float64}  # μ(ρ) 极大值
    mu_spinodal_quark::Union{Nothing, Float64}   # μ(ρ) 极小值
    rho_spinodal_hadron::Union{Nothing, Float64} # 强子相侧 spinodal
    rho_spinodal_quark::Union{Nothing, Float64}  # 夸克相侧 spinodal
    derivative_sign_changes::Int
end

struct MaxwellResult
    converged::Bool
    mu_transition::Union{Nothing, Float64}  # 相变化学势
    rho_hadron::Union{Nothing, Float64}     # 强子相密度
    rho_quark::Union{Nothing, Float64}      # 夸克相密度
    area_residual::Union{Nothing, Float64}
    iterations::Int
    details::Dict{Symbol, Any}
end
```

### 3.2 初值策略模块

**位置**: `src/pnjl/solver/SeedStrategies.jl`

**相关策略**:
- `PhaseAwareSeed`: 单点求解，根据相边界自动选择初值
- `PhaseAwareContinuitySeed`: 连续扫描，结合相边界感知和连续性追踪

---

## 四、数据格式说明

### 4.1 boundary.csv 列说明

| 列名 | 含义 | 单位 |
|------|------|------|
| xi | 各向异性参数 | - |
| T_MeV | 温度 | MeV |
| mu_transition_MeV | 相变化学势（Maxwell 构造） | MeV |
| rho_hadron | 强子相密度（低密度侧） | ρ/ρ₀ |
| rho_quark | 夸克相密度（高密度侧） | ρ/ρ₀ |

### 4.2 spinodals.csv 列说明

| 列名 | 含义 | 单位 |
|------|------|------|
| xi | 各向异性参数 | - |
| T_MeV | 温度 | MeV |
| mu_spinodal_hadron_MeV | μ(ρ) 曲线极大值（强子相侧 spinodal） | MeV |
| mu_spinodal_quark_MeV | μ(ρ) 曲线极小值（夸克相侧 spinodal） | MeV |
| rho_spinodal_hadron | 强子相侧 spinodal 的密度 | ρ/ρ₀ |
| rho_spinodal_quark | 夸克相侧 spinodal 的密度 | ρ/ρ₀ |

### 4.3 物理意义

```
μ(ρ) 曲线示意图（一阶相变区域）:

    μ
    │     ╭──╮ ← spinodal_hadron (极大值)
    │    ╱    ╲
    │   ╱      ╲
    │──╱────────╲── μ_transition (Maxwell 构造)
    │ ╱          ╲
    │╱            ╰──╮ ← spinodal_quark (极小值)
    │                 ╲
    └─────────────────────→ ρ
       ↑              ↑
    rho_hadron    rho_quark
```

- **spinodal_hadron**: 强子相的亚稳态边界，超过此点强子相不稳定
- **spinodal_quark**: 夸克相的亚稳态边界，低于此点夸克相不稳定
- **mu_transition**: 两相共存的化学势（等面积法则）

---

## 五、当前数据状态

已生成的数据（`data/reference/pnjl/`）:

| ξ | T 范围 | boundary 点数 | spinodal 点数 | CEP |
|---|--------|--------------|--------------|-----|
| 0.0 | 50-131 MeV | 25 | 19 | T=130.94, μ=291.21 |
| 0.2 | 50-123 MeV | 21 | 20 | T=123.00, μ=316.06 |
| 0.4 | 50-113 MeV | 19 | 18 | T=113.02, μ=338.18 |

---

## 六、后续工作

### 待实现
- [ ] Crossover 线计算
- [ ] T-μ 扫描脚本（使用 PhaseAwareContinuitySeed）
- [ ] 热力学导数计算脚本
- [ ] GitHub Actions 自动化

### 可能的改进
- [ ] 将脚本中的重复代码迁移到模块
- [ ] 添加更多 ξ 值的数据
- [ ] CEP 附近的细化扫描

---

## 七、相关文档

- `docs/notes/pnjl/一阶相变方法对比.md` - 方法对比
- `docs/notes/pnjl/一阶相变双分支扫描策略.md` - 双分支策略
- `docs/dev/PNJL_Solver_Refactoring_Notes.md` - 重构笔记

---

*更新日期：2025-12-30*
