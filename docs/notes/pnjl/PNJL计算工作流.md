# PNJL 模型计算工作流

> 文档版本：2025-12-29
> 状态：规划中

---

## 一、工作流概述

PNJL 模型的完整计算流程需要按特定顺序执行，因为后续步骤依赖前序步骤的输出数据。

```
┌─────────────────────────────────────────────────────────────────┐
│                      PNJL 计算工作流                             │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  Step 1: 相结构计算                                              │
│  ├── CEP 搜索 (临界终点)                                         │
│  ├── 一阶相变线 (T-ρ + Maxwell)                                  │
│  ├── Crossover 线                                               │
│  └── 输出: cep.csv, boundary.csv, crossover.csv                 │
│           ↓                                                     │
│  Step 2: 相图绘制                                                │
│  ├── T-μ 相图                                                   │
│  ├── T-ρ 相图                                                   │
│  └── 输出: phase_diagram_Tmu.png, phase_diagram_Trho.png        │
│           ↓                                                     │
│  Step 3: 参数空间扫描                                            │
│  ├── T-μ 扫描 (使用 boundary.csv 选择初解)                       │
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

## 二、各步骤详细说明

### Step 1: 相结构计算

**目的**：确定 QCD 相图的基本结构

**输入参数**：
- `xi`: 各向异性参数（默认 0.0）
- `T_range`: 温度范围（如 50-200 MeV）
- `rho_range`: 密度范围（如 0-4 ρ₀）

**输出文件**：

| 文件 | 内容 | 格式 |
|------|------|------|
| `cep.csv` | 临界终点 (T_CEP, μ_CEP) | T_MeV, mu_MeV, xi |
| `boundary.csv` | 一阶相变线 | T_MeV, mu_coex_MeV, rho_gas, rho_liquid |
| `crossover.csv` | Crossover 线 | T_MeV, mu_crossover_MeV |

**计算方法**：
- CEP：搜索 ∂²Ω/∂φ² = 0 的点
- 一阶相变线：T-ρ 扫描 + Maxwell 等面积构造
- Crossover 线：手征序参量变化最快的位置

**依赖**：无（第一步）

---

### Step 2: 相图绘制

**目的**：可视化相结构

**输入**：Step 1 的输出文件

**输出文件**：
- `phase_diagram_Tmu.png`: T-μ 平面相图
- `phase_diagram_Trho.png`: T-ρ 平面相图

**图中元素**：
- CEP 点（标记）
- 一阶相变线（实线）
- Crossover 线（虚线）
- Spinodal 线（可选，点线）

---

### Step 3: 参数空间扫描

**目的**：生成热力学量数据表

#### 3.1 T-μ 扫描

**输入参数**：
- `T_range`: 温度范围
- `mu_range`: 化学势范围
- `boundary.csv`: 一阶相变线数据（用于初解选择）

**初解选择策略**：
```julia
function select_seed(T, μ, boundary_data)
    μ_c = interpolate_mu_c(boundary_data, T)
    if T > T_CEP || isnan(μ_c)
        return :auto  # Crossover 区域
    end
    return μ < μ_c ? :hadron : :quark
end
```

**输出**：`tmu_scan.csv`
- 列：T_MeV, mu_MeV, xi, omega, pressure, rho, entropy, energy, phi_u, ..., M_u_MeV, ...

#### 3.2 T-ρ 扫描

**输入参数**：
- `T_range`: 温度范围
- `rho_range`: 密度范围

**输出**：`trho_scan.csv`
- 列：T_MeV, rho, xi, mu_MeV, omega, pressure, ...

---

### Step 4: 热力学导数计算

**目的**：计算响应函数和输运系数

**输入**：Step 3 的扫描数据或直接计算

**输出**：
- `derivatives.csv`: 质量导数、热力学响应函数
- `bulk_viscosity.csv`: 体粘滞系数

**计算内容**：
- 质量导数：∂M_f/∂T, ∂M_f/∂μ
- 比热：C_V = T(∂s/∂T)
- 压缩率：κ = -V⁻¹(∂V/∂P)
- 体粘滞系数：ζ/s

---

## 三、数据依赖关系

```
cep.csv ─────────────────────────────────┐
                                         │
boundary.csv ──┬── phase_diagram_*.png   │
               │                         │
               └── tmu_scan.csv ─────────┼── derivatives.csv
                                         │
crossover.csv ───────────────────────────┘
                                         
trho_scan.csv ───────────────────────────── bulk_viscosity.csv
```

**关键依赖**：
- `tmu_scan.csv` 依赖 `boundary.csv`（初解选择）
- `derivatives.csv` 可能依赖扫描数据（数值微分）或直接计算（解析微分）

---

## 四、GitHub Actions 自动化（规划）

### 4.1 设计目标

- 用户可在 GitHub 网站手动触发计算
- 通过 `workflow_dispatch` 输入扫描参数
- 无需本地 Julia 环境
- 结果自动上传为 Artifacts

### 4.2 工作流配置示例

```yaml
name: PNJL Calculation Pipeline

on:
  workflow_dispatch:
    inputs:
      xi:
        description: 'Anisotropy parameter'
        required: true
        default: '0.0'
      T_min:
        description: 'Minimum temperature (MeV)'
        required: true
        default: '50'
      T_max:
        description: 'Maximum temperature (MeV)'
        required: true
        default: '200'
      step:
        description: 'Which step to run (1-4 or all)'
        required: true
        default: 'all'

jobs:
  phase-structure:
    if: ${{ inputs.step == '1' || inputs.step == 'all' }}
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: julia-actions/setup-julia@v2
      - name: Calculate phase structure
        run: julia scripts/calculate_phase_structure.jl ...
      - uses: actions/upload-artifact@v4
        with:
          name: phase-data
          path: data/processed/results/pnjl/

  tmu-scan:
    needs: phase-structure
    if: ${{ inputs.step == '3' || inputs.step == 'all' }}
    runs-on: ubuntu-latest
    steps:
      - uses: actions/download-artifact@v4
        with:
          name: phase-data
      - name: Run T-mu scan
        run: julia scripts/run_tmu_scan.jl ...
```

### 4.3 参数化运行

用户可通过 GitHub Actions UI 输入：
- 各向异性参数 ξ
- 温度范围 [T_min, T_max]
- 化学势/密度范围
- 网格密度
- 要执行的步骤

### 4.4 输出管理

- 计算结果作为 Artifacts 保存
- 可设置保留期限（如 90 天）
- 支持下载或在后续 workflow 中使用

---

## 五、本地运行指南

### 5.1 完整流程

```bash
# Step 1: 相结构
julia scripts/calculate_phase_structure.jl --xi=0.0 --T_range=50:10:200

# Step 2: 相图（可选）
julia scripts/plot_phase_diagram.jl

# Step 3: 扫描
julia scripts/run_tmu_scan.jl --xi=0.0 --T_range=50:10:200 --mu_range=0:10:400
julia scripts/run_trho_scan.jl --xi=0.0 --T_range=50:10:200 --rho_range=0:0.1:4

# Step 4: 导数
julia scripts/calculate_derivatives.jl
```

### 5.2 单步运行

如果只需要更新某一步：
```bash
# 只重新计算 T-μ 扫描（假设 boundary.csv 已存在）
julia scripts/run_tmu_scan.jl --use-existing-boundary
```

---

## 六、相关文件

### 脚本（待创建）
- `scripts/calculate_phase_structure.jl` - 相结构计算
- `scripts/plot_phase_diagram.jl` - 相图绘制
- `scripts/run_tmu_scan.jl` - T-μ 扫描
- `scripts/run_trho_scan.jl` - T-ρ 扫描
- `scripts/calculate_derivatives.jl` - 导数计算

### 数据目录
- `data/processed/results/pnjl/` - 计算结果
- `data/raw/pnjl/` - 原始数据/种子表

### 文档
- `docs/notes/pnjl/一阶相变方法对比.md` - 方法对比
- `docs/notes/pnjl/一阶相变双分支扫描策略.md` - 双分支策略（备用）

---

## 七、注意事项

1. **数据一致性**：修改 ξ 参数后需要重新从 Step 1 开始
2. **CEP 附近**：需要更细的网格（1-2 MeV 步长）
3. **计算时间**：完整流程可能需要数小时，建议使用 GitHub Actions
4. **版本控制**：结果文件不应提交到 git，使用 Artifacts 管理

---

*创建日期：2025-12-29*
