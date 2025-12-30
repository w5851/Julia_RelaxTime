"""
    PNJL

PNJL 模型主模块，提供统一的接口访问所有子模块功能。

## 子模块
- `Core`: 核心计算（积分、热力学量）
- `Solver`: 求解器（约束模式、条件函数、初值策略）
- `Derivatives`: 热力学导数计算

## 兼容性模块（旧接口）
- `AnisoGapSolver`: 旧版求解器（保留兼容性）
- `SeedCache`: 旧版种子缓存（已弃用）

## 使用示例
```julia
using PNJL

# 新接口
result = solve(FixedMu(), T_fm, μ_fm)
md = mass_derivatives(T_fm, μ_fm)

# 旧接口（兼容）
result = AnisoGapSolver.solve_fixed_mu(T_mev, μ_fm)
```
"""

# 确保常量模块已加载
if !isdefined(Main, :Constants_PNJL)
    include(joinpath(@__DIR__, "..", "Constants_PNJL.jl"))
end

module PNJL

using ..Constants_PNJL

# ============================================================================
# 新架构模块
# ============================================================================

# Core 模块
include(joinpath("core", "Integrals.jl"))
include(joinpath("core", "Thermodynamics.jl"))

# Solver 模块
include(joinpath("solver", "ConstraintModes.jl"))
include(joinpath("solver", "SeedStrategies.jl"))
include(joinpath("solver", "Conditions.jl"))
include(joinpath("solver", "ImplicitSolver.jl"))

# Derivatives 模块
include(joinpath("derivatives", "ThermoDerivatives.jl"))

# 使用新模块
using .Integrals
using .Thermodynamics
using .ConstraintModes
using .SeedStrategies
using .Conditions
using .ImplicitSolver
using .ThermoDerivatives

# 导出 Core 功能
export cached_nodes, vacuum_integral, calculate_energy_sum, calculate_log_sum
export DEFAULT_THETA_COUNT, DEFAULT_MOMENTUM_COUNT
export calculate_mass_vec, calculate_chiral, calculate_U
export calculate_pressure, calculate_omega, calculate_rho, calculate_thermo
export ρ0

# 导出 Solver 功能
export ConstraintMode, FixedMu, FixedRho, FixedEntropy, FixedSigma
export state_dim, param_dim, constraint_description
export SeedStrategy, DefaultSeed, MultiSeed, ContinuitySeed, PhaseAwareSeed, PhaseAwareContinuitySeed
export get_seed, update!, reset!, get_all_seeds, set_phase!
export HADRON_SEED_5, QUARK_SEED_5, HADRON_SEED_8, QUARK_SEED_8
export PhaseBoundaryData, load_phase_boundary, interpolate_mu_c, get_phase_hint
export gap_conditions, build_conditions, build_residual!, GapParams
export solve, solve_multi, SolverResult
export create_implicit_solver, solve_with_derivatives

# 导出 Derivatives 功能
export mass_derivatives, thermo_derivatives, bulk_derivative_coeffs
export bulk_viscosity_coefficients, compute_B_bracket
export dP_dT, dP_dmu

# ============================================================================
# 扫描模块（使用新架构）
# ============================================================================

include(joinpath("scans", "TmuScan.jl"))
include(joinpath("scans", "TrhoScan.jl"))
include(joinpath("scans", "DualBranchScan.jl"))

using .TmuScan
using .TrhoScan
using .DualBranchScan

# 导出扫描功能
export run_tmu_scan, run_trho_scan
export build_default_rho_grid
# 双分支扫描
export run_dual_branch_scan, find_phase_transition, merge_branches, scan_phase_diagram
export DualBranchResult, BranchPoint, PhaseTransitionInfo

# ============================================================================
# 分析模块
# ============================================================================

include(joinpath("analysis", "PhaseTransition.jl"))

using .PhaseTransition

# 导出分析功能
export SShapeResult, detect_s_shape
export MaxwellResult, maxwell_construction
export group_curves_by_temperature

# ============================================================================
# 兼容性模块（旧接口）
# ============================================================================

# 旧版求解器（保留兼容性）
include(joinpath("solvers", "AnisoGapSolver.jl"))
using .AnisoGapSolver
export AnisoGapSolver

# 注意：以下旧模块暂时禁用，需要时可手动加载
# - SeedCache.jl（已弃用）
# - TrhoSeedChain.jl（依赖 LineSearches）
# - SinglePointSolver.jl
# - AdaptiveRhoRefinement.jl
# - PhaseTransition.jl, CEPFinder.jl, MaxwellRhoMu.jl
# - analysis/ThermoDerivatives.jl（旧版，新版在 derivatives/）

end # module PNJL
