if !isdefined(Main, :PNJL)
    @eval begin
        """
            PNJL

        PNJL 模型主模块，提供统一的接口访问所有子模块功能。

        ## 子模块
        - `Core`: 核心计算（积分、热力学量）
        - `Solver`: 求解器（约束模式、条件函数、初值策略）
        - `Derivatives`: 热力学导数计算
        - `Scans`: T-μ/T-ρ 扫描
        - `Analysis`: 相变分析（Maxwell 构造、S 形检测）

        ## 使用示例
        ```julia
        using PNJL

        # 固定化学势求解
        result = solve(FixedMu(), T_fm, μ_fm)

        # 固定密度求解
        result = solve(FixedRho(1.0), T_fm)

        # 热力学导数
        md = mass_derivatives(T_fm, μ_fm)
        ```
        """
        module PNJL

        const _INCLUDE_ONCE_PATH = normpath(joinpath(@__DIR__, "..", "utils", "IncludeOnce.jl"))
        if !isdefined(Main, :IncludeOnce)
            Base.include(Main, _INCLUDE_ONCE_PATH)
        end
        const IncludeOnce = Main.IncludeOnce

        # 确保常量模块已加载
        const _CONSTANTS_PATH = normpath(joinpath(@__DIR__, "..", "Constants_PNJL.jl"))
        IncludeOnce.include_once!(Main, :Constants_PNJL, _CONSTANTS_PATH)
        using Main.Constants_PNJL

# ============================================================================
# 新架构模块
# ============================================================================

# Core 模块
include(joinpath(@__DIR__, "core", "Integrals.jl"))
include(joinpath(@__DIR__, "core", "Thermodynamics.jl"))

# Solver 模块
include(joinpath(@__DIR__, "solver", "ConstraintModes.jl"))
include(joinpath(@__DIR__, "solver", "SeedStrategies.jl"))
include(joinpath(@__DIR__, "solver", "Conditions.jl"))
include(joinpath(@__DIR__, "solver", "ImplicitSolver.jl"))

# Derivatives 模块
include(joinpath(@__DIR__, "derivatives", "ThermoDerivatives.jl"))

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
export calculate_number_densities
export ρ0

# 导出 Solver 功能
export ConstraintMode, FixedMu, FixedRho, FixedAsymmetricRho, FixedEntropy, FixedSigma
export state_dim, param_dim, constraint_description
export SeedStrategy, DefaultSeed, MultiSeed, ContinuitySeed, HybridContinuitySeed, PhaseAwareSeed, PhaseAwareContinuitySeed
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

include(joinpath(@__DIR__, "scans", "ScanCommon.jl"))
include(joinpath(@__DIR__, "scans", "ScanResultFinalize.jl"))
include(joinpath(@__DIR__, "scans", "TmuScan.jl"))
include(joinpath(@__DIR__, "scans", "TrhoScan.jl"))
include(joinpath(@__DIR__, "scans", "AdaptiveRhoRefinement.jl"))

using .TmuScan
using .TrhoScan
using .AdaptiveRhoRefinement

"""load_dual_branch_scan!() -> Module

按需加载一阶相变分析用的双分支扫描模块 `PNJL.DualBranchScan`。

说明：该模块不进入主线默认导出/默认加载；需要时显式调用本函数。
"""
function load_dual_branch_scan!()
    if !isdefined(@__MODULE__, :DualBranchScan)
        include(joinpath(@__DIR__, "scans", "DualBranchScan.jl"))
    end
    return DualBranchScan
end

export load_dual_branch_scan!

# 导出扫描功能
export run_tmu_scan, run_trho_scan
export build_default_rho_grid

# 自适应 ρ 网格加密
export AdaptiveRhoConfig, suggest_refinement_points, merge_rho_values

# ============================================================================
# 分析模块
# ============================================================================

include(joinpath(@__DIR__, "analysis", "PhaseTransition.jl"))

using .PhaseTransition

# 导出分析功能
export SShapeResult, detect_s_shape
export MaxwellResult, maxwell_construction
export group_curves_by_temperature
export CrossoverResult, detect_crossover, scan_crossover_line

# ============================================================================
# 兼容性模块（旧接口）- 已弃用
# ============================================================================

# 注意：以下旧模块已删除，功能已整合到新架构
# - AnisoGapSolver.jl（已弃用，功能在 solver/ImplicitSolver.jl）
# - SeedCache.jl（已弃用）
# - TrhoSeedChain.jl（依赖 LineSearches）
# - SinglePointSolver.jl
# - AdaptiveRhoRefinement.jl
# - CEPFinder.jl, MaxwellRhoMu.jl（已整合到 PhaseTransition.jl）
# - analysis/ThermoDerivatives.jl（旧版，新版在 derivatives/）

        end # module PNJL
    end
end
