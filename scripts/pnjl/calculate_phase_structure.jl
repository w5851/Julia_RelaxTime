#!/usr/bin/env julia
"""
PNJL 相结构计算脚本（薄 CLI）

职责：
- 仅做环境激活、可选 warmup、调用共享 phase CLI support 层
- 真正的参数解析与执行逻辑位于 `phase_cli_support.jl`

示例：
    julia scripts/pnjl/calculate_phase_structure.jl --model_kind=PNJL --verbose
    julia scripts/pnjl/calculate_phase_structure.jl --T_min=120 --T_max=180 --T_step=20 --rho_max=0.8 --rho_step=0.2
    julia scripts/pnjl/calculate_phase_structure.jl --promote_reference
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

include(joinpath(@__DIR__, "..", "..", "src", "models", "Models.jl"))
using .Models

include(joinpath(@__DIR__, "phase_cli_support.jl"))
using .PhaseCliSupport

if get(ENV, "PHASE_PRECOMPILE_WARMUP", "1") in ("1", "true", "TRUE", "yes", "YES")
    profile = Symbol(lowercase(get(ENV, "PHASE_PRECOMPILE_PROFILE", "scan")))
    Models.run_precompile_profile(profile)
end

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

if abspath(PROGRAM_FILE) == @__FILE__
    exit(PhaseCliSupport.main(Models, PROJECT_ROOT, collect(String.(ARGS))))
end
