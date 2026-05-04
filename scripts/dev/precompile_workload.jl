#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(@__DIR__, "..", "..", "src", "models", "Models.jl"))
using .Models

# Representative compile workload for AD-heavy solver/derivative paths.
Models.run_precompile_profile(:full)

# Real scan pipeline path that dominates current cold-start cost.
scan_output_dir = mktempdir()
Models.run_scan_pipeline(
    :tmu;
    model_kind=:PNJL,
    T_values=[150.0],
    mu_values=[0.0, 100.0],
    xi_values=[0.0],
    output_path=joinpath(scan_output_dir, "tmu.csv"),
    overwrite=true,
)
Models.run_scan_pipeline(
    :trho;
    model_kind=:PNJL,
    T_values=[150.0],
    rho_values=[0.1, 0.2],
    xi_values=[0.0],
    output_path=joinpath(scan_output_dir, "trho.csv"),
    overwrite=true,
)

# Light phase pipeline path used in integration smoke/core.
phase_output_dir = mktempdir()
Models.run_phase_pipeline(
    :PNJL;
    mode=:research,
    T_grid=[150.0],
    rho_grid=[0.1, 0.2, 0.3],
    xi=0.0,
    output_dir=phase_output_dir,
    profile=:smoke,
    solver_backend=:models,
    reverse_rho=true,
    seed_policy=:hybrid_continuity,
    p_num=6,
    t_num=3,
    iterations=4,
)
