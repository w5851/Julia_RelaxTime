using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const BASELINE_PATH = joinpath(PROJECT_ROOT, "tests", "baselines", "relaxtime", "baseline_meson_density_regimes_v1.csv")

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

function _load_rows(path::String)
    isfile(path) || error("baseline CSV not found: $path")
    rows = NamedTuple[]
    for line in readlines(path)[2:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        cols = split(s, ','; keepempty=true)
        getf(i) = (i <= length(cols) ? strip(cols[i]) : "")
        parse_or_empty(::Type{Float64}, s::AbstractString) = isempty(s) ? nothing : parse(Float64, s)
        parse_or_empty(::Type{Int}, s::AbstractString) = isempty(s) ? nothing : parse(Int, s)
        parse_or_empty(::Type{Bool}, s::AbstractString) = isempty(s) ? nothing : lowercase(String(s)) == "true"
        push!(rows, (
            label=getf(1),
            T_fm=parse(Float64, getf(2)),
            mu_fm=parse(Float64, getf(3)),
            xi=parse(Float64, getf(4)),
            regime=Symbol(getf(5)),
            scheme_or_stage=getf(6),
            q_nodes=parse(Int, getf(7)),
            qmax=parse_or_empty(Float64, getf(8)),
            omega_min=parse_or_empty(Float64, getf(9)),
            omega_max=parse_or_empty(Float64, getf(10)),
            omega_nodes=parse_or_empty(Int, getf(11)),
            eta=parse_or_empty(Float64, getf(12)),
            gamma_zero_tol=parse_or_empty(Float64, getf(13)),
            solver_iterations=parse_or_empty(Int, getf(14)),
            pole_residual_norm_max=parse_or_empty(Float64, getf(15)),
            pole_require_converged=parse_or_empty(Bool, getf(16)),
            n_pi=parse(Float64, getf(17)),
            n_K=parse(Float64, getf(18)),
            kpi_ratio=parse(Float64, getf(19)),
        ))
    end
    return rows
end

function _run_row(row)
    common_kwargs = (
        xi=row.xi,
        mesons=(:pi, :K),
        p_num=8,
        t_num=4,
        solver_kwargs=(; iterations=20),
        mass_kwargs=(; iterations=20),
    )

    if row.regime === :stable
        res = Models.solve_gap_and_meson_density_point(
            row.T_fm,
            row.mu_fm;
            common_kwargs...,
            density_kwargs=(; num_q_nodes=row.q_nodes),
        )
        return res.meson_density
    elseif row.regime === :phase_shift
        scheme = row.scheme_or_stage == "gbu_reference" ? :gbu_reference : :current
        res = Models.solve_gap_and_phase_shift_meson_density_point(
            row.T_fm,
            row.mu_fm;
            common_kwargs...,
            density_kwargs=(;
                scheme=scheme,
                qmax=row.qmax,
                q_nodes=row.q_nodes,
                omega_min=row.omega_min,
                omega_max=row.omega_max,
                omega_nodes=row.omega_nodes,
                eta=row.eta,
            ),
        )
        return res.phase_shift_meson_density
    elseif row.regime === :strict_bw
        stage = row.scheme_or_stage == "stage2_qpole" ? :stage2_qpole : :stage1_reduced
        res = Models.solve_gap_and_strict_bw_meson_density_point(
            row.T_fm,
            row.mu_fm;
            common_kwargs...,
            density_kwargs=(;
                stage=stage,
                qmax=row.qmax,
                q_nodes=row.q_nodes,
                omega_max=row.omega_max,
                omega_nodes=row.omega_nodes,
                gamma_zero_tol=row.gamma_zero_tol,
                solver_iterations=row.solver_iterations,
                pole_residual_norm_max=row.pole_residual_norm_max,
                pole_require_converged=row.pole_require_converged,
            ),
        )
        return res.strict_bw_meson_density
    else
        error("unsupported regime: $(row.regime)")
    end
end

@testset "Meson density regime regression" begin
    rows = _load_rows(BASELINE_PATH)
    rtol = 5e-3
    atol = 1e-8

    for row in rows
        @testset "$(row.label)" begin
            res = _run_row(row)
            @test isfinite(res.n_pi)
            @test isfinite(res.n_K)
            @test isfinite(res.kpi_ratio)
            @test isapprox(res.n_pi, row.n_pi; rtol=rtol, atol=atol)
            @test isapprox(res.n_K, row.n_K; rtol=rtol, atol=atol)
            @test isapprox(res.kpi_ratio, row.kpi_ratio; rtol=rtol, atol=atol)
        end
    end
end
