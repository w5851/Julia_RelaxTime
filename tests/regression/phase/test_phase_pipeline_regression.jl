using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const PHASE_BASELINE_PATH = joinpath(PROJECT_ROOT, "tests", "baselines", "phase", "baseline_phase_pipeline_v1.csv")

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

function _load_rows(path::String)
    isfile(path) || error("baseline CSV not found: $path")
    rows = NamedTuple[]
    for line in readlines(path)[2:end]
        s = strip(line)
        isempty(s) && continue
        cols = split(s, ',')
        push!(rows, (
            kind=strip(cols[1]), idx=parse(Int, strip(cols[2])), flag=parse(Int, strip(cols[3])) == 1,
            T_MeV=parse(Float64, strip(cols[4])), mu_MeV=parse(Float64, strip(cols[5])),
            rho_a=parse(Float64, strip(cols[6])), rho_b=parse(Float64, strip(cols[7])), aux=parse(Float64, strip(cols[8])),
        ))
    end
    return rows
end

@testset "Phase pipeline regression" begin
    rows = _load_rows(PHASE_BASELINE_PATH)
    tmp = mktempdir()
    result = Models.run_phase_pipeline(
        :PNJL;
        T_grid=[120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0],
        rho_grid=collect(0.1:0.1:3.0),
        xi=0.0,
        output_dir=tmp,
        profile=:regression,
        solver_backend=:legacy,
        reverse_rho=true,
        seed_policy=:hybrid_continuity,
        p_num=8,
        t_num=4,
        iterations=40,
        compute_crossover=true,
        crossover_method=:inflection,
        crossover_variable=:phi_u,
        crossover_n_mu=8,
        cep_strategy=:interpolate,
        promote_reference=false,
    )

    cep_row = only(filter(r -> r.kind == "cep", rows))
    @test result.cep.found == cep_row.flag
    @test isapprox(result.cep.T_cep_MeV, cep_row.T_MeV; rtol=1e-6, atol=1e-10)
    @test isapprox(result.cep.mu_cep_MeV, cep_row.mu_MeV; rtol=1e-6, atol=1e-10)

    boundary_rows = filter(r -> r.kind == "boundary", rows)
    @test length(result.first_order_boundary) == length(boundary_rows)
    for (row, actual) in zip(boundary_rows, result.first_order_boundary)
        @test actual.converged == row.flag
        @test isapprox(actual.T_MeV, row.T_MeV; rtol=1e-6, atol=1e-10)
        @test isapprox(actual.mu_transition_MeV, row.mu_MeV; rtol=1e-6, atol=1e-10)
        @test isapprox(actual.rho_hadron, row.rho_a; rtol=1e-6, atol=1e-10)
        @test isapprox(actual.rho_quark, row.rho_b; rtol=1e-6, atol=1e-10)
        @test isapprox(actual.area_residual, row.aux; rtol=1e-6, atol=1e-10)
    end

    crossover_rows = filter(r -> r.kind == "crossover", rows)
    @test length(result.crossover_line) == length(crossover_rows)
    for (row, actual) in zip(crossover_rows, result.crossover_line)
        @test actual.converged == row.flag
        @test isapprox(actual.T_crossover_MeV, row.T_MeV; rtol=1e-6, atol=1e-10, nans=true)
        @test isapprox(actual.mu_MeV, row.mu_MeV; rtol=1e-6, atol=1e-10)
        @test isapprox(actual.rho, row.rho_a; rtol=1e-6, atol=1e-10, nans=true)
        @test isapprox(actual.derivative, row.rho_b; rtol=1e-6, atol=1e-10, nans=true)
    end
end