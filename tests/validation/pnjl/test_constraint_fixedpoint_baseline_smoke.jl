using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const BASELINE_PATH = joinpath(PROJECT_ROOT, "tests", "baselines", "pnjl", "baseline_pnjl_constraint_fixedpoints_v1.csv")
const HBARC_MEV_FM = 197.327

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const PNJL = Models.pnjl_module()

function _to_fm_inv(x_mev::Real)
    return Float64(x_mev) / HBARC_MEV_FM
end

function _load_baseline(path::String)
    isfile(path) || error("baseline CSV not found: $path")
    lines = readlines(path)
    isempty(lines) && error("baseline CSV is empty: $path")

    rows = NamedTuple[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue

        cols = split(s, ',')
        length(cols) == 16 || error("invalid baseline row: $line")

        parse_float = value -> parse(Float64, strip(value))
        parse_bool = value -> parse(Int, strip(value)) == 1

        push!(rows, (
            kind=strip(cols[1]),
            T_MeV=parse_float(cols[2]),
            mu_MeV=parse_float(cols[3]),
            rho_target=parse_float(cols[4]),
            pnjl_converged=parse_bool(cols[5]),
            models_converged=parse_bool(cols[6]),
            pnjl_pressure=parse_float(cols[7]),
            models_pressure=parse_float(cols[8]),
            pnjl_rho_norm=parse_float(cols[9]),
            models_rho_norm=parse_float(cols[10]),
            pnjl_entropy=parse_float(cols[11]),
            models_entropy=parse_float(cols[12]),
            pnjl_energy=parse_float(cols[13]),
            models_energy=parse_float(cols[14]),
            pnjl_residual_norm=parse_float(cols[15]),
            models_residual_norm=parse_float(cols[16]),
        ))
    end

    return rows
end

@testset "PNJL fixedpoint constraint baseline smoke" begin
    baseline = _load_baseline(BASELINE_PATH)
    model = Models.create_model(:PNJL)

    rtol = 1e-8
    atol = 1e-10

    for row in baseline
        T_fm = _to_fm_inv(row.T_MeV)

        if row.kind == "fixed_mu"
            mu_fm = _to_fm_inv(row.mu_MeV)

            pnjl_res = PNJL.solve(PNJL.FixedMu(), T_fm, mu_fm; p_num=8, t_num=4, iterations=80)
            models_res = Models.solve_constraint(model, Models.FixedMu(), T_fm;
                μ_fm=mu_fm,
                seed_guess=PNJL.HADRON_SEED_5,
                p_num=8,
                t_num=4,
            )

            @test pnjl_res.converged == row.pnjl_converged
            @test models_res.converged == row.models_converged

            @test isapprox(pnjl_res.pressure, row.pnjl_pressure; rtol=rtol, atol=atol)
            @test isapprox(models_res.pressure, row.models_pressure; rtol=rtol, atol=atol)
            @test isapprox(pnjl_res.rho_norm, row.pnjl_rho_norm; rtol=rtol, atol=atol)
            @test isapprox(models_res.rho_norm, row.models_rho_norm; rtol=rtol, atol=atol)
            @test isapprox(pnjl_res.entropy, row.pnjl_entropy; rtol=rtol, atol=atol)
            @test isapprox(models_res.entropy, row.models_entropy; rtol=rtol, atol=atol)
            @test isapprox(pnjl_res.energy, row.pnjl_energy; rtol=rtol, atol=atol)
            @test isapprox(models_res.energy, row.models_energy; rtol=rtol, atol=atol)
            @test isapprox(pnjl_res.residual_norm, row.pnjl_residual_norm; rtol=rtol, atol=atol)
            @test isapprox(models_res.residual_norm, row.models_residual_norm; rtol=rtol, atol=atol)
        elseif row.kind == "fixed_rho"
            pnjl_res = PNJL.solve(PNJL.FixedRho(row.rho_target), T_fm; p_num=8, t_num=4, iterations=80)
            models_res = Models.solve_constraint(model, Models.FixedRho(row.rho_target), T_fm;
                seed_guess=PNJL.HADRON_SEED_5,
                p_num=8,
                t_num=4,
            )

            @test pnjl_res.converged == row.pnjl_converged
            @test models_res.converged == row.models_converged

            @test isapprox(pnjl_res.pressure, row.pnjl_pressure; rtol=rtol, atol=atol)
            @test isapprox(models_res.pressure, row.models_pressure; rtol=rtol, atol=atol)
            @test isapprox(pnjl_res.rho_norm, row.pnjl_rho_norm; rtol=rtol, atol=atol)
            @test isapprox(models_res.rho_norm, row.models_rho_norm; rtol=rtol, atol=atol)
            @test isapprox(pnjl_res.entropy, row.pnjl_entropy; rtol=rtol, atol=atol)
            @test isapprox(models_res.entropy, row.models_entropy; rtol=rtol, atol=atol)
            @test isapprox(pnjl_res.energy, row.pnjl_energy; rtol=rtol, atol=atol)
            @test isapprox(models_res.energy, row.models_energy; rtol=rtol, atol=atol)
            @test isapprox(pnjl_res.residual_norm, row.pnjl_residual_norm; rtol=rtol, atol=atol)
            @test isapprox(models_res.residual_norm, row.models_residual_norm; rtol=rtol, atol=atol)
        else
            error("unknown baseline kind: $(row.kind)")
        end
    end
end
