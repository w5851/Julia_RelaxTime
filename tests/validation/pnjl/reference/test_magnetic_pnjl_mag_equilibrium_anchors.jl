using Test
using StaticArrays
using LinearAlgebra: norm

if !isdefined(Main, :Models)
    include(joinpath(@__DIR__, "..", "..", "..", "..", "src", "models", "Models.jl"))
end
if !isdefined(Main, :validation_targets_path)
    include(joinpath(@__DIR__, "..", "..", "common", "data_paths.jl"))
end

const MAGNETIC_EQUILIBRIUM_TARGET_PATH = validation_targets_path(
    "pnjl", "reference", "pnjl_magnetic_external_equilibrium_anchor_targets_v1.csv",
)

function _load_magnetic_equilibrium_targets(path::String)
    isfile(path) || error("magnetic equilibrium target CSV not found: $path")
    lines = readlines(path)
    length(lines) == 4 || error("expected three magnetic equilibrium anchor rows, got $(length(lines) - 1)")
    header = split(lines[1], ',')
    rows = Dict{String, String}[]
    for line in lines[2:end]
        values = split(line, ',')
        length(header) == length(values) || error("invalid magnetic equilibrium target row")
        push!(rows, Dict(header[i] => values[i] for i in eachindex(header)))
    end
    return rows
end

@inline _magnetic_anchor_f(row, key) = parse(Float64, row[key])
@inline _magnetic_anchor_i(row, key) = parse(Int, row[key])

function _external_anchor_params(ext_hc::Float64)
    Λ_inv_fm = 602.3 / ext_hc
    return Models.PNJLCore.PNJLParams(
        profile="pnjl_mag_external_equilibrium_anchor",
        physics_profile="pnjl_mag_external_equilibrium_anchor",
        hbarc_MeV_fm=ext_hc,
        alpha_em=1 / 137.035999084,
        N_color=3,
        N_flavor=3,
        rho0_fm3=0.16,
        Λ_inv_fm=Λ_inv_fm,
        m_ud0_inv_fm=5.5 / ext_hc,
        m_s0_inv_fm=140.7 / ext_hc,
        G_fm2=1.835 / Λ_inv_fm^2,
        K_fm5=12.36 / Λ_inv_fm^5,
        thermal_p_max_inv_fm=40.0,
        T0_inv_fm=210.0 / ext_hc,
        polyakov_scheme=:log,
        a0=3.51,
        a1=-2.47,
        a2=15.2,
        a3=0.0,
        b3=-1.75,
        b4=7.555,
    )
end

function _external_anchor_imc(ext_hc::Float64)
    return Models.MagneticThermodynamics.MagneticIMCParams(
        a=0.0108805,
        b=-1.0133e-4,
        c=0.02228,
        d=1.84558e-4,
        Λ_QCD_MeV=300.0 * Main.Constants_PNJL.ħc_MeV_fm / ext_hc,
    )
end

function _anchor_seeds(state::SVector{5, Float64})
    return [
        state,
        state + SVector{5, Float64}(0.20, 0.20, 0.20, 0.0, 0.0),
        SVector{5, Float64}(-0.03, -0.03, -0.04, 0.50, 0.50),
    ]
end

@testset "pnjl_mag external equilibrium anchors" begin
    rows = _load_magnetic_equilibrium_targets(MAGNETIC_EQUILIBRIUM_TARGET_PATH)
    @test length(rows) == 3
    for row in rows
        @test row["source_commit"] == "e1fc81d3c3c9d220c49972e54307b66a604cb9db"
        @test row["acceptance_status"] == "accepted_equilibrium_anchor"

        ext_hc = _magnetic_anchor_f(row, "external_hc_MeV_fm")
        params = _external_anchor_params(ext_hc)
        base = Models.PNJLModel(params)
        conf = Models.MagneticThermodynamics.MagneticConfig(
            eB_fm2=_magnetic_anchor_f(row, "eB_fm_minus2"),
            n_max=_magnetic_anchor_i(row, "n_max"),
            p_num=_magnetic_anchor_i(row, "p_num"),
            pz_max=_magnetic_anchor_f(row, "pz_max_fm_inv"),
            imc=_external_anchor_imc(ext_hc),
            route=:mfir,
            zeta_num=_magnetic_anchor_i(row, "zeta_num"),
            params=params,
        )
        model = Models.PNJLMagneticModel(base, conf)
        @test model.magnetic.route == :mfir

        T_fm = _magnetic_anchor_f(row, "T_MeV") / ext_hc
        μ_fm = _magnetic_anchor_f(row, "muB_MeV") / (3.0 * ext_hc)
        mu_vec = SVector{3, Float64}(μ_fm, μ_fm, μ_fm)
        target_state = SVector{5, Float64}(
            _magnetic_anchor_f(row, "phi_u"),
            _magnetic_anchor_f(row, "phi_d"),
            _magnetic_anchor_f(row, "phi_s"),
            _magnetic_anchor_f(row, "Phi"),
            _magnetic_anchor_f(row, "PhiBar"),
        )
        seeds = _anchor_seeds(target_state)
        result = Models.solve_magnetic_gap(
            model,
            T_fm,
            mu_vec;
            p_num=_magnetic_anchor_i(row, "p_num"),
            pz_max=_magnetic_anchor_f(row, "pz_max_fm_inv"),
            n_max=_magnetic_anchor_i(row, "n_max"),
            initial_guess=seeds[1],
            seed_candidates=seeds[2:end],
            include_default_seeds=false,
            method=:trust_region,
            fallback_method=:newton,
            iterations=120,
            residual_norm_max=_magnetic_anchor_f(row, "residual_norm_max"),
            root_merge_tol=1e-5,
        )

        @test result.converged
        @test result.attempt_count >= length(seeds)
        distances = [norm(candidate.x_state - target_state) for candidate in result.candidates]
        @test !isempty(distances)
        candidate = result.candidates[argmin(distances)]
        @test candidate.converged
        @test candidate.physical
        @test candidate.n_max == _magnetic_anchor_i(row, "n_max")
        @test minimum(distances) <= _magnetic_anchor_f(row, "state_atol")
        @test candidate.residual_norm <= _magnetic_anchor_f(row, "residual_norm_max")
        @test isapprox(
            candidate.omega,
            _magnetic_anchor_f(row, "omega_fm4");
            rtol=_magnetic_anchor_f(row, "omega_rtol"),
            atol=_magnetic_anchor_f(row, "omega_atol"),
        )
    end
end
