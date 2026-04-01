using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const HBARC_MEV_FM = 197.327

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

function _to_fm_inv(x_mev::Real)
    return Float64(x_mev) / HBARC_MEV_FM
end

@testset "Wave-D cleanup stability" begin
    model = Main.Models.create_model(:PNJL)

    T_fm = _to_fm_inv(150.0)
    mu_fm = _to_fm_inv(0.0)

    unified = Main.Models.solve_constraint(
        model,
        Main.Models.FixedMu(),
        T_fm;
        μ_fm=mu_fm,
        seed_guess=Main.Models.HADRON_SEED_5,
        p_num=8,
        t_num=4,
    )

    legacy_err = try
        Main.Models.solve_fixedmu_constraint(
            model,
            T_fm,
            mu_fm;
            seed_guess=Main.Models.HADRON_SEED_5,
            p_num=8,
            t_num=4,
        )
        nothing
    catch exc
        exc
    end

    @test unified.converged
    @test isfinite(unified.pressure)
    @test isfinite(unified.rho_norm)
    @test legacy_err isa ArgumentError
    @test occursin("hard-deprecated", sprint(showerror, legacy_err))
    @test occursin("Models.solve_constraint", sprint(showerror, legacy_err))
end
