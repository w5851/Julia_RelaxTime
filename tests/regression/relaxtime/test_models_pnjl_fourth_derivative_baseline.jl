using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const HBARC_MEV_FM = 197.327

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const FIXED_POINTS = [
    (label="T160_mu20", T_fm=160.0 / HBARC_MEV_FM, muB_fm=20.0 / HBARC_MEV_FM),
    (label="T190_mu0", T_fm=190.0 / HBARC_MEV_FM, muB_fm=0.0),
]

const CHI4_B_BASELINE = Dict(
    "T160_mu20" => (expected=0.01269147425644905, rtol=2e-5, atol=1e-9),
    "T190_mu0" => (expected=0.0633065946585104, rtol=2e-5, atol=1e-9),
)

@testset "Models PNJL chi4_B fixed-point baseline" begin
    kwargs = (; p_num=16, t_num=6)
    chi4_by_label = Dict{String, Float64}()

    for point in FIXED_POINTS
        chi4 = Models.chi4_B(point.T_fm, point.muB_fm; kwargs...)
        @test isfinite(chi4)
        chi4_by_label[point.label] = chi4
    end

    for point in FIXED_POINTS
        chi4 = chi4_by_label[point.label]
        ref = CHI4_B_BASELINE[point.label]
        @test isapprox(chi4, ref.expected; rtol=ref.rtol, atol=ref.atol)
    end
end
