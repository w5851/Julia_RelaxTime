using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const _REFERENCE_CSV = joinpath(
    PROJECT_ROOT,
    "data",
    "outputs",
    "results",
    "relaxtime",
    "mott_phase",
    "reference_100_300_fine",
    "mott_phase_scan.csv",
)

function _load_reference_row(T_target::Float64, xi_target::Float64)
    found = nothing
    open(_REFERENCE_CSV, "r") do io
        for line in eachline(io)
            stripped = strip(line)
            isempty(stripped) && continue
            startswith(stripped, "#") && continue
            if startswith(stripped, "run_id,")
                continue
            end
            cols = split(stripped, ',')
            length(cols) >= 7 || continue
            T_val = parse(Float64, cols[2])
            xi_val = parse(Float64, cols[4])
            if isapprox(T_val, T_target; atol=1e-9) && isapprox(xi_val, xi_target; atol=1e-9)
                found = (
                    M_pi=parse(Float64, cols[5]),
                    M_K=parse(Float64, cols[6]),
                )
                break
            end
        end
    end
    found === nothing && error("reference row not found for T=$(T_target), xi=$(xi_target)")
    return found
end

@testset "Meson density workflow matches reference Mott continuation masses" begin
    continuation_state = nothing
    result_190 = nothing

    for T_MeV in (150.0, 160.0, 170.0, 180.0, 190.0)
        T_fm = T_MeV / Main.Constants_PNJL.ħc_MeV_fm
        result_190 = Models.solve_gap_and_meson_density_point(
            T_fm,
            0.0;
            xi=0.0,
            mesons=(:pi, :K),
            continuation_state=continuation_state,
            mixed_branch_align=:strict_sign_binding,
            p_num=8,
            t_num=4,
            solver_kwargs=(; iterations=20),
            mass_kwargs=(; iterations=20),
            density_kwargs=(; num_q_nodes=64),
        )
        continuation_state = result_190.continuation_state
    end

    ref = _load_reference_row(190.0, 0.0)
    ref_m_pi = ref.M_pi
    ref_m_K = ref.M_K

    @test result_190 !== nothing
    @test result_190.meson_density.m_pi ≈ ref_m_pi rtol=5e-3
    @test result_190.meson_density.m_K ≈ ref_m_K rtol=5e-3
    @test result_190.meson_density.m_K > 0.0
end
