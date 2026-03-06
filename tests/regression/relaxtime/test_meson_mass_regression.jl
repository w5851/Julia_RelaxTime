using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const MESON_BASELINE_PATH = joinpath(PROJECT_ROOT, "tests", "baselines", "relaxtime", "baseline_meson_mass_fixedpoints_v1.csv")

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const MesonMassWorkflow = Main.Models.meson_workflow_module()
using .MesonMassWorkflow

function _load_rows(path::String)
    isfile(path) || error("baseline CSV not found: $path")
    rows = NamedTuple[]
    for line in readlines(path)[2:end]
        s = strip(line)
        isempty(s) && continue
        cols = split(s, ',')
        push!(rows, (
            label=strip(cols[1]), T_fm=parse(Float64, strip(cols[2])), mu_fm=parse(Float64, strip(cols[3])), xi=parse(Float64, strip(cols[4])),
            masses=(parse(Float64, strip(cols[5])), parse(Float64, strip(cols[6])), parse(Float64, strip(cols[7]))),
            pi_converged=parse(Int, strip(cols[8])) == 1, pi_mass=parse(Float64, strip(cols[9])), pi_threshold=parse(Float64, strip(cols[10])), pi_gap=parse(Float64, strip(cols[11])),
            k_converged=parse(Int, strip(cols[12])) == 1, k_mass=parse(Float64, strip(cols[13])), k_threshold=parse(Float64, strip(cols[14])), k_gap=parse(Float64, strip(cols[15])),
        ))
    end
    return rows
end

@testset "Meson mass workflow regression" begin
    rows = _load_rows(MESON_BASELINE_PATH)
    for row in rows
        res = MesonMassWorkflow.solve_gap_and_meson_point(row.T_fm, row.mu_fm; xi=row.xi, mesons=(:pi, :K), p_num=8, t_num=4, solver_kwargs=(iterations=30,), mass_kwargs=(iterations=30,))
        rpi = res.meson_results[:pi]
        rk = res.meson_results[:K]
        @test isapprox(res.equilibrium.masses[1], row.masses[1]; rtol=1e-6, atol=1e-10)
        @test isapprox(res.equilibrium.masses[2], row.masses[2]; rtol=1e-6, atol=1e-10)
        @test isapprox(res.equilibrium.masses[3], row.masses[3]; rtol=1e-6, atol=1e-10)
        @test rpi.converged == row.pi_converged
        @test rk.converged == row.k_converged
        @test isapprox(rpi.mass, row.pi_mass; rtol=1e-6, atol=1e-10)
        @test isapprox(rpi.threshold, row.pi_threshold; rtol=1e-6, atol=1e-10)
        @test isapprox(rpi.gap, row.pi_gap; rtol=1e-6, atol=1e-10)
        @test isapprox(rk.mass, row.k_mass; rtol=1e-6, atol=1e-10)
        @test isapprox(rk.threshold, row.k_threshold; rtol=1e-6, atol=1e-10)
        @test isapprox(rk.gap, row.k_gap; rtol=1e-6, atol=1e-10)
    end
end