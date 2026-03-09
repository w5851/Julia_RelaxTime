using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const FLAVOR_IMPLICIT_BASELINE_PATH = joinpath(PROJECT_ROOT, "tests", "baselines", "pnjl", "baseline_pnjl_flavor_mu_implicit_diff_fixedpoints_v1.csv")

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

function _load_rows(path::String)
    isfile(path) || error("baseline CSV not found: $path")
    lines = readlines(path)
    isempty(lines) && error("baseline CSV is empty: $path")

    rows = NamedTuple[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, '#') && continue
        cols = split(s, ',')
        length(cols) == 30 || error("invalid baseline row: $line")
        push!(rows, (
            label=strip(cols[1]),
            T_fm=parse(Float64, strip(cols[2])),
            mu_vec=SVector(
                parse(Float64, strip(cols[3])),
                parse(Float64, strip(cols[4])),
                parse(Float64, strip(cols[5])),
            ),
            x=ntuple(i -> parse(Float64, strip(cols[5 + i])), 5),
            dx_dT=ntuple(i -> parse(Float64, strip(cols[10 + i])), 5),
            dx_dmu_u=ntuple(i -> parse(Float64, strip(cols[15 + i])), 5),
            dx_dmu_d=ntuple(i -> parse(Float64, strip(cols[20 + i])), 5),
            dx_dmu_s=ntuple(i -> parse(Float64, strip(cols[25 + i])), 5),
        ))
    end
    return rows
end

@testset "PNJL flavor-mu implicit diff regression" begin
    rows = _load_rows(FLAVOR_IMPLICIT_BASELINE_PATH)
    rtol = 1e-6
    atol = 1e-10
    p_num = 24
    t_num = 6
    xi = 0.0

    for row in rows
        @testset "$(row.label)" begin
            result = Models.solve_pnjl_with_flavor_mu_derivatives(row.T_fm, row.mu_vec; order=1, xi=xi, p_num=p_num, t_num=t_num)
            for i in 1:5
                @test isapprox(result.x[i], row.x[i]; rtol=rtol, atol=atol)
                @test isapprox(result.dx_dT[i], row.dx_dT[i]; rtol=rtol, atol=atol)
                @test isapprox(result.dx_dmu_vec[i, 1], row.dx_dmu_u[i]; rtol=rtol, atol=atol)
                @test isapprox(result.dx_dmu_vec[i, 2], row.dx_dmu_d[i]; rtol=rtol, atol=atol)
                @test isapprox(result.dx_dmu_vec[i, 3], row.dx_dmu_s[i]; rtol=rtol, atol=atol)
            end
        end
    end
end