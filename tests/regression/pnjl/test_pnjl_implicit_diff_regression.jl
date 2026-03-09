using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const IMPLICIT_BASELINE_PATH = joinpath(PROJECT_ROOT, "tests", "baselines", "pnjl", "baseline_pnjl_implicit_diff_fixedpoints_v1.csv")

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
        startswith(s, "#") && continue
        cols = split(s, ',')
        length(cols) == 18 || error("invalid baseline row: $line")
        push!(rows, (
            label=strip(cols[1]),
            T_fm=parse(Float64, strip(cols[2])),
            mu_fm=parse(Float64, strip(cols[3])),
            x=ntuple(i -> parse(Float64, strip(cols[3 + i])), 5),
            dx_dT=ntuple(i -> parse(Float64, strip(cols[8 + i])), 5),
            dx_dmu=ntuple(i -> parse(Float64, strip(cols[13 + i])), 5),
        ))
    end
    return rows
end

@testset "PNJL implicit diff regression" begin
    rows = _load_rows(IMPLICIT_BASELINE_PATH)
    PNJL = Models.pnjl_module()
    p_num = 24
    t_num = 6
    xi = 0.0
    rtol = 1e-6
    atol = 1e-10

    for row in rows
        @testset "$(row.label)" begin
            result = PNJL.solve_with_derivatives(row.T_fm, row.mu_fm; order=1, xi=xi, p_num=p_num, t_num=t_num)
            for i in 1:5
                @test isapprox(result.x[i], row.x[i]; rtol=rtol, atol=atol)
                @test isapprox(result.dx_dT[i], row.dx_dT[i]; rtol=rtol, atol=atol)
                @test isapprox(result.dx_dμ[i], row.dx_dmu[i]; rtol=rtol, atol=atol)
            end
        end
    end
end