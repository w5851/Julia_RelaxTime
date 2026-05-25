using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const MIXEDJET_BASELINE_PATH = joinpath(PROJECT_ROOT, "tests", "baselines", "pnjl", "baseline_pnjl_chi_bqs_mixedjet_v1.csv")

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

function _load_mixedjet_rows(path::String)
    isfile(path) || error("baseline CSV not found: $path")
    lines = readlines(path)
    isempty(lines) && error("baseline CSV is empty: $path")

    rows = NamedTuple[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, '#') && continue
        cols = split(s, ',')
        length(cols) == 12 || error("invalid baseline row: $line")
        push!(rows, (
            label=strip(cols[1]),
            T_fm=parse(Float64, strip(cols[2])),
            muB_fm=parse(Float64, strip(cols[3])),
            muQ_fm=parse(Float64, strip(cols[4])),
            muS_fm=parse(Float64, strip(cols[5])),
            xi=parse(Float64, strip(cols[6])),
            p_num=parse(Int, strip(cols[7])),
            t_num=parse(Int, strip(cols[8])),
            orders=(
                parse(Int, strip(cols[9])),
                parse(Int, strip(cols[10])),
                parse(Int, strip(cols[11])),
            ),
            chi_BQS_mixedjet=parse(Float64, strip(cols[12])),
        ))
    end
    return rows
end

@testset "PNJL chi_BQS mixed Taylor jet fixed-point regression" begin
    PNJL = Models.pnjl_module()
    rows = _load_mixedjet_rows(MIXEDJET_BASELINE_PATH)
    @test length(rows) >= 3

    for row in rows
        @testset "$(row.label) orders=$(row.orders)" begin
            actual = PNJL.chi_BQS(
                row.T_fm,
                row.muB_fm,
                row.muQ_fm,
                row.muS_fm;
                orders=row.orders,
                xi=row.xi,
                p_num=row.p_num,
                t_num=row.t_num,
                derivative_backend=:mixedjet,
            )
            @test isapprox(actual, row.chi_BQS_mixedjet; rtol=1e-9, atol=1e-11)
        end
    end
end
