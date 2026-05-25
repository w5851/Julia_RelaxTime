using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const _MODELS_PATH = normpath(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
const CHI_B_FD_REFERENCE_BASELINE_PATH = joinpath(PROJECT_ROOT, "tests", "baselines", "pnjl", "baseline_pnjl_chi_b_taylordiff_fd_reference_v1.csv")

if !isdefined(Main, :Models)
    Base.include(Main, _MODELS_PATH)
end

using .Models

const TD = Models.PNJLChiBTaylorDiff

function _load_chi_b_fd_reference_rows(path::String)
    isfile(path) || error("baseline CSV not found: $path")
    lines = readlines(path)
    isempty(lines) && error("baseline CSV is empty: $path")

    rows = NamedTuple[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, '#') && continue
        cols = split(s, ',')
        length(cols) == 8 || error("invalid baseline row: $line")
        push!(rows, (
            label=strip(cols[1]),
            T_fm=parse(Float64, strip(cols[2])),
            muB_fm=parse(Float64, strip(cols[3])),
            xi=parse(Float64, strip(cols[4])),
            p_num=parse(Int, strip(cols[5])),
            t_num=parse(Int, strip(cols[6])),
            order=parse(Int, strip(cols[7])),
            chi_B_forwarddiff=parse(Float64, strip(cols[8])),
        ))
    end
    return rows
end

@testset "PNJL chi_B TaylorDiff fixed-point FD-baseline parity" begin
    rows = filter(row -> row.label == "fixed_p8_t4", _load_chi_b_fd_reference_rows(CHI_B_FD_REFERENCE_BASELINE_PATH))
    @test length(rows) == 4

    by_order = Dict(row.order => row for row in rows)
    first_row = by_order[1]
    vals = TD.chi_B_taylordiff_all(
        first_row.T_fm,
        first_row.muB_fm;
        max_order=4,
        xi=first_row.xi,
        p_num=first_row.p_num,
        t_num=first_row.t_num,
    )

    for order in 1:4
        expected = by_order[order].chi_B_forwarddiff
        @test isapprox(vals[order], expected; rtol=1e-9, atol=1e-11)
    end
end
