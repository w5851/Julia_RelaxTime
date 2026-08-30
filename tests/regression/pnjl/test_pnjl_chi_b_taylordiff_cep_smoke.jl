using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const _MODELS_PATH = normpath(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
const CEP_REFERENCE_PATH = joinpath(
    PROJECT_ROOT, "data", "reference", "pnjl", "issue130_phase_reference_v2",
    "accepted", "tables", "cep_boundary_accepted_phase_map_v1.csv",
)
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

function _load_xi0_cep(path::String)
    lines = readlines(path)
    isempty(lines) && error("CEP reference is empty: $path")
    header = split(strip(lines[1]), ',')
    idx_xi = findfirst(==("xi"), header)
    idx_T = findfirst(==("T_CEP_MeV"), header)
    idx_T_mid = findfirst(==("T_midpoint_MeV"), header)
    idx_muB = findfirst(==("muB_CEP_MeV"), header)
    idx_muq_proxy = findfirst(==("mu_CEP_proxy_MeV"), header)
    for line in lines[2:end]
        cols = split(strip(line), ',')
        parse(Float64, cols[idx_xi]) == 0.0 || continue
        return (
            T_CEP_MeV=parse(Float64, cols[idx_T === nothing ? idx_T_mid : idx_T]),
            muB_CEP_MeV=idx_muB === nothing ? 3.0 * parse(Float64, cols[idx_muq_proxy]) : parse(Float64, cols[idx_muB]),
        )
    end
    error("xi=0 CEP reference row not found: $path")
end

@testset "PNJL chi_B TaylorDiff CEP-neighborhood parity p8/t4" begin
    cep = _load_xi0_cep(CEP_REFERENCE_PATH)
    rows = _load_chi_b_fd_reference_rows(CHI_B_FD_REFERENCE_BASELINE_PATH)
    labels = ("cep_true_p8_t4", "cep_minus_0p5_p8_t4", "cep_plus_0p5_p8_t4")

    @test cep.muB_CEP_MeV > 800.0

    for label in labels
        subset = filter(row -> row.label == label, rows)
        @test length(subset) == 4
        by_order = Dict(row.order => row for row in subset)
        first_row = by_order[1]
        vals = TD.chi_B_taylordiff_all(
            first_row.T_fm,
            first_row.muB_fm;
            max_order=4,
            xi=first_row.xi,
            p_num=first_row.p_num,
            t_num=first_row.t_num,
            series_residual_tol=1e-6,
        )
        @testset "$(label)" begin
            for order in 1:4
                @test isapprox(vals[order], by_order[order].chi_B_forwarddiff; rtol=1e-8, atol=1e-9)
            end
        end
    end
end

@testset "PNJL chi_B TaylorDiff p24/t8 CEP returns orders 1..4" begin
    cep = _load_xi0_cep(CEP_REFERENCE_PATH)
    rows = filter(row -> row.label == "cep_true_p24_t8", _load_chi_b_fd_reference_rows(CHI_B_FD_REFERENCE_BASELINE_PATH))
    @test length(rows) == 3
    by_order = Dict(row.order => row for row in rows)
    first_row = by_order[1]

    vals = TD.chi_B_taylordiff_all(
        first_row.T_fm,
        first_row.muB_fm;
        max_order=4,
        xi=first_row.xi,
        p_num=first_row.p_num,
        t_num=first_row.t_num,
        series_residual_tol=1e-6,
    )

    @test cep.muB_CEP_MeV > 800.0
    @test length(vals) == 4
    @test all(isfinite, vals)

    for order in 1:3
        @test isapprox(vals[order], by_order[order].chi_B_forwarddiff; rtol=1e-8, atol=1e-8)
    end
end
