using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const CONSERVED_CHARGE_BASELINE_PATH = joinpath(PROJECT_ROOT, "tests", "baselines", "pnjl", "baseline_pnjl_conserved_charge_fixedpoints_v1.csv")

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
        length(cols) == 29 || error("invalid baseline row: $line")
        push!(rows, (
            label=strip(cols[1]),
            T_fm=parse(Float64, strip(cols[2])),
            muB_fm=parse(Float64, strip(cols[3])),
            muQ_fm=parse(Float64, strip(cols[4])),
            muS_fm=parse(Float64, strip(cols[5])),
            V=parse(Float64, strip(cols[6])),
            chi1_B=parse(Float64, strip(cols[7])),
            chi2_B=parse(Float64, strip(cols[8])),
            chi3_B=parse(Float64, strip(cols[9])),
            chi4_B=parse(Float64, strip(cols[10])),
            chi1_Q=parse(Float64, strip(cols[11])),
            chi2_Q=parse(Float64, strip(cols[12])),
            chi3_Q=parse(Float64, strip(cols[13])),
            chi4_Q=parse(Float64, strip(cols[14])),
            chi1_S=parse(Float64, strip(cols[15])),
            chi2_S=parse(Float64, strip(cols[16])),
            chi3_S=parse(Float64, strip(cols[17])),
            chi4_S=parse(Float64, strip(cols[18])),
            chi11_BQ=parse(Float64, strip(cols[19])),
            chi11_BS=parse(Float64, strip(cols[20])),
            chi11_QS=parse(Float64, strip(cols[21])),
            C200=parse(Float64, strip(cols[22])),
            C020=parse(Float64, strip(cols[23])),
            C002=parse(Float64, strip(cols[24])),
            C110=parse(Float64, strip(cols[25])),
            C101=parse(Float64, strip(cols[26])),
            C011=parse(Float64, strip(cols[27])),
            Ssigma=parse(Float64, strip(cols[28])),
            kappa_sigma2=parse(Float64, strip(cols[29])),
        ))
    end
    return rows
end

@testset "PNJL conserved-charge regression" begin
    rows = _load_rows(CONSERVED_CHARGE_BASELINE_PATH)
    PNJL = Models.pnjl_module()
    p_num = 16
    t_num = 6
    xi = 0.0
    rtol = 1e-6
    atol = 1e-10

    for row in rows
        @testset "$(row.label)" begin
            kwargs = (; xi=xi, p_num=p_num, t_num=t_num)
            @test isapprox(PNJL.chi1_B(row.T_fm, row.muB_fm; kwargs...), row.chi1_B; rtol=rtol, atol=atol)
            @test isapprox(PNJL.chi2_B(row.T_fm, row.muB_fm; kwargs...), row.chi2_B; rtol=rtol, atol=atol)
            @test isapprox(PNJL.chi3_B(row.T_fm, row.muB_fm; kwargs...), row.chi3_B; rtol=rtol, atol=atol)
            @test isapprox(PNJL.chi4_B(row.T_fm, row.muB_fm; kwargs...), row.chi4_B; rtol=rtol, atol=atol)
            @test isapprox(PNJL.chi1_Q(row.T_fm, row.muB_fm, row.muQ_fm, row.muS_fm; kwargs...), row.chi1_Q; rtol=rtol, atol=atol)
            @test isapprox(PNJL.chi2_Q(row.T_fm, row.muB_fm, row.muQ_fm, row.muS_fm; kwargs...), row.chi2_Q; rtol=rtol, atol=atol)
            @test isapprox(PNJL.chi3_Q(row.T_fm, row.muB_fm, row.muQ_fm, row.muS_fm; kwargs...), row.chi3_Q; rtol=rtol, atol=atol)
            @test isapprox(PNJL.chi4_Q(row.T_fm, row.muB_fm, row.muQ_fm, row.muS_fm; kwargs...), row.chi4_Q; rtol=rtol, atol=atol)
            @test isapprox(PNJL.chi1_S(row.T_fm, row.muB_fm, row.muQ_fm, row.muS_fm; kwargs...), row.chi1_S; rtol=rtol, atol=atol)
            @test isapprox(PNJL.chi2_S(row.T_fm, row.muB_fm, row.muQ_fm, row.muS_fm; kwargs...), row.chi2_S; rtol=rtol, atol=atol)
            @test isapprox(PNJL.chi3_S(row.T_fm, row.muB_fm, row.muQ_fm, row.muS_fm; kwargs...), row.chi3_S; rtol=rtol, atol=atol)
            @test isapprox(PNJL.chi4_S(row.T_fm, row.muB_fm, row.muQ_fm, row.muS_fm; kwargs...), row.chi4_S; rtol=rtol, atol=atol)
            @test isapprox(PNJL.chi11_BQ(row.T_fm, row.muB_fm, row.muQ_fm, row.muS_fm; kwargs...), row.chi11_BQ; rtol=rtol, atol=atol)
            @test isapprox(PNJL.chi11_BS(row.T_fm, row.muB_fm, row.muQ_fm, row.muS_fm; kwargs...), row.chi11_BS; rtol=rtol, atol=atol)
            @test isapprox(PNJL.chi11_QS(row.T_fm, row.muB_fm, row.muQ_fm, row.muS_fm; kwargs...), row.chi11_QS; rtol=rtol, atol=atol)
            @test isapprox(PNJL.cumulant_BQS(row.T_fm, row.muB_fm, row.muQ_fm, row.muS_fm, row.V; orders=(2, 0, 0), kwargs...), row.C200; rtol=rtol, atol=atol)
            @test isapprox(PNJL.cumulant_BQS(row.T_fm, row.muB_fm, row.muQ_fm, row.muS_fm, row.V; orders=(0, 2, 0), kwargs...), row.C020; rtol=rtol, atol=atol)
            @test isapprox(PNJL.cumulant_BQS(row.T_fm, row.muB_fm, row.muQ_fm, row.muS_fm, row.V; orders=(0, 0, 2), kwargs...), row.C002; rtol=rtol, atol=atol)
            @test isapprox(PNJL.cumulant_BQS(row.T_fm, row.muB_fm, row.muQ_fm, row.muS_fm, row.V; orders=(1, 1, 0), kwargs...), row.C110; rtol=rtol, atol=atol)
            @test isapprox(PNJL.cumulant_BQS(row.T_fm, row.muB_fm, row.muQ_fm, row.muS_fm, row.V; orders=(1, 0, 1), kwargs...), row.C101; rtol=rtol, atol=atol)
            @test isapprox(PNJL.cumulant_BQS(row.T_fm, row.muB_fm, row.muQ_fm, row.muS_fm, row.V; orders=(0, 1, 1), kwargs...), row.C011; rtol=rtol, atol=atol)
            @test isapprox(PNJL.baryon_Ssigma(row.T_fm, row.muB_fm; kwargs...), row.Ssigma; rtol=rtol, atol=atol)
            @test isapprox(PNJL.baryon_kappa_sigma2(row.T_fm, row.muB_fm; kwargs...), row.kappa_sigma2; rtol=rtol, atol=atol)
        end
    end
end