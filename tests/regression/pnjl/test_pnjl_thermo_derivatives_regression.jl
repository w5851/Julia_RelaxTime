using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const THERMO_BASELINE_PATH = joinpath(PROJECT_ROOT, "tests", "baselines", "pnjl", "baseline_pnjl_thermo_derivatives_fixedpoints_v1.csv")

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
        length(cols) == 24 || error("invalid baseline row: $line")
        push!(rows, (
            label=strip(cols[1]),
            T_fm=parse(Float64, strip(cols[2])),
            mu_fm=parse(Float64, strip(cols[3])),
            pressure=parse(Float64, strip(cols[4])),
            energy=parse(Float64, strip(cols[5])),
            rho=parse(Float64, strip(cols[6])),
            entropy=parse(Float64, strip(cols[7])),
            dP_dT=parse(Float64, strip(cols[8])),
            dP_dmu=parse(Float64, strip(cols[9])),
            dEpsilon_dT=parse(Float64, strip(cols[10])),
            dEpsilon_dmu=parse(Float64, strip(cols[11])),
            dn_dT=parse(Float64, strip(cols[12])),
            dn_dmu=parse(Float64, strip(cols[13])),
            dP_depsilon_n=parse(Float64, strip(cols[14])),
            dP_dn_epsilon=parse(Float64, strip(cols[15])),
            dM_dT=(parse(Float64, strip(cols[16])), parse(Float64, strip(cols[17])), parse(Float64, strip(cols[18]))),
            dM_dmu=(parse(Float64, strip(cols[19])), parse(Float64, strip(cols[20])), parse(Float64, strip(cols[21]))),
            v_n_sq=parse(Float64, strip(cols[22])),
            dmuB_dT_sigma=parse(Float64, strip(cols[23])),
            n_B=parse(Float64, strip(cols[24])),
        ))
    end
    return rows
end

@testset "PNJL thermo derivatives regression" begin
    rows = _load_rows(THERMO_BASELINE_PATH)
    PNJL = Models.pnjl_module()
    p_num = 24
    t_num = 6
    xi = 0.0
    rtol = 1e-6
    atol = 1e-10

    for row in rows
        @testset "$(row.label)" begin
            td = PNJL.thermo_derivatives(row.T_fm, row.mu_fm; xi=xi, p_num=p_num, t_num=t_num)
            bv = PNJL.bulk_viscosity_coefficients(row.T_fm, row.mu_fm; xi=xi, p_num=p_num, t_num=t_num)

            @test isapprox(td.pressure, row.pressure; rtol=rtol, atol=atol)
            @test isapprox(td.energy, row.energy; rtol=rtol, atol=atol)
            @test isapprox(td.rho, row.rho; rtol=rtol, atol=atol)
            @test isapprox(td.entropy, row.entropy; rtol=rtol, atol=atol)
            @test isapprox(td.dP_dT, row.dP_dT; rtol=rtol, atol=atol)
            @test isapprox(td.dP_dmu, row.dP_dmu; rtol=rtol, atol=atol)
            @test isapprox(td.dEpsilon_dT, row.dEpsilon_dT; rtol=rtol, atol=atol)
            @test isapprox(td.dEpsilon_dmu, row.dEpsilon_dmu; rtol=rtol, atol=atol)
            @test isapprox(td.dn_dT, row.dn_dT; rtol=rtol, atol=atol)
            @test isapprox(td.dn_dmu, row.dn_dmu; rtol=rtol, atol=atol)
            @test isapprox(td.dP_depsilon_n, row.dP_depsilon_n; rtol=rtol, atol=atol)
            @test isapprox(td.dP_dn_epsilon, row.dP_dn_epsilon; rtol=rtol, atol=atol)
            @test isapprox(td.dM_dT[1], row.dM_dT[1]; rtol=rtol, atol=atol)
            @test isapprox(td.dM_dT[2], row.dM_dT[2]; rtol=rtol, atol=atol)
            @test isapprox(td.dM_dT[3], row.dM_dT[3]; rtol=rtol, atol=atol)
            @test isapprox(td.dM_dmu[1], row.dM_dmu[1]; rtol=rtol, atol=atol)
            @test isapprox(td.dM_dmu[2], row.dM_dmu[2]; rtol=rtol, atol=atol)
            @test isapprox(td.dM_dmu[3], row.dM_dmu[3]; rtol=rtol, atol=atol)
            @test isapprox(bv.v_n_sq, row.v_n_sq; rtol=rtol, atol=atol)
            @test isapprox(bv.dμB_dT_sigma, row.dmuB_dT_sigma; rtol=rtol, atol=atol)
            @test isapprox(bv.n_B, row.n_B; rtol=rtol, atol=atol)
        end
    end
end