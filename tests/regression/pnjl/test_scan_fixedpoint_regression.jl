using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SCAN_BASELINE_PATH = joinpath(PROJECT_ROOT, "tests", "baselines", "pnjl", "baseline_pnjl_scan_fixedpoints_v1.csv")

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const PNJL = Models.pnjl_module()
const pnjl_run_tmu_scan = getproperty(PNJL, :run_tmu_scan)
const pnjl_run_trho_scan = getproperty(PNJL, :run_trho_scan)

function _load_baseline(path::String)
    isfile(path) || error("baseline CSV not found: $path")
    lines = readlines(path)
    isempty(lines) && error("baseline CSV is empty: $path")

    rows = NamedTuple[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        cols = split(s, ',')
        length(cols) == 10 || error("invalid baseline row: $line")
        push!(rows, (
            kind=strip(cols[1]),
            T=parse(Float64, strip(cols[2])),
            mu=parse(Float64, strip(cols[3])),
            rho=parse(Float64, strip(cols[4])),
            xi=parse(Float64, strip(cols[5])),
            pressure=parse(Float64, strip(cols[6])),
            rho_out=parse(Float64, strip(cols[7])),
            mu_avg=parse(Float64, strip(cols[8])),
            entropy=parse(Float64, strip(cols[9])),
            energy=parse(Float64, strip(cols[10])),
        ))
    end
    return rows
end

function _parse_single_row(path::String)
    lines = readlines(path)
    length(lines) >= 2 || error("invalid csv output: $path")
    header = split(lines[1], ',')
    cols = split(lines[2], ',')
    return Dict(header[i] => cols[i] for i in eachindex(header))
end

function _f(row::AbstractDict, key::String)
    haskey(row, key) || error("missing column '$key'")
    return parse(Float64, row[key])
end

@testset "PNJL scan fixedpoint regression" begin
    baseline = _load_baseline(SCAN_BASELINE_PATH)
    rtol = 8e-2
    atol = 1e-6

    for row in baseline
        tmp_dir = mktempdir()
        out = joinpath(tmp_dir, "single.csv")

        if row.kind == "tmu"
            stats = pnjl_run_tmu_scan(
                T_values=[row.T],
                mu_values=[row.mu],
                xi_values=[row.xi],
                output_path=out,
                overwrite=true,
                resume=false,
                use_phase_aware=false,
                p_num=8,
                t_num=4,
                iterations=80,
            )
            @test stats.total == 1
            @test stats.success == 1

            result = _parse_single_row(out)
            @test isapprox(_f(result, "pressure_fm4"), row.pressure; rtol=rtol, atol=atol)
            @test isapprox(_f(result, "entropy_fm3"), row.entropy; rtol=rtol, atol=atol)
            @test isapprox(_f(result, "energy_fm4"), row.energy; rtol=rtol, atol=atol)
            @test isapprox(_f(result, "rho"), row.rho_out; rtol=rtol, atol=atol)
        elseif row.kind == "trho"
            stats = pnjl_run_trho_scan(
                T_values=[row.T],
                rho_values=[row.rho],
                xi_values=[row.xi],
                output_path=out,
                overwrite=true,
                resume=false,
                reverse_rho=false,
                p_num=8,
                t_num=4,
                iterations=80,
            )
            @test stats.total == 1
            @test stats.success == 1

            result = _parse_single_row(out)
            @test isapprox(_f(result, "pressure_fm4"), row.pressure; rtol=rtol, atol=atol)
            @test isapprox(_f(result, "entropy_fm3"), row.entropy; rtol=rtol, atol=atol)
            @test isapprox(_f(result, "energy_fm4"), row.energy; rtol=rtol, atol=atol)
            @test isapprox(_f(result, "mu_avg_MeV"), row.mu_avg; rtol=rtol, atol=atol)
        else
            error("unknown baseline kind: $(row.kind)")
        end
    end
end