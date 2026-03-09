using Test

if !isdefined(Main, :_compute_relaxtime_literature_transport_point)
    include(joinpath(@__DIR__, "literature_validation_helpers.jl"))
end

const RELAXTIME_USBAR_LEGACY_SIGMA_DATA_PATH = validation_targets_path(
    "relaxtime",
    "legacy",
    "sigma",
    "relaxtime_sigma_usbar_mu0_legacy_consensus_targets_v1.csv",
)

function _load_usbar_legacy_sigma_targets(path::String)
    isfile(path) || error("legacy sigma validation target CSV not found: $path")
    lines = readlines(path)
    isempty(lines) && error("legacy sigma validation target CSV is empty: $path")

    rows = NamedTuple[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        cols = split(s, ',')
        length(cols) == 12 || error("invalid legacy sigma validation target row: $line")
        push!(rows, (
            target_id=strip(cols[1]),
            series=strip(cols[2]),
            process=Symbol(strip(cols[3])),
            T_MeV=parse(Float64, strip(cols[4])),
            muB_MeV=parse(Float64, strip(cols[5])),
            sqrt_s_GeV=parse(Float64, strip(cols[6])),
            fortran_sigma_mb=parse(Float64, strip(cols[7])),
            cpp_sigma_mb=parse(Float64, strip(cols[8])),
            consensus_sigma_mb=parse(Float64, strip(cols[9])),
            relative_delta=parse(Float64, strip(cols[10])),
            rtol=parse(Float64, strip(cols[11])),
            source=strip(cols[12]),
        ))
    end
    return rows
end

@testset "RelaxTime legacy consensus usbar sigma targets" begin
    targets = _load_usbar_legacy_sigma_targets(RELAXTIME_USBAR_LEGACY_SIGMA_DATA_PATH)
    @test !isempty(targets)

    for row in targets
        actual = _compute_relaxtime_literature_sigma_mb(
            row.T_MeV,
            row.muB_MeV,
            row.process,
            row.sqrt_s_GeV,
        )
        @test isfinite(actual)
        @test isapprox(actual, row.consensus_sigma_mb; rtol=row.rtol, atol=0.0)
    end
end