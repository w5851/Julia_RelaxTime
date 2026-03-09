using Test

if !isdefined(Main, :_compute_relaxtime_literature_transport_point)
    include(joinpath(@__DIR__, "literature_validation_helpers.jl"))
end

const RELAXTIME_SIGMA_VALIDATION_DATA_PATH = validation_targets_path(
    "relaxtime",
    "literature",
    "sigma",
    "relaxtime_sigma_literature_targets_v1.csv",
)

function _load_sigma_targets(path::String)
    isfile(path) || error("validation target CSV not found: $path")
    lines = readlines(path)
    isempty(lines) && error("validation target CSV is empty: $path")

    rows = NamedTuple[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        cols = split(s, ',')
        length(cols) == 9 || error("invalid sigma validation target row: $line")
        push!(rows, (
            target_id=strip(cols[1]),
            series=strip(cols[2]),
            process=Symbol(strip(cols[3])),
            T_MeV=parse(Float64, strip(cols[4])),
            muB_MeV=parse(Float64, strip(cols[5])),
            sqrt_s_GeV=parse(Float64, strip(cols[6])),
            expected_sigma_mb=parse(Float64, strip(cols[7])),
            rtol=parse(Float64, strip(cols[8])),
            source=strip(cols[9]),
        ))
    end
    return rows
end

@testset "RelaxTime literature digitized sigma targets" begin
    targets = _load_sigma_targets(RELAXTIME_SIGMA_VALIDATION_DATA_PATH)
    @test !isempty(targets)

    for row in targets
        actual = _compute_relaxtime_literature_sigma_mb(
            row.T_MeV,
            row.muB_MeV,
            row.process,
            row.sqrt_s_GeV,
        )
        @test isfinite(actual)
        @test isapprox(actual, row.expected_sigma_mb; rtol=row.rtol, atol=0.0)
    end
end