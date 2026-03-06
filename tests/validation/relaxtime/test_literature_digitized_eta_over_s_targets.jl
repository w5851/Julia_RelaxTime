using Test

if !isdefined(Main, :_compute_relaxtime_literature_transport_point)
    include(joinpath(@__DIR__, "literature_validation_helpers.jl"))
end

const RELAXTIME_ETA_VALIDATION_DATA_PATH = joinpath(
    RELAXTIME_LITERATURE_VALIDATION_PROJECT_ROOT,
    "tests",
    "validation",
    "data",
    "relaxtime_eta_over_s_literature_targets_v1.csv",
)

function _load_eta_targets(path::String)
    isfile(path) || error("validation target CSV not found: $path")
    lines = readlines(path)
    isempty(lines) && error("validation target CSV is empty: $path")

    rows = NamedTuple[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        cols = split(s, ',')
        length(cols) == 7 || error("invalid eta validation target row: $line")
        push!(rows, (
            target_id=strip(cols[1]),
            series=strip(cols[2]),
            muB_MeV=parse(Float64, strip(cols[3])),
            T_MeV=parse(Float64, strip(cols[4])),
            expected_eta_over_s=parse(Float64, strip(cols[5])),
            rtol=parse(Float64, strip(cols[6])),
            source=strip(cols[7]),
        ))
    end
    return rows
end

@testset "RelaxTime literature digitized eta_over_s targets" begin
    targets = _load_eta_targets(RELAXTIME_ETA_VALIDATION_DATA_PATH)
    @test !isempty(targets)

    for row in targets
        actual = _compute_relaxtime_literature_transport_point(row.T_MeV, row.muB_MeV).eta_over_s
        @test isfinite(actual)
        @test isapprox(actual, row.expected_eta_over_s; rtol=row.rtol, atol=0.0)
    end
end