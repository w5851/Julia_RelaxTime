using Test

if !isdefined(Main, :validation_targets_path)
    Base.include(Main, joinpath(@__DIR__, "..", "..", "common", "data_paths.jl"))
end

const PNJL_CROSSOVER_LEGACY_SOURCE_CONSISTENCY_PATH = validation_targets_path(
    "pnjl",
    "reference",
    "pnjl_crossover_legacy_source_consistency_targets_v1.csv",
)

function _load_crossover_source_consistency_targets(path::String)
    isfile(path) || error("crossover source-consistency target CSV not found: $path")
    lines = readlines(path)
    isempty(lines) && error("crossover source-consistency target CSV is empty: $path")

    rows = NamedTuple[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        cols = split(s, ',')
        length(cols) == 7 || error("invalid crossover source-consistency row: $line")
        push!(rows, (
            record_id=strip(cols[1]),
            muB_MeV=parse(Float64, strip(cols[2])),
            mu_q_MeV=parse(Float64, strip(cols[3])),
            fortran_T_MeV=parse(Float64, strip(cols[4])),
            cpp_T_MeV=parse(Float64, strip(cols[5])),
            max_abs_delta_T_MeV=parse(Float64, strip(cols[6])),
            source=strip(cols[7]),
        ))
    end
    return rows
end

@testset "PNJL crossover legacy source consistency" begin
    targets = _load_crossover_source_consistency_targets(PNJL_CROSSOVER_LEGACY_SOURCE_CONSISTENCY_PATH)
    @test !isempty(targets)

    for row in targets
        delta_T = abs(row.cpp_T_MeV - row.fortran_T_MeV)
        @test delta_T <= row.max_abs_delta_T_MeV
    end
end