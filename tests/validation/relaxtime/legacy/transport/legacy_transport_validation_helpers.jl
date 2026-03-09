if !isdefined(Main, :_compute_relaxtime_literature_transport_point)
    Base.include(Main, joinpath(@__DIR__, "..", "..", "literature_validation_helpers.jl"))
end

if !isdefined(Main, :_compute_relaxtime_tau_point)
    Base.include(Main, joinpath(@__DIR__, "..", "..", "common", "tau_validation_helpers.jl"))
end

function _load_legacy_transport_scalar_targets(path::String)
    isfile(path) || error("legacy transport target CSV not found: $path")
    lines = readlines(path)
    isempty(lines) && error("legacy transport target CSV is empty: $path")

    rows = NamedTuple[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        cols = split(s, ',')
        length(cols) == 7 || error("invalid legacy transport scalar target row: $line")
        push!(rows, (
            record_id=strip(cols[1]),
            T_MeV=parse(Float64, strip(cols[2])),
            muB_MeV=parse(Float64, strip(cols[3])),
            field=strip(cols[4]),
            expected_value=parse(Float64, strip(cols[5])),
            rtol=parse(Float64, strip(cols[6])),
            source=strip(cols[7]),
        ))
    end
    return rows
end