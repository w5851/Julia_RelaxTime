if !isdefined(Main, :validation_targets_path)
    Base.include(Main, joinpath(@__DIR__, "..", "common", "data_paths.jl"))
end

if !isdefined(Main, :Models)
    Base.include(Main, joinpath(VALIDATION_PROJECT_ROOT, "src", "models", "Models.jl"))
end

function _load_phase_targets(path::String)
    isfile(path) || error("validation target CSV not found: $path")
    lines = readlines(path)
    isempty(lines) && error("validation target CSV is empty: $path")

    rows = NamedTuple[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        cols = split(s, ',')
        length(cols) == 9 || error("invalid validation target row: $line")
        push!(rows, (
            target_id=strip(cols[1]),
            observable=strip(cols[2]),
            expected_T_MeV=parse(Float64, strip(cols[3])),
            expected_mu_MeV=parse(Float64, strip(cols[4])),
            lower_T_MeV=parse(Float64, strip(cols[5])),
            upper_T_MeV=parse(Float64, strip(cols[6])),
            lower_mu_MeV=parse(Float64, strip(cols[7])),
            upper_mu_MeV=parse(Float64, strip(cols[8])),
            source=strip(cols[9]),
        ))
    end
    return rows
end

function _target_by_id(rows, target_id::String)
    for row in rows
        row.target_id == target_id && return row
    end
    error("validation target not found: $target_id")
end

function _load_crossover_reference_targets(path::String)
    isfile(path) || error("crossover validation target CSV not found: $path")
    lines = readlines(path)
    isempty(lines) && error("crossover validation target CSV is empty: $path")

    rows = NamedTuple[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        cols = split(s, ',')
        length(cols) == 10 || error("invalid crossover validation target row: $line")
        push!(rows, (
            target_id=strip(cols[1]),
            observable=strip(cols[2]),
            xi=parse(Float64, strip(cols[3])),
            mu_MeV=parse(Float64, strip(cols[4])),
            expected_T_MeV=parse(Float64, strip(cols[5])),
            lower_T_MeV=parse(Float64, strip(cols[6])),
            upper_T_MeV=parse(Float64, strip(cols[7])),
            method=strip(cols[8]),
            variable=strip(cols[9]),
            source=strip(cols[10]),
        ))
    end
    return rows
end