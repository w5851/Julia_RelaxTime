module XiSmoothnessSampling

using Random

export SampleRow
export sample_params

struct SampleRow
    sample_id::String
    source::String
    T_MeV::Float64
    muq_MeV::Float64
    muB_MeV::Float64
    anchor_type::String
    anchor_T_MeV::Float64
    anchor_muq_MeV::Float64
    delta_T::Float64
    delta_muq::Float64
    rng_seed::Int
end

struct AnchorPoint
    anchor_type::String
    T_MeV::Float64
    muq_MeV::Float64
end

@inline _is_comment(line::AbstractString) = isempty(strip(line)) || startswith(strip(line), "#")

function _first_non_comment_line(io)
    for line in eachline(io)
        _is_comment(line) && continue
        return line
    end
    return nothing
end

function _parse_header_map(header_line::AbstractString)
    cols = split(strip(header_line), ',')
    return Dict{String, Int}(String(strip(c)) => i for (i, c) in enumerate(cols))
end

function _parse_float(parts::Vector{SubString{String}}, idx::Int)
    idx <= length(parts) || return nothing
    v = tryparse(Float64, strip(parts[idx]))
    return v === nothing ? nothing : Float64(v)
end

function _parse_required_float(parts::Vector{SubString{String}}, idx::Int;
    path::String,
    line_no::Int,
    col::String,
)
    idx <= length(parts) || throw(ArgumentError("$(path):$(line_no): missing column $(col) value"))
    raw = strip(parts[idx])
    v = tryparse(Float64, raw)
    v === nothing && throw(ArgumentError("$(path):$(line_no): invalid Float64 in column $(col): $(repr(raw))"))
    return Float64(v)
end

@inline _xi_is_zero(x::Float64) = abs(x) <= 1e-10

function _read_boundary_anchors(path::String)
    isfile(path) || throw(ArgumentError("boundary_csv does not exist: $(path)"))
    anchors = AnchorPoint[]
    open(path, "r") do io
        header = nothing
        colmap = Dict{String, Int}()
        header_seen = false
        line_no = 0
        idx_xi = 0
        idx_T = 0
        idx_mu = 0

        for raw in eachline(io)
            line_no += 1
            if !header_seen
                _is_comment(raw) && continue
                header = raw
                colmap = _parse_header_map(header)
                for col in ("xi", "T_MeV", "mu_transition_MeV")
                    haskey(colmap, col) || throw(ArgumentError("$(path): missing required column in boundary_csv: $(col)"))
                end
                idx_xi = colmap["xi"]
                idx_T = colmap["T_MeV"]
                idx_mu = colmap["mu_transition_MeV"]
                header_seen = true
                continue
            end

            _is_comment(raw) && continue
            parts = split(raw, ',')
            xi = _parse_required_float(parts, idx_xi; path=path, line_no=line_no, col="xi")
            T = _parse_required_float(parts, idx_T; path=path, line_no=line_no, col="T_MeV")
            mu_transition = _parse_required_float(parts, idx_mu; path=path, line_no=line_no, col="mu_transition_MeV")
            _xi_is_zero(xi) || continue
            push!(anchors, AnchorPoint("boundary", T, mu_transition / 3.0))
        end

        header_seen || return anchors
    end
    return anchors
end

function _read_crossover_anchors(path::String)
    isfile(path) || throw(ArgumentError("crossover_csv does not exist: $(path)"))
    anchors = AnchorPoint[]
    open(path, "r") do io
        header = nothing
        colmap = Dict{String, Int}()
        header_seen = false
        line_no = 0
        idx_xi = 0
        idx_mu = 0
        idx_T = 0

        for raw in eachline(io)
            line_no += 1
            if !header_seen
                _is_comment(raw) && continue
                header = raw
                colmap = _parse_header_map(header)
                for col in ("xi", "mu_MeV")
                    haskey(colmap, col) || throw(ArgumentError("$(path): missing required column in crossover_csv: $(col)"))
                end
                haskey(colmap, "T_crossover_chiral_MeV") || throw(ArgumentError("$(path): missing required column in crossover_csv: T_crossover_chiral_MeV"))
                idx_xi = colmap["xi"]
                idx_mu = colmap["mu_MeV"]
                idx_T = colmap["T_crossover_chiral_MeV"]
                header_seen = true
                continue
            end

            _is_comment(raw) && continue
            parts = split(raw, ',')
            xi = _parse_required_float(parts, idx_xi; path=path, line_no=line_no, col="xi")
            mu = _parse_required_float(parts, idx_mu; path=path, line_no=line_no, col="mu_MeV")
            T = _parse_required_float(parts, idx_T; path=path, line_no=line_no, col="T_crossover_chiral_MeV")
            _xi_is_zero(xi) || continue
            push!(anchors, AnchorPoint("crossover", T, mu / 3.0))
        end

        header_seen || return anchors
    end
    return anchors
end

function _validate_range(name::String, range::Tuple{Float64, Float64})
    lo, hi = range
    isfinite(lo) || throw(ArgumentError("$(name) lower bound must be finite"))
    isfinite(hi) || throw(ArgumentError("$(name) upper bound must be finite"))
    lo <= hi || throw(ArgumentError("$(name) lower bound must not exceed upper bound"))
end

@inline function _clamp_range(x::Float64, range::Tuple{Float64, Float64})
    lo, hi = range
    return clamp(x, lo, hi)
end

@inline function _uniform_in(rng::AbstractRNG, range::Tuple{Float64, Float64})
    lo, hi = range
    return lo + rand(rng) * (hi - lo)
end

function sample_params(total::Int;
    seed::Int=42,
    random_count::Int=12,
    near_count::Int=12,
    T_range::Tuple{Float64, Float64}=(50.0, 270.0),
    muq_range::Tuple{Float64, Float64}=(0.0, 360.0),
    boundary_csv::String=joinpath("data", "reference", "pnjl", "boundary.csv"),
    crossover_csv::String=joinpath("data", "reference", "pnjl", "crossover.csv"),
)
    total > 0 || throw(ArgumentError("total must be positive"))
    random_count >= 0 || throw(ArgumentError("random_count must be non-negative"))
    near_count >= 0 || throw(ArgumentError("near_count must be non-negative"))
    random_count + near_count == total || throw(ArgumentError("random_count + near_count must equal total"))

    _validate_range("T_range", T_range)
    _validate_range("muq_range", muq_range)

    rng = MersenneTwister(seed)

    rows = NamedTuple{(:seq, :source, :T_MeV, :muq_MeV, :anchor_type, :anchor_T_MeV, :anchor_muq_MeV, :delta_T, :delta_muq),
                      Tuple{Int, String, Float64, Float64, String, Float64, Float64, Float64, Float64}}[]

    seq = 0
    for _ in 1:random_count
        seq += 1
        T = _uniform_in(rng, T_range)
        muq = _uniform_in(rng, muq_range)
        push!(rows, (
            seq = seq,
            source = "random_uniform",
            T_MeV = T,
            muq_MeV = muq,
            anchor_type = "none",
            anchor_T_MeV = NaN,
            anchor_muq_MeV = NaN,
            delta_T = NaN,
            delta_muq = NaN,
        ))
    end

    if near_count > 0
        anchors = vcat(_read_boundary_anchors(boundary_csv), _read_crossover_anchors(crossover_csv))
        sort!(anchors, by=a -> (a.anchor_type, a.T_MeV, a.muq_MeV))
        isempty(anchors) && throw(ArgumentError("no xi=0 anchors found in boundary/crossover csv"))

        for i in 1:near_count
            seq += 1
            anchor = anchors[mod1(i, length(anchors))]
            raw_delta_T = (2.0 * rand(rng) - 1.0) * 18.0
            raw_delta_muq = (2.0 * rand(rng) - 1.0) * 36.0
            T = _clamp_range(anchor.T_MeV + raw_delta_T, T_range)
            muq = _clamp_range(anchor.muq_MeV + raw_delta_muq, muq_range)
            push!(rows, (
                seq = seq,
                source = "near_phase_line",
                T_MeV = T,
                muq_MeV = muq,
                anchor_type = anchor.anchor_type,
                anchor_T_MeV = anchor.T_MeV,
                anchor_muq_MeV = anchor.muq_MeV,
                delta_T = T - anchor.T_MeV,
                delta_muq = muq - anchor.muq_MeV,
            ))
        end
    end

    sort!(rows, by=r -> (
        r.source == "random_uniform" ? 0 : 1,
        r.T_MeV,
        r.muq_MeV,
        r.anchor_type,
        r.anchor_T_MeV,
        r.anchor_muq_MeV,
        r.delta_T,
        r.delta_muq,
        r.seq,
    ))

    out = SampleRow[]
    for (i, r) in enumerate(rows)
        sid = "S" * lpad(string(i), 3, '0')
        push!(out, SampleRow(
            sid,
            r.source,
            r.T_MeV,
            r.muq_MeV,
            3.0 * r.muq_MeV,
            r.anchor_type,
            r.anchor_T_MeV,
            r.anchor_muq_MeV,
            r.delta_T,
            r.delta_muq,
            seed,
        ))
    end

    return out
end

end
