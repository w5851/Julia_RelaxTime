using Dates
using SHA
using UUIDs: uuid4

const _PNJL_SCAN_JOBS_LOCK = ReentrantLock()
const _PNJL_SCAN_JOBS = Dict{String, Dict{String, Any}}()
const _PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const _SCAN_OUTPUT_ROOT = abspath(normpath(joinpath(_PROJECT_ROOT, "data", "outputs")))
const _PNJL_SCAN_JOB_SEQ = Ref(0)
const _PNJL_SCAN_MAX_RUNNING = 2
const _PNJL_SCAN_MAX_PENDING = 32
const _PNJL_SCAN_MAX_XI_POINTS = 64
const _PNJL_SCAN_SEMAPHORE = Base.Semaphore(_PNJL_SCAN_MAX_RUNNING)
const _PNJL_SCAN_DEFAULT_TIMEOUT_SECONDS = 0
const _PNJL_SCAN_FINISHED_KEEP_MAX = 64
const _PNJL_SCAN_FINISHED_TTL_SECONDS = 86400
const _SCAN_JOB_STATUS_QUEUED = "queued"
const _SCAN_JOB_STATUS_RUNNING = "running"
const _SCAN_JOB_STATUS_SUCCEEDED = "succeeded"
const _SCAN_JOB_STATUS_FAILED = "failed"
const _SCAN_JOB_STATUS_CANCELLED = "cancelled"
const _TERMINAL_SCAN_JOB_STATUSES = Set([
    _SCAN_JOB_STATUS_SUCCEEDED,
    _SCAN_JOB_STATUS_FAILED,
    _SCAN_JOB_STATUS_CANCELLED,
])
const _ALLOWED_SCAN_JOB_TRANSITIONS = Dict(
    _SCAN_JOB_STATUS_QUEUED => Set([_SCAN_JOB_STATUS_RUNNING, _SCAN_JOB_STATUS_CANCELLED]),
    _SCAN_JOB_STATUS_RUNNING => Set([
        _SCAN_JOB_STATUS_RUNNING,
        _SCAN_JOB_STATUS_SUCCEEDED,
        _SCAN_JOB_STATUS_FAILED,
        _SCAN_JOB_STATUS_CANCELLED,
    ]),
    _SCAN_JOB_STATUS_SUCCEEDED => Set{String}(),
    _SCAN_JOB_STATUS_FAILED => Set{String}(),
    _SCAN_JOB_STATUS_CANCELLED => Set{String}(),
)
const _PNJL_SCAN_PRUNE_METRICS = Dict{String, Int}(
    "total" => 0,
    "ttl" => 0,
    "keep_max" => 0,
)
const _PNJL_SCAN_IDEMPOTENCY_CACHE = Dict{String, Dict{String, Any}}()
const _PNJL_SCAN_RUNTIME_METRICS = Dict{String, Int}(
    "job_succeeded" => 0,
    "job_failed" => 0,
    "job_cancelled" => 0,
    "duration_le_10s" => 0,
    "duration_10s_60s" => 0,
    "duration_gt_60s" => 0,
)

@inline function _json_response(status::Int, payload::AbstractDict)
    headers = [
        "Content-Type" => "application/json; charset=utf-8",
        "Access-Control-Allow-Origin" => "*",
    ]
    return HTTP.Response(status, headers, JSON3.write(payload))
end

@inline function _error_response(status::Int, code::String, message::String; extras::AbstractDict=Dict{String, Any}())
    return _json_response(status, _error_payload(code, message; extras=extras))
end

@inline function _env_nonneg_int(name::String, default::Int)
    raw = get(ENV, name, nothing)
    raw === nothing && return default
    value = try
        parse(Int, strip(raw))
    catch
        @warn "Invalid env integer, fallback to default" name=name raw=raw default=default
        return default
    end
    if value < 0
        @warn "Negative env integer, fallback to default" name=name value=value default=default
        return default
    end
    return value
end

@inline function _default_scan_job_timeout_seconds()
    return _env_nonneg_int("PNJL_SCAN_JOB_TIMEOUT_SECONDS", _PNJL_SCAN_DEFAULT_TIMEOUT_SECONDS)
end

@inline function _scan_jobs_policy()
    return (
        keep_max=_env_nonneg_int("PNJL_SCAN_FINISHED_KEEP_MAX", _PNJL_SCAN_FINISHED_KEEP_MAX),
        ttl_seconds=_env_nonneg_int("PNJL_SCAN_FINISHED_TTL_SECONDS", _PNJL_SCAN_FINISHED_TTL_SECONDS),
    )
end

@inline _timestamp_now() = Dates.format(Dates.now(), dateformat"yyyy-mm-ddTHH:MM:SS")

@inline function _is_valid_scan_job_status(status::String)
    return haskey(_ALLOWED_SCAN_JOB_TRANSITIONS, status)
end

function _new_job_event(job::Dict{String, Any}, event_code::String; extras::AbstractDict=Dict{String, Any}())
    event = Dict{String, Any}(
        "job_id" => get(job, "job_id", nothing),
        "kind" => get(job, "kind", nothing),
        "state" => get(job, "status", nothing),
        "timestamp" => _timestamp_now(),
        "event_code" => event_code,
    )
    for (k, v) in pairs(extras)
        event[string(k)] = v
    end
    return event
end

function _append_job_event!(job::Dict{String, Any}, event_code::String; extras::AbstractDict=Dict{String, Any}())
    events = get!(job, "events") do
        Vector{Dict{String, Any}}()
    end
    event = _new_job_event(job, event_code; extras=extras)
    push!(events, event)
    @info "PNJL_SCAN_JOB_EVENT" event
    return event
end

function _is_valid_job_status_transition(from_status::String, to_status::String)
    _is_valid_scan_job_status(from_status) || return false
    _is_valid_scan_job_status(to_status) || return false
    return to_status in _ALLOWED_SCAN_JOB_TRANSITIONS[from_status]
end

function _set_job_status!(job::Dict{String, Any}, to_status::String)
    _is_valid_scan_job_status(to_status) || throw(ArgumentError("invalid target job status: $(to_status)"))

    from_status = String(get(job, "status", ""))
    _is_valid_scan_job_status(from_status) || throw(ArgumentError("invalid current job status: $(from_status)"))
    _is_valid_job_status_transition(from_status, to_status) || throw(ArgumentError("illegal job status transition: $(from_status) -> $(to_status)"))

    job["status"] = to_status
    if to_status == _SCAN_JOB_STATUS_RUNNING && get(job, "started_at", nothing) === nothing
        job["started_at"] = _timestamp_now()
    end
    if to_status in _TERMINAL_SCAN_JOB_STATUSES
        job["ended_at"] = _timestamp_now()
    end
    return nothing
end

function _job_elapsed_seconds(job::Dict{String, Any}; now_dt::DateTime=Dates.now())
    started_at = _parse_job_ended_at(get(job, "started_at", nothing))
    started_at === nothing && return nothing
    return max(0, Dates.value(now_dt - started_at) ÷ 1000)
end

function _update_scan_runtime_metrics_on_terminal!(job::Dict{String, Any})
    status = String(get(job, "status", ""))
    if status == _SCAN_JOB_STATUS_SUCCEEDED
        _PNJL_SCAN_RUNTIME_METRICS["job_succeeded"] = get(_PNJL_SCAN_RUNTIME_METRICS, "job_succeeded", 0) + 1
    elseif status == _SCAN_JOB_STATUS_FAILED
        _PNJL_SCAN_RUNTIME_METRICS["job_failed"] = get(_PNJL_SCAN_RUNTIME_METRICS, "job_failed", 0) + 1
    elseif status == _SCAN_JOB_STATUS_CANCELLED
        _PNJL_SCAN_RUNTIME_METRICS["job_cancelled"] = get(_PNJL_SCAN_RUNTIME_METRICS, "job_cancelled", 0) + 1
    end

    started_dt = _parse_job_ended_at(get(job, "started_at", nothing))
    ended_dt = _parse_job_ended_at(get(job, "ended_at", nothing))
    if started_dt !== nothing && ended_dt !== nothing && ended_dt >= started_dt
        duration_seconds = Dates.value(ended_dt - started_dt) ÷ 1000
        if duration_seconds <= 10
            _PNJL_SCAN_RUNTIME_METRICS["duration_le_10s"] = get(_PNJL_SCAN_RUNTIME_METRICS, "duration_le_10s", 0) + 1
        elseif duration_seconds <= 60
            _PNJL_SCAN_RUNTIME_METRICS["duration_10s_60s"] = get(_PNJL_SCAN_RUNTIME_METRICS, "duration_10s_60s", 0) + 1
        else
            _PNJL_SCAN_RUNTIME_METRICS["duration_gt_60s"] = get(_PNJL_SCAN_RUNTIME_METRICS, "duration_gt_60s", 0) + 1
        end
    end
    return nothing
end

function _runtime_metrics_snapshot(queue_snapshot)
    return Dict{String, Any}(
        "terminal" => Dict(
            "succeeded" => get(_PNJL_SCAN_RUNTIME_METRICS, "job_succeeded", 0),
            "failed" => get(_PNJL_SCAN_RUNTIME_METRICS, "job_failed", 0),
            "cancelled" => get(_PNJL_SCAN_RUNTIME_METRICS, "job_cancelled", 0),
        ),
        "duration_buckets" => Dict(
            "le_10s" => get(_PNJL_SCAN_RUNTIME_METRICS, "duration_le_10s", 0),
            "between_10s_60s" => get(_PNJL_SCAN_RUNTIME_METRICS, "duration_10s_60s", 0),
            "gt_60s" => get(_PNJL_SCAN_RUNTIME_METRICS, "duration_gt_60s", 0),
        ),
        "queue" => Dict(
            "queued" => queue_snapshot.queued,
            "running" => queue_snapshot.running,
        ),
    )
end

function _mark_job_timeout!(job::Dict{String, Any})
    _set_job_status!(job, _SCAN_JOB_STATUS_FAILED)
    job["error"] = "scan execution timeout"
    job["reason_code"] = "TIMEOUT"
    job["internal_error"] = nothing
    _update_scan_runtime_metrics_on_terminal!(job)
    return nothing
end

function _maybe_mark_job_timeout!(job_id::String; now_ts::Union{Nothing, String}=nothing)
    now_dt = now_ts === nothing ? Dates.now() : DateTime(now_ts, dateformat"yyyy-mm-ddTHH:MM:SS")
    return lock(_PNJL_SCAN_JOBS_LOCK) do
        haskey(_PNJL_SCAN_JOBS, job_id) || return false
        job = _PNJL_SCAN_JOBS[job_id]
        status = String(get(job, "status", ""))
        status == _SCAN_JOB_STATUS_RUNNING || return false

        policy = get(job, "policy", Dict{String, Any}())
        timeout_seconds = Int(get(policy, "timeout_seconds", _default_scan_job_timeout_seconds()))
        timeout_seconds > 0 || return false

        elapsed = _job_elapsed_seconds(job; now_dt=now_dt)
        elapsed === nothing && return false
        elapsed >= timeout_seconds || return false

        _mark_job_timeout!(job)
        return true
    end
end

@inline function _to_bool(x, default::Bool)
    x === nothing && return default
    if x isa Bool
        return x
    end
    if x isa Number
        return Int(x) != 0
    end
    s = lowercase(String(x))
    if s in ("1", "true", "yes", "y", "on")
        return true
    elseif s in ("0", "false", "no", "n", "off")
        return false
    end
    error("Invalid boolean value: $(x)")
end

function _request_header(req::HTTP.Request, name::String)
    target = lowercase(name)
    for (k, v) in req.headers
        if lowercase(String(k)) == target
            return String(v)
        end
    end
    return nothing
end

function _collect_fingerprint_lines!(acc::Vector{String}, path::String, value)
    if value === nothing || value isa Number || value isa Bool || value isa AbstractString
        push!(acc, string(path, "=", repr(value)))
        return
    end

    if value isa AbstractVector || value isa Tuple
        for (idx, item) in pairs(value)
            _collect_fingerprint_lines!(acc, string(path, "[", idx, "]"), item)
        end
        return
    end

    if value isa AbstractDict
        key_list = sort!(collect(keys(value)); by=k -> String(k))
        for key in key_list
            child = get(value, key, nothing)
            child_path = isempty(path) ? String(key) : string(path, ".", String(key))
            _collect_fingerprint_lines!(acc, child_path, child)
        end
        return
    end

    if value isa JSON3.Object
        dict = _to_symbol_dict(value)
        _collect_fingerprint_lines!(acc, path, dict)
        return
    end

    push!(acc, string(path, "=", repr(string(value))))
end

function _request_fingerprint(kind::String, params::Dict{Symbol, Any})
    acc = String[]
    _collect_fingerprint_lines!(acc, "kind", lowercase(kind))
    _collect_fingerprint_lines!(acc, "params", params)
    payload = join(acc, '\n')
    return bytes2hex(sha1(payload))
end

function _extract_idempotency_key(req::HTTP.Request, body_dict::Dict{Symbol, Any}, params::Dict{Symbol, Any})
    key = _request_header(req, "Idempotency-Key")
    if key === nothing || isempty(strip(key))
        body_key = get(body_dict, :idempotency_key, get(params, :idempotency_key, nothing))
        body_key === nothing && return nothing
        key = String(body_key)
    end
    trimmed = strip(key)
    isempty(trimmed) && return nothing
    return String(trimmed)
end

function _check_idempotency_conflict_or_replay(idempotency_key::AbstractString, fingerprint::AbstractString)
    lock(_PNJL_SCAN_JOBS_LOCK) do
        if !haskey(_PNJL_SCAN_IDEMPOTENCY_CACHE, idempotency_key)
            return (mode=:new, job_id=nothing)
        end

        cached = _PNJL_SCAN_IDEMPOTENCY_CACHE[idempotency_key]
        cached_fp = String(get(cached, "fingerprint", ""))
        cached_job_id = String(get(cached, "job_id", ""))
        if cached_fp == fingerprint
            haskey(_PNJL_SCAN_JOBS, cached_job_id) || return (mode=:new, job_id=nothing)
            return (mode=:replay, job_id=cached_job_id)
        end
        return (mode=:conflict, job_id=cached_job_id)
    end
end

function _record_idempotency_entry!(idempotency_key::AbstractString, fingerprint::AbstractString, job_id::AbstractString)
    lock(_PNJL_SCAN_JOBS_LOCK) do
        _PNJL_SCAN_IDEMPOTENCY_CACHE[idempotency_key] = Dict{String, Any}(
            "fingerprint" => fingerprint,
            "job_id" => job_id,
            "created_at" => _timestamp_now(),
        )
    end
    return nothing
end

function _to_float64_vec(x, name::String)
    x === nothing && return nothing
    vals = Float64[]
    for item in x
        push!(vals, _to_float64(item))
    end
    isempty(vals) && error("$(name) cannot be empty")
    return vals
end

function _to_string(x, default::String)
    x === nothing && return default
    return String(x)
end

function _to_nonneg_int(x, default::Int; name::String)
    v = _to_int(x, default)
    v >= 0 || error("$(name) must be >= 0")
    return v
end

@inline function _normalized_path_for_compare(path::AbstractString)
    return replace(lowercase(abspath(normpath(path))), '/' => '\\')
end

function _ensure_safe_output_path(path_raw)
    output_raw = _to_string(path_raw, "")
    isempty(strip(output_raw)) && throw(ArgumentError("output_path cannot be empty"))

    candidate = if isabspath(output_raw)
        abspath(normpath(output_raw))
    else
        abspath(normpath(joinpath(_PROJECT_ROOT, output_raw)))
    end

    base_cmp = _normalized_path_for_compare(_SCAN_OUTPUT_ROOT)
    candidate_cmp = _normalized_path_for_compare(candidate)
    safe_prefix = string(base_cmp, "\\")
    if !(candidate_cmp == base_cmp || startswith(candidate_cmp, safe_prefix))
        throw(ArgumentError("output_path must be under data/outputs"))
    end

    return candidate
end

function _normalize_grid_dict(obj)
    raw = obj isa Dict ? obj : _to_symbol_dict(obj)
    return raw isa Dict ? raw : _to_symbol_dict(raw)
end

function _expand_xi_grid(grid_raw)
    grid = _normalize_grid_dict(grid_raw)
    start = _to_float64(get(grid, :start, get(grid, :min, nothing)))
    stop = _to_float64(get(grid, :stop, get(grid, :max, nothing)))
    step = _to_float64(get(grid, :step, nothing))

    step > 0 || error("xi_grid.step must be > 0")
    stop >= start || error("xi_grid.stop must be >= xi_grid.start")

    vals = Float64[]
    x = start
    while x <= stop + 1e-12
        push!(vals, x)
        x += step
        length(vals) <= _PNJL_SCAN_MAX_XI_POINTS || error("xi points exceed limit: $(_PNJL_SCAN_MAX_XI_POINTS)")
    end
    isempty(vals) && error("xi_grid generated no points")
    return vals
end

function _resolve_xi_values(params::Dict{Symbol, Any})
    has_xi = haskey(params, :xi)
    has_xi_values = haskey(params, :xi_values)
    has_xi_grid = haskey(params, :xi_grid)

    strategy_count = (has_xi ? 1 : 0) + (has_xi_values ? 1 : 0) + (has_xi_grid ? 1 : 0)
    strategy_count <= 1 || error("Use only one xi strategy: xi | xi_values | xi_grid")

    vals = if has_xi_values
        _to_float64_vec(params[:xi_values], "xi_values")
    elseif has_xi_grid
        _expand_xi_grid(params[:xi_grid])
    elseif has_xi
        [_to_float64(params[:xi])]
    else
        [0.0]
    end

    length(vals) <= _PNJL_SCAN_MAX_XI_POINTS || error("xi points exceed limit: $(_PNJL_SCAN_MAX_XI_POINTS)")
    return unique(vals)
end

function _resolve_scan_mode(params::Dict{Symbol, Any})
    mode_raw = get(params, :mode, "scan")
    mode = lowercase(String(mode_raw))
    mode in ("scan", "point") || error("mode must be scan or point")
    return mode
end

function _resolve_tmu_axes(params::Dict{Symbol, Any})
    mode = _resolve_scan_mode(params)
    if mode == "point"
        T_mev = _to_float64(get(params, :T_mev, nothing))
        mu_mev = _to_float64(get(params, :mu_mev, nothing))
        return [T_mev], [mu_mev]
    end

    t_values = haskey(params, :T_values) ? _to_float64_vec(params[:T_values], "T_values") : nothing
    mu_values = haskey(params, :mu_values) ? _to_float64_vec(params[:mu_values], "mu_values") : nothing
    return t_values, mu_values
end

function _resolve_trho_axes(params::Dict{Symbol, Any})
    mode = _resolve_scan_mode(params)
    if mode == "point"
        T_mev = _to_float64(get(params, :T_mev, nothing))
        rho_value = _to_float64(get(params, :rho_value, nothing))
        return [T_mev], [rho_value]
    end

    t_values = haskey(params, :T_values) ? _to_float64_vec(params[:T_values], "T_values") : nothing
    rho_values = haskey(params, :rho_values) ? _to_float64_vec(params[:rho_values], "rho_values") : nothing
    return t_values, rho_values
end

function _extract_params_dict(req::HTTP.Request)
    body = isempty(req.body) ? Dict{Symbol, Any}() : JSON3.read(String(req.body))
    body_dict = body isa Dict ? body : _to_symbol_dict(body)
    params_obj = haskey(body_dict, :params) ? body_dict[:params] : Dict{Symbol, Any}()
    params_dict = params_obj isa Dict ? params_obj : _to_symbol_dict(params_obj)
    return body_dict, params_dict
end

function _new_job(kind::String, request_snapshot::Dict{String, Any}; total_points::Union{Int, Nothing}=nothing)
    job_id = string(uuid4())
    seq = lock(_PNJL_SCAN_JOBS_LOCK) do
        _PNJL_SCAN_JOB_SEQ[] += 1
        _PNJL_SCAN_JOB_SEQ[]
    end
    job = Dict{String, Any}(
        "job_id" => job_id,
        "seq" => seq,
        "kind" => kind,
        "status" => _SCAN_JOB_STATUS_QUEUED,
        "created_at" => _timestamp_now(),
        "started_at" => nothing,
        "ended_at" => nothing,
        "progress" => Dict{String, Any}(
            "total" => total_points,
            "completed" => 0,
            "percent" => total_points === nothing ? nothing : 0.0,
        ),
        "result" => nothing,
        "error" => nothing,
        "internal_error" => nothing,
        "reason_code" => nothing,
        "request" => request_snapshot,
        "events" => Vector{Dict{String, Any}}(),
    )
    lock(_PNJL_SCAN_JOBS_LOCK) do
        _PNJL_SCAN_JOBS[job_id] = job
        _append_job_event!(job, "created")
        policy = _scan_jobs_policy()
        _prune_finished_jobs_locked!(policy.keep_max; ttl_seconds=policy.ttl_seconds)
    end
    return job_id
end

@inline function _parse_job_ended_at(raw)
    raw === nothing && return nothing
    if raw isa DateTime
        return raw
    end
    if raw isa AbstractString
        s = strip(raw)
        isempty(s) && return nothing
        try
            return DateTime(s, dateformat"yyyy-mm-ddTHH:MM:SS")
        catch
            try
                return DateTime(s)
            catch
                return nothing
            end
        end
    end
    return nothing
end

function _prune_finished_jobs_locked!(
    keep_max::Int=_PNJL_SCAN_FINISHED_KEEP_MAX;
    now::DateTime=Dates.now(),
    ttl_seconds::Int=_PNJL_SCAN_FINISHED_TTL_SECONDS,
)
    keep_max >= 0 || throw(ArgumentError("keep_max must be >= 0"))
    ttl_seconds >= 0 || throw(ArgumentError("ttl_seconds must be >= 0"))

    finished = Tuple{String, Int}[]
    expired_ids = String[]
    for (job_id, job) in _PNJL_SCAN_JOBS
        st = get(job, "status", "")
        if st in _TERMINAL_SCAN_JOB_STATUSES
            ended_at = _parse_job_ended_at(get(job, "ended_at", nothing))
            if ended_at !== nothing && (now - ended_at) > Dates.Second(ttl_seconds)
                push!(expired_ids, job_id)
                continue
            end
            push!(finished, (job_id, Int(get(job, "seq", typemax(Int)))))
        end
    end

    removed = 0
    removed_ttl = 0
    for job_id in expired_ids
        if haskey(_PNJL_SCAN_JOBS, job_id)
            delete!(_PNJL_SCAN_JOBS, job_id)
            removed += 1
            removed_ttl += 1
        end
    end

    excess = length(finished) - keep_max
    if excess <= 0
        _PNJL_SCAN_PRUNE_METRICS["ttl"] = get(_PNJL_SCAN_PRUNE_METRICS, "ttl", 0) + removed_ttl
        _PNJL_SCAN_PRUNE_METRICS["total"] = get(_PNJL_SCAN_PRUNE_METRICS, "total", 0) + removed
        return removed
    end

    sort!(finished; by=last)
    removed_keep_max = 0
    for i in 1:excess
        job_id = finished[i][1]
        if haskey(_PNJL_SCAN_JOBS, job_id)
            delete!(_PNJL_SCAN_JOBS, job_id)
            removed += 1
            removed_keep_max += 1
        end
    end
    _PNJL_SCAN_PRUNE_METRICS["ttl"] = get(_PNJL_SCAN_PRUNE_METRICS, "ttl", 0) + removed_ttl
    _PNJL_SCAN_PRUNE_METRICS["keep_max"] = get(_PNJL_SCAN_PRUNE_METRICS, "keep_max", 0) + removed_keep_max
    _PNJL_SCAN_PRUNE_METRICS["total"] = get(_PNJL_SCAN_PRUNE_METRICS, "total", 0) + removed
    return removed
end

function _queue_snapshot()
    lock(_PNJL_SCAN_JOBS_LOCK) do
        queued = 0
        running = 0
        for job in values(_PNJL_SCAN_JOBS)
            st = get(job, "status", "")
            st == "queued" && (queued += 1)
            st == "running" && (running += 1)
        end
        return (queued=queued, running=running)
    end
end

function _queue_position(job_id::String)
    lock(_PNJL_SCAN_JOBS_LOCK) do
        haskey(_PNJL_SCAN_JOBS, job_id) || return nothing
        job = _PNJL_SCAN_JOBS[job_id]
        get(job, "status", "") == "queued" || return nothing
        seq = Int(get(job, "seq", typemax(Int)))
        pos = 1
        for other in values(_PNJL_SCAN_JOBS)
            if get(other, "status", "") == "queued" && Int(get(other, "seq", typemax(Int))) < seq
                pos += 1
            end
        end
        return pos
    end
end

function _with_job(f::Function, job_id::String)
    lock(_PNJL_SCAN_JOBS_LOCK) do
        haskey(_PNJL_SCAN_JOBS, job_id) || return nothing
        return f(_PNJL_SCAN_JOBS[job_id])
    end
end

function _update_job_progress!(job::Dict{String, Any}, completed::Int)
    progress = job["progress"]::Dict{String, Any}
    total = progress["total"]
    progress["completed"] = completed
    if total isa Int && total > 0
        progress["percent"] = min(100.0, 100.0 * completed / total)
    else
        progress["percent"] = nothing
    end
    _append_job_event!(job, "progress"; extras=Dict(
        "completed" => progress["completed"],
        "total" => progress["total"],
        "percent" => progress["percent"],
    ))
end

function _safe_output_path(kind::String, params::Dict{Symbol, Any}, job_id::String)
    if haskey(params, :output_path)
        return _ensure_safe_output_path(params[:output_path])
    end

    ts = Dates.format(Dates.now(), dateformat"yyyymmdd_HHMMSS")
    if kind == "tmu"
        return joinpath(_PROJECT_ROOT, "data", "outputs", "results", "pnjl", "scan", "tmu", "tmu_scan_job_$(ts)_$(job_id[1:8]).csv")
    elseif kind == "trho"
        return joinpath(_PROJECT_ROOT, "data", "outputs", "results", "pnjl", "scan", "trho", "trho_scan_job_$(ts)_$(job_id[1:8]).csv")
    else
        error("Unsupported scan kind: $(kind)")
    end
end

@inline function _build_job_diagnostics(job_id, kind, job_status; err=nothing)
    payload = Dict{String, Any}(
        "job_id" => job_id,
        "kind" => kind,
        "job_status" => job_status,
    )
    if err !== nothing
        payload["error_type"] = string(typeof(err))
    end
    if err === nothing && job_id isa AbstractString
        reason_code = _with_job(String(job_id)) do job
            get(job, "reason_code", nothing)
        end
        reason_code !== nothing && (payload["reason_code"] = reason_code)
    end
    return payload
end

function _estimate_total_points(kind::String, params::Dict{Symbol, Any})
    xi_values = _resolve_xi_values(params)
    n_xi = length(xi_values)
    if kind == "tmu"
        t_values, mu_values = _resolve_tmu_axes(params)
        if t_values === nothing || mu_values === nothing
            return nothing
        end
        return length(t_values) * length(mu_values) * n_xi
    elseif kind == "trho"
        t_values, rho_values = _resolve_trho_axes(params)
        if t_values === nothing || rho_values === nothing
            return nothing
        end
        return length(t_values) * length(rho_values) * n_xi
    end
    return nothing
end

function _run_tmu_scan_job!(job_id::String, params::Dict{Symbol, Any})
    output_path = _safe_output_path("tmu", params, job_id)
    kwargs = Dict{Symbol, Any}(
        :output_path => output_path,
        :overwrite => haskey(params, :overwrite) ? _to_bool(params[:overwrite], false) : false,
        :resume => haskey(params, :resume) ? _to_bool(params[:resume], true) : true,
        :use_phase_aware => haskey(params, :use_phase_aware) ? _to_bool(params[:use_phase_aware], true) : true,
        :bootstrap_multiseed => haskey(params, :bootstrap_multiseed) ? _to_bool(params[:bootstrap_multiseed], true) : true,
        :p_num => haskey(params, :p_num) ? _to_int(params[:p_num], 24) : 24,
        :t_num => haskey(params, :t_num) ? _to_int(params[:t_num], 8) : 8,
    )

    t_values, mu_values = _resolve_tmu_axes(params)
    t_values !== nothing && (kwargs[:T_values] = t_values)
    mu_values !== nothing && (kwargs[:mu_values] = mu_values)
    kwargs[:xi_values] = _resolve_xi_values(params)

    processed_ref = Ref(0)
    progress_cb = (_point, _result) -> begin
        processed_ref[] += 1
        _with_job(job_id) do job
            _update_job_progress!(job, processed_ref[])
        end
        nothing
    end

    stats = Models.run_tmu_scan(; kwargs..., progress_cb=progress_cb)
    return Dict{String, Any}(
        "output_path" => stats.output,
        "stats" => Dict(
            "total" => stats.total,
            "success" => stats.success,
            "failure" => stats.failure,
            "skipped" => stats.skipped,
        ),
    )
end

function _run_trho_scan_job!(job_id::String, params::Dict{Symbol, Any})
    output_path = _safe_output_path("trho", params, job_id)
    kwargs = Dict{Symbol, Any}(
        :output_path => output_path,
        :overwrite => haskey(params, :overwrite) ? _to_bool(params[:overwrite], false) : false,
        :resume => haskey(params, :resume) ? _to_bool(params[:resume], true) : true,
        :reverse_rho => haskey(params, :reverse_rho) ? _to_bool(params[:reverse_rho], true) : true,
        :p_num => haskey(params, :p_num) ? _to_int(params[:p_num], 24) : 24,
        :t_num => haskey(params, :t_num) ? _to_int(params[:t_num], 8) : 8,
    )

    t_values, rho_values = _resolve_trho_axes(params)
    t_values !== nothing && (kwargs[:T_values] = t_values)
    rho_values !== nothing && (kwargs[:rho_values] = rho_values)
    kwargs[:xi_values] = _resolve_xi_values(params)

    processed_ref = Ref(0)
    progress_cb = (_point, _result) -> begin
        processed_ref[] += 1
        _with_job(job_id) do job
            _update_job_progress!(job, processed_ref[])
        end
        nothing
    end

    stats = Models.run_trho_scan(; kwargs..., progress_cb=progress_cb)
    return Dict{String, Any}(
        "output_path" => stats.output,
        "stats" => Dict(
            "total" => stats.total,
            "success" => stats.success,
            "failure" => stats.failure,
            "skipped" => stats.skipped,
        ),
    )
end

function _execute_scan_job!(job_id::String, kind::String, params::Dict{Symbol, Any})
    if kind == "tmu"
        return _run_tmu_scan_job!(job_id, params)
    elseif kind == "trho"
        return _run_trho_scan_job!(job_id, params)
    end
    error("Unsupported scan kind: $(kind), expected tmu or trho")
end

function _launch_scan_job(job_id::String, kind::String, params::Dict{Symbol, Any})
    Threads.@spawn begin
        Base.acquire(_PNJL_SCAN_SEMAPHORE)
        try
            started = _with_job(job_id) do job
                if String(get(job, "status", "")) != _SCAN_JOB_STATUS_QUEUED
                    return false
                end
                _set_job_status!(job, _SCAN_JOB_STATUS_RUNNING)
                _append_job_event!(job, "started")
                return true
            end
            started === true || return

            attempt = 0
            max_retries = _with_job(job_id) do job
                policy = get(job, "policy", Dict{String, Any}())
                Int(get(policy, "max_retries", 0))
            end
            max_retries === nothing && (max_retries = 0)

            while true
                _maybe_mark_job_timeout!(job_id) && break
                try
                    result = _execute_scan_job!(job_id, kind, params)
                    _with_job(job_id) do job
                        String(get(job, "status", "")) == _SCAN_JOB_STATUS_RUNNING || return
                        _set_job_status!(job, _SCAN_JOB_STATUS_SUCCEEDED)
                        job["result"] = result
                        job["reason_code"] = nothing
                        stats = get(result, "stats", Dict{String, Any}())
                        total = get(stats, "total", nothing)
                        if total isa Int
                            _update_job_progress!(job, total)
                        end
                        _update_scan_runtime_metrics_on_terminal!(job)
                        _append_job_event!(job, "ended"; extras=Dict("outcome" => "succeeded"))
                    end
                    break
                catch e
                    if attempt < max_retries
                        attempt += 1
                        _with_job(job_id) do job
                            String(get(job, "status", "")) == _SCAN_JOB_STATUS_RUNNING || return
                            _set_job_status!(job, _SCAN_JOB_STATUS_RUNNING)
                            job["error"] = "scan execution retrying"
                            job["internal_error"] = sprint(showerror, e, catch_backtrace())
                            _append_job_event!(job, "retry"; extras=Dict("attempt" => attempt))
                        end
                    else
                        _with_job(job_id) do job
                            String(get(job, "status", "")) == _SCAN_JOB_STATUS_RUNNING || return
                            _set_job_status!(job, _SCAN_JOB_STATUS_FAILED)
                            job["error"] = "scan execution failed"
                            job["reason_code"] = "EXECUTION_ERROR"
                            job["internal_error"] = sprint(showerror, e, catch_backtrace())
                            _update_scan_runtime_metrics_on_terminal!(job)
                            _append_job_event!(job, "ended"; extras=Dict("outcome" => "failed"))
                        end
                        break
                    end
                end
            end
        finally
            Base.release(_PNJL_SCAN_SEMAPHORE)
        end
    end
end

function handle_pnjl_scan_job_create(req::HTTP.Request)
    req.method == "POST" || return HTTP.Response(405, ["Content-Type" => "text/plain"], "Method Not Allowed")

    try
        body_dict, params = _extract_params_dict(req)
        kind_raw = get(body_dict, :kind, get(body_dict, :scan_kind, ""))
        kind = lowercase(String(kind_raw))
        kind in ("tmu", "trho") || error("Missing/invalid kind; expected tmu or trho")

        idempotency_key = _extract_idempotency_key(req, body_dict, params)
        req_fingerprint = _request_fingerprint(kind, params)
        replay_mode = nothing
        replay_job_id = nothing
        if idempotency_key !== nothing
            idem = _check_idempotency_conflict_or_replay(idempotency_key, req_fingerprint)
            replay_mode = idem.mode
            replay_job_id = idem.job_id
            if replay_mode == :conflict
                return _error_response(
                    409,
                    "IDEMPOTENCY_KEY_CONFLICT",
                    "idempotency key reused with different payload";
                    extras=Dict(
                        "idempotency" => Dict(
                            "key" => idempotency_key,
                            "replayed" => false,
                            "conflict" => true,
                        ),
                        "job_id" => replay_job_id,
                    ),
                )
            elseif replay_mode == :replay && replay_job_id !== nothing
                pos = _queue_position(replay_job_id)
                return _json_response(202, Dict(
                    "status" => "accepted",
                    "job_id" => replay_job_id,
                    "kind" => kind,
                    "status_url" => "/api/modules/pnjl-scan/jobs/$(replay_job_id)",
                    "result_url" => "/api/modules/pnjl-scan/jobs/$(replay_job_id)/result",
                    "queue" => Dict(
                        "position" => pos,
                        "max_running" => _PNJL_SCAN_MAX_RUNNING,
                        "max_pending" => _PNJL_SCAN_MAX_PENDING,
                    ),
                    "idempotency" => Dict(
                        "key" => idempotency_key,
                        "replayed" => true,
                        "conflict" => false,
                    ),
                    "diagnostics" => _build_job_diagnostics(replay_job_id, kind, _SCAN_JOB_STATUS_QUEUED),
                ))
            end
        end

        snap = _queue_snapshot()
        if snap.queued >= _PNJL_SCAN_MAX_PENDING
            return _error_response(
                429,
                "QUEUE_FULL",
                "queue is full";
                extras=Dict(
                    "policy" => Dict(
                        "max_running" => _PNJL_SCAN_MAX_RUNNING,
                        "max_pending" => _PNJL_SCAN_MAX_PENDING,
                    ),
                ),
            )
        end

        max_retries = _to_nonneg_int(get(params, :max_retries, 0), 0; name="max_retries")
        max_retries <= 3 || error("max_retries must be <= 3")
        timeout_default = _default_scan_job_timeout_seconds()
        timeout_seconds = _to_nonneg_int(get(params, :timeout_seconds, timeout_default), timeout_default; name="timeout_seconds")

        total_points = _estimate_total_points(kind, params)
        request_snapshot = Dict{String, Any}(
            "kind" => kind,
            "params" => Dict(string(k) => v for (k, v) in params),
        )
        job_id = _new_job(kind, request_snapshot; total_points=total_points)
        _with_job(job_id) do job
            job["policy"] = Dict(
                "max_retries" => max_retries,
                "timeout_seconds" => timeout_seconds,
                "max_running" => _PNJL_SCAN_MAX_RUNNING,
                "max_pending" => _PNJL_SCAN_MAX_PENDING,
                "mode" => _resolve_scan_mode(params),
                "xi_strategy" => haskey(params, :xi_grid) ? "xi_grid" : (haskey(params, :xi_values) ? "xi_values" : (haskey(params, :xi) ? "xi" : "default")),
            )
        end
        if idempotency_key !== nothing
            _record_idempotency_entry!(idempotency_key, req_fingerprint, job_id)
        end
        _launch_scan_job(job_id, kind, params)

        pos = _queue_position(job_id)

        return _json_response(202, Dict(
            "status" => "accepted",
            "job_id" => job_id,
            "kind" => kind,
            "status_url" => "/api/modules/pnjl-scan/jobs/$(job_id)",
            "result_url" => "/api/modules/pnjl-scan/jobs/$(job_id)/result",
            "queue" => Dict(
                "position" => pos,
                "max_running" => _PNJL_SCAN_MAX_RUNNING,
                "max_pending" => _PNJL_SCAN_MAX_PENDING,
            ),
            "idempotency" => Dict(
                "key" => idempotency_key,
                "replayed" => false,
                "conflict" => false,
            ),
            "diagnostics" => _build_job_diagnostics(job_id, kind, _SCAN_JOB_STATUS_QUEUED),
        ))
    catch e
        return _error_response(
            400,
            "INVALID_REQUEST",
            sprint(showerror, e);
            extras=Dict(
                "diagnostics" => _build_job_diagnostics(nothing, nothing, "rejected"; err=e),
            ),
        )
    end
end

function handle_pnjl_scan_job_cancel(job_id::String)
    payload = _with_job(job_id) do job
        status = String(get(job, "status", ""))
        if status in _TERMINAL_SCAN_JOB_STATUSES
            return _error_payload(
                "JOB_NOT_CANCELLABLE",
                "job already in terminal state";
                extras=Dict(
                    "job_id" => job_id,
                    "job_status" => status,
                ),
            )
        end

        if !(status == _SCAN_JOB_STATUS_QUEUED || status == _SCAN_JOB_STATUS_RUNNING)
            return _error_payload(
                "JOB_NOT_CANCELLABLE",
                "job is not cancellable";
                extras=Dict(
                    "job_id" => job_id,
                    "job_status" => status,
                ),
            )
        end

        _set_job_status!(job, _SCAN_JOB_STATUS_CANCELLED)
        job["error"] = "scan job cancelled"
        job["reason_code"] = "CANCELLED"
        job["internal_error"] = nothing
        _update_scan_runtime_metrics_on_terminal!(job)
        _append_job_event!(job, "ended"; extras=Dict("outcome" => "cancelled"))

        return Dict{String, Any}(
            "status" => "ok",
            "job_id" => job_id,
            "kind" => job["kind"],
            "job_status" => job["status"],
            "diagnostics" => _build_job_diagnostics(job_id, job["kind"], job["status"]),
        )
    end

    payload === nothing && return _error_response(404, "JOB_NOT_FOUND", "job not found"; extras=Dict("job_id" => job_id))
    if payload["status"] == "error"
        return _json_response(409, payload)
    end
    return _json_response(200, payload)
end

function handle_pnjl_scan_job_status(job_id::String)
    pos = _queue_position(job_id)
    snap = _queue_snapshot()
    policy = _scan_jobs_policy()
    payload = _with_job(job_id) do job
        Dict{String, Any}(
            "status" => "ok",
            "job_id" => job["job_id"],
            "kind" => job["kind"],
            "job_status" => job["status"],
            "created_at" => job["created_at"],
            "started_at" => job["started_at"],
            "ended_at" => job["ended_at"],
            "progress" => job["progress"],
            "error" => (job["status"] == "failed" ? job["error"] : nothing),
            "queue" => Dict(
                "position" => pos,
                "queued" => snap.queued,
                "running" => snap.running,
                "max_running" => _PNJL_SCAN_MAX_RUNNING,
                "max_pending" => _PNJL_SCAN_MAX_PENDING,
            ),
            "policy" => get(job, "policy", Dict{String, Any}()),
            "events" => get(job, "events", Vector{Dict{String, Any}}()),
            "governance" => Dict(
                "finished_keep_max" => policy.keep_max,
                "finished_ttl_seconds" => policy.ttl_seconds,
                "pruned_total" => get(_PNJL_SCAN_PRUNE_METRICS, "total", 0),
                "pruned_ttl" => get(_PNJL_SCAN_PRUNE_METRICS, "ttl", 0),
                "pruned_keep_max" => get(_PNJL_SCAN_PRUNE_METRICS, "keep_max", 0),
            ),
            "metrics" => _runtime_metrics_snapshot(snap),
            "diagnostics" => _build_job_diagnostics(job["job_id"], job["kind"], job["status"]),
        )
    end

    payload === nothing && return _error_response(404, "JOB_NOT_FOUND", "job not found"; extras=Dict("job_id" => job_id))
    return _json_response(200, payload)
end

function handle_pnjl_scan_job_result(job_id::String)
    payload = _with_job(job_id) do job
        status = job["status"]
        result = job["result"]
        if status != "succeeded"
            return _error_payload(
                "JOB_NOT_SUCCEEDED",
                "job not succeeded";
                extras=Dict(
                    "job_id" => job_id,
                    "job_status" => status,
                    "diagnostics" => _build_job_diagnostics(job_id, job["kind"], status),
                ),
            )
        end
        output_path = get(result, "output_path", nothing)
        Dict{String, Any}(
            "status" => "ok",
            "job_id" => job_id,
            "job_status" => status,
            "result" => result,
            "metadata" => Dict(
                "output_exists" => output_path isa AbstractString ? isfile(output_path) : false,
                "output_mtime" => output_path isa AbstractString && isfile(output_path) ? string(Dates.unix2datetime(mtime(output_path))) : nothing,
            ),
            "diagnostics" => _build_job_diagnostics(job_id, job["kind"], status),
        )
    end

    payload === nothing && return _error_response(404, "JOB_NOT_FOUND", "job not found"; extras=Dict("job_id" => job_id))
    if payload["status"] == "error"
        return _json_response(409, payload)
    end
    return _json_response(200, payload)
end
