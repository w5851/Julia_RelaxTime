using Dates
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
const _PNJL_SCAN_FINISHED_KEEP_MAX = 64

@inline function _json_response(status::Int, payload::AbstractDict)
    headers = [
        "Content-Type" => "application/json; charset=utf-8",
        "Access-Control-Allow-Origin" => "*",
    ]
    return HTTP.Response(status, headers, JSON3.write(payload))
end

@inline _timestamp_now() = Dates.format(Dates.now(), dateformat"yyyy-mm-ddTHH:MM:SS")

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
        "status" => "queued",
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
        "request" => request_snapshot,
    )
    lock(_PNJL_SCAN_JOBS_LOCK) do
        _PNJL_SCAN_JOBS[job_id] = job
        _prune_finished_jobs_locked!(_PNJL_SCAN_FINISHED_KEEP_MAX)
    end
    return job_id
end

function _prune_finished_jobs_locked!(keep_max::Int=_PNJL_SCAN_FINISHED_KEEP_MAX)
    keep_max >= 0 || throw(ArgumentError("keep_max must be >= 0"))

    finished = Tuple{String, Int}[]
    for (job_id, job) in _PNJL_SCAN_JOBS
        st = get(job, "status", "")
        if st == "succeeded" || st == "failed"
            push!(finished, (job_id, Int(get(job, "seq", typemax(Int)))))
        end
    end

    excess = length(finished) - keep_max
    excess <= 0 && return 0

    sort!(finished; by=last)
    removed = 0
    for i in 1:excess
        job_id = finished[i][1]
        if haskey(_PNJL_SCAN_JOBS, job_id)
            delete!(_PNJL_SCAN_JOBS, job_id)
            removed += 1
        end
    end
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
    return payload
end

function _estimate_total_points(kind::String, params::Dict{Symbol, Any})
    xi_values = _resolve_xi_values(params)
    n_xi = length(xi_values)
    if kind == "tmu"
        t_values = haskey(params, :T_values) ? _to_float64_vec(params[:T_values], "T_values") : nothing
        mu_values = haskey(params, :mu_values) ? _to_float64_vec(params[:mu_values], "mu_values") : nothing
        if t_values === nothing || mu_values === nothing
            return nothing
        end
        return length(t_values) * length(mu_values) * n_xi
    elseif kind == "trho"
        t_values = haskey(params, :T_values) ? _to_float64_vec(params[:T_values], "T_values") : nothing
        rho_values = haskey(params, :rho_values) ? _to_float64_vec(params[:rho_values], "rho_values") : nothing
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

    haskey(params, :T_values) && (kwargs[:T_values] = _to_float64_vec(params[:T_values], "T_values"))
    haskey(params, :mu_values) && (kwargs[:mu_values] = _to_float64_vec(params[:mu_values], "mu_values"))
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

    haskey(params, :T_values) && (kwargs[:T_values] = _to_float64_vec(params[:T_values], "T_values"))
    haskey(params, :rho_values) && (kwargs[:rho_values] = _to_float64_vec(params[:rho_values], "rho_values"))
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
            _with_job(job_id) do job
                job["status"] = "running"
                job["started_at"] = _timestamp_now()
            end

            attempt = 0
            max_retries = _with_job(job_id) do job
                policy = get(job, "policy", Dict{String, Any}())
                Int(get(policy, "max_retries", 0))
            end
            max_retries === nothing && (max_retries = 0)

            while true
                try
                    result = _execute_scan_job!(job_id, kind, params)
                    _with_job(job_id) do job
                        job["status"] = "succeeded"
                        job["result"] = result
                        job["ended_at"] = _timestamp_now()
                        stats = get(result, "stats", Dict{String, Any}())
                        total = get(stats, "total", nothing)
                        if total isa Int
                            _update_job_progress!(job, total)
                        end
                    end
                    break
                catch e
                    if attempt < max_retries
                        attempt += 1
                        _with_job(job_id) do job
                            job["status"] = "running"
                            job["error"] = sprint(showerror, e)
                        end
                    else
                        _with_job(job_id) do job
                            job["status"] = "failed"
                            job["error"] = sprint(showerror, e, catch_backtrace())
                            job["ended_at"] = _timestamp_now()
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

        snap = _queue_snapshot()
        if snap.queued >= _PNJL_SCAN_MAX_PENDING
            return _json_response(429, Dict(
                "status" => "error",
                "error" => "queue is full",
                "policy" => Dict(
                    "max_running" => _PNJL_SCAN_MAX_RUNNING,
                    "max_pending" => _PNJL_SCAN_MAX_PENDING,
                ),
            ))
        end

        max_retries = _to_nonneg_int(get(params, :max_retries, 0), 0; name="max_retries")
        max_retries <= 3 || error("max_retries must be <= 3")

        total_points = _estimate_total_points(kind, params)
        request_snapshot = Dict{String, Any}(
            "kind" => kind,
            "params" => Dict(string(k) => v for (k, v) in params),
        )
        job_id = _new_job(kind, request_snapshot; total_points=total_points)
        _with_job(job_id) do job
            job["policy"] = Dict(
                "max_retries" => max_retries,
                "max_running" => _PNJL_SCAN_MAX_RUNNING,
                "max_pending" => _PNJL_SCAN_MAX_PENDING,
                "xi_strategy" => haskey(params, :xi_grid) ? "xi_grid" : (haskey(params, :xi_values) ? "xi_values" : (haskey(params, :xi) ? "xi" : "default")),
            )
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
            "diagnostics" => _build_job_diagnostics(job_id, kind, "queued"),
        ))
    catch e
        return _json_response(400, Dict(
            "status" => "error",
            "error" => sprint(showerror, e),
            "diagnostics" => _build_job_diagnostics(nothing, nothing, "rejected"; err=e),
        ))
    end
end

function handle_pnjl_scan_job_status(job_id::String)
    pos = _queue_position(job_id)
    snap = _queue_snapshot()
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
            "error" => job["error"],
            "queue" => Dict(
                "position" => pos,
                "queued" => snap.queued,
                "running" => snap.running,
                "max_running" => _PNJL_SCAN_MAX_RUNNING,
                "max_pending" => _PNJL_SCAN_MAX_PENDING,
            ),
            "policy" => get(job, "policy", Dict{String, Any}()),
            "diagnostics" => _build_job_diagnostics(job["job_id"], job["kind"], job["status"]),
        )
    end

    payload === nothing && return _json_response(404, Dict("status" => "error", "error" => "job not found", "job_id" => job_id))
    return _json_response(200, payload)
end

function handle_pnjl_scan_job_result(job_id::String)
    payload = _with_job(job_id) do job
        status = job["status"]
        result = job["result"]
        if status != "succeeded"
            return Dict{String, Any}(
                "status" => "error",
                "job_id" => job_id,
                "error" => "job not succeeded",
                "job_status" => status,
                "diagnostics" => _build_job_diagnostics(job_id, job["kind"], status),
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

    payload === nothing && return _json_response(404, Dict("status" => "error", "error" => "job not found", "job_id" => job_id))
    if payload["status"] == "error"
        return _json_response(409, payload)
    end
    return _json_response(200, payload)
end
