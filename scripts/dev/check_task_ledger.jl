#!/usr/bin/env julia

module TaskLedgerGovernance

using Dates
using TOML

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_LEDGER_REL = joinpath("config", "governance", "task_tracks.toml")
const STATUSES = Set([
    "inbox", "triaged", "ready", "active", "blocked", "review",
    "accepted", "promoted", "archived", "deferred", "cancelled",
])
const CLASSIFICATIONS = Set(["blocker", "required_follow_up", "independent", "research"])
const AUTHOR_REVIEW_STATUSES = Set(["not_required", "pending", "accepted"])
const ACTION_STATUSES = Set(["ready", "active", "blocked", "review", "accepted", "promoted"])
const TOP_LEVEL_FIELDS = Set(["schema_version", "primary_track", "updated_at", "tracks", "items"])
const TRACK_FIELDS = Set([
    "id", "title", "status", "current_task", "blocked_by", "unlocks", "next_action",
    "required_author_review", "author_review_status", "promotion_required", "current_branch",
    "current_sha", "run_ids", "evidence", "updated_at",
])
const ITEM_FIELDS = Set([
    "id", "track_id", "classification", "status", "parent", "task_file", "reason",
    "next_action", "blocked_by", "backlog_file", "evidence",
])
const SHA_RE = r"^[0-9a-f]{40}$"
const BRANCH_RE = r"^[A-Za-z0-9][A-Za-z0-9._/-]*$"
const RUN_RE = r"^[0-9]+$"
const ID_RE = r"^[A-Za-z0-9][A-Za-z0-9_-]*$"
const REF_RE = r"^(track|item):([A-Za-z0-9][A-Za-z0-9_-]*)$"

const ALLOWED_TRANSITIONS = Dict(
    "inbox" => Set(["triaged"]),
    "triaged" => Set(["ready", "deferred", "cancelled"]),
    "ready" => Set(["active", "deferred", "cancelled"]),
    "active" => Set(["blocked", "review", "deferred", "cancelled"]),
    "blocked" => Set(["ready", "active", "deferred", "cancelled"]),
    "review" => Set(["active", "blocked", "accepted", "cancelled"]),
    "accepted" => Set(["promoted", "archived"]),
    "promoted" => Set(["archived"]),
    "deferred" => Set(["triaged"]),
    "archived" => Set{String}(),
    "cancelled" => Set{String}(),
)

export DEFAULT_LEDGER_REL, STATUSES, CLASSIFICATIONS, ALLOWED_TRANSITIONS
export validate_ledger, validate_transition, summarize_porcelain, preflight_report, main

_label(kind, id, field) = "$(kind)[$(id)].$(field)"

function _string(table, key::String, label::String, violations::Vector{String}; required::Bool=true)
    value = get(table, key, nothing)
    if value === nothing
        required && push!(violations, "$(label): missing required field")
        return ""
    elseif !(value isa AbstractString)
        push!(violations, "$(label): expected string")
        return ""
    end
    value = String(value)
    required && isempty(strip(value)) && push!(violations, "$(label): must not be empty")
    return value
end

function _bool(table, key::String, label::String, violations::Vector{String}; default=nothing)
    value = get(table, key, default)
    if value === nothing
        push!(violations, "$(label): missing required field")
        return false
    elseif !(value isa Bool)
        push!(violations, "$(label): expected boolean")
        return false
    end
    return value
end

function _string_vector(table, key::String, label::String, violations::Vector{String}; default=String[])
    value = get(table, key, default)
    if !(value isa AbstractVector)
        push!(violations, "$(label): expected an array of strings")
        return String[]
    end
    result = String[]
    for (index, item) in enumerate(value)
        if !(item isa AbstractString)
            push!(violations, "$(label)[$index]: expected string")
        else
            push!(result, String(item))
        end
    end
    return result
end

function _table_vector(parsed, key::String, label::String, violations::Vector{String})
    value = get(parsed, key, nothing)
    if value === nothing
        push!(violations, "$(label): missing required array of tables")
        return Vector{Dict{String,Any}}()
    elseif !(value isa AbstractVector)
        push!(violations, "$(label): expected an array of tables")
        return Vector{Dict{String,Any}}()
    end
    tables = Vector{Dict{String,Any}}()
    for (index, item) in enumerate(value)
        if !(item isa AbstractDict)
            push!(violations, "$(label)[$index]: expected table")
        else
            push!(tables, Dict{String,Any}(String(k) => v for (k, v) in item))
        end
    end
    return tables
end

function _validate_keys(table, allowed::Set{String}, label::String, violations::Vector{String})
    keys_seen = Set(String(key) for key in keys(table))
    for key in setdiff(allowed, keys_seen)
        push!(violations, "$(label).$(key): missing required field")
    end
    for key in setdiff(keys_seen, allowed)
        push!(violations, "$(label): unsupported field $(key)")
    end
end

function _valid_date(value::String)
    try
        Date(value, dateformat"yyyy-mm-dd")
        return true
    catch
        return false
    end
end

function _validate_date(value::String, label::String, violations::Vector{String})
    isempty(value) || _valid_date(value) || push!(violations, "$(label): expected YYYY-MM-DD")
end

function _repository_path(root::String, value::String)
    isempty(value) && return nothing
    isabspath(value) && return nothing
    full = normpath(joinpath(root, value))
    relative = relpath(full, root)
    startswith(relative, "..") && return nothing
    return full
end

function _validate_id(value::String, label::String, violations::Vector{String})
    isempty(value) || occursin(ID_RE, value) || push!(violations, "$(label): invalid identifier")
end

function _validate_refs(values::Vector{String}, label::String, violations::Vector{String})
    for value in values
        occursin(REF_RE, value) || push!(violations, "$(label): invalid reference $(value)")
    end
end

function _validate_evidence(root::String, values::Vector{String}, label::String, violations::Vector{String})
    for (index, value) in enumerate(values)
        entry_label = "$(label)[$index]"
        if startswith(value, "file:")
            target = value[6:end]
            isempty(target) && push!(violations, "$(entry_label): empty file reference")
            repository_path = _repository_path(root, target)
            repository_path === nothing && !isempty(target) && push!(violations, "$(entry_label): file reference must stay inside repository")
            repository_path === nothing || ispath(repository_path) || push!(violations, "$(entry_label): file does not exist: $(target)")
        elseif startswith(value, "run:")
            occursin(RUN_RE, value[5:end]) || push!(violations, "$(entry_label): invalid run reference")
        elseif startswith(value, "sha:")
            occursin(SHA_RE, value[5:end]) || push!(violations, "$(entry_label): invalid SHA reference")
        elseif startswith(value, "pr:")
            occursin(RUN_RE, value[4:end]) || push!(violations, "$(entry_label): invalid PR reference")
        elseif startswith(value, "review:") || startswith(value, "promotion:") || startswith(value, "gate:") || startswith(value, "note:")
            isempty(strip(value[findfirst(':', value)+1:end])) && push!(violations, "$(entry_label): empty evidence detail")
        else
            push!(violations, "$(entry_label): unsupported evidence prefix")
        end
    end
end

function _has_prefix(values::Vector{String}, prefixes::Tuple)
    any(value -> any(startswith(value, prefix) for prefix in prefixes), values)
end

function _validate_status_contract(status::String, table, label::String, evidence::Vector{String}, blocked_by::Vector{String}, violations::Vector{String})
    status in STATUSES || push!(violations, "$(label).status: unsupported status $(status)")
    status in ACTION_STATUSES && isempty(strip(String(get(table, "next_action", "")))) && push!(violations, "$(label).next_action: required for $(status)")
    status == "blocked" && isempty(blocked_by) && push!(violations, "$(label).blocked_by: blocked state requires at least one dependency")
    status in ("accepted", "promoted") && isempty(evidence) && push!(violations, "$(label).evidence: required for $(status)")
    status == "promoted" && !_has_prefix(evidence, ("promotion:", "gate:")) && push!(violations, "$(label).evidence: promoted state requires promotion:/gate: evidence")
end

function _dependency_graph(tracks, items, track_ids, item_ids, violations)
    nodes = Set{String}()
    graph = Dict{String,Vector{String}}()
    for id in track_ids
        node = "track:" * id
        push!(nodes, node)
        graph[node] = String[]
    end
    for id in item_ids
        node = "item:" * id
        push!(nodes, node)
        graph[node] = String[]
    end

    function add_edges(table, node, field, label)
        refs = _string_vector(table, field, label, violations)
        for ref in refs
            ref in nodes || push!(violations, "$(label): unknown reference $(ref)")
            push!(get!(graph, node, String[]), ref)
        end
    end

    function check_refs(table, field, label)
        refs = _string_vector(table, field, label, violations)
        for ref in refs
            ref in nodes || push!(violations, "$(label): unknown reference $(ref)")
        end
    end

    for (id, table) in tracks
        add_edges(table, "track:" * id, "blocked_by", _label("track", id, "blocked_by"))
        check_refs(table, "unlocks", _label("track", id, "unlocks"))
    end
    for (id, table) in items
        add_edges(table, "item:" * id, "blocked_by", _label("item", id, "blocked_by"))
        parent = String(get(table, "parent", ""))
        if !isempty(parent) && !(parent in nodes)
            push!(violations, "$(_label("item", id, "parent")): unknown reference $(parent)")
        end
    end

    colors = Dict(node => 0 for node in keys(graph))
    function visit(node)
        color = get(colors, node, 0)
        color == 1 && return true
        color == 2 && return false
        colors[node] = 1
        for dependency in get(graph, node, String[])
            visit(dependency) && return true
        end
        colors[node] = 2
        return false
    end
    for node in keys(graph)
        if visit(node)
            push!(violations, "blocked_by dependency graph contains a cycle")
            break
        end
    end
    return graph
end

function _parse_ledger(root::String, ledger_rel::String, violations::Vector{String})
    path = isabspath(ledger_rel) ? ledger_rel : joinpath(root, ledger_rel)
    isfile(path) || begin
        push!(violations, "ledger file not found: $(path)")
        return nothing
    end
    try
        return TOML.parsefile(path)
    catch err
        push!(violations, "ledger TOML parse failed: $(sprint(showerror, err))")
        return nothing
    end
end

function validate_ledger(root::AbstractString=ROOT; ledger_rel::AbstractString=DEFAULT_LEDGER_REL)
    root = normpath(abspath(String(root)))
    violations = String[]
    parsed = _parse_ledger(root, String(ledger_rel), violations)
    parsed === nothing && return violations

    _validate_keys(parsed, TOP_LEVEL_FIELDS, "ledger", violations)
    get(parsed, "schema_version", nothing) == 1 || push!(violations, "schema_version must be integer 1")
    primary_track = _string(parsed, "primary_track", "ledger", violations)
    updated_at = _string(parsed, "updated_at", "ledger", violations)
    _validate_date(updated_at, "ledger.updated_at", violations)

    track_tables = _table_vector(parsed, "tracks", "tracks", violations)
    item_tables = _table_vector(parsed, "items", "items", violations)
    tracks = Dict{String,Dict{String,Any}}()
    items = Dict{String,Dict{String,Any}}()
    task_files = Dict{String,String}()

    for table in track_tables
        _validate_keys(table, TRACK_FIELDS, "track", violations)
        id = _string(table, "id", "track", violations; required=false)
        _validate_id(id, "track.id", violations)
        isempty(id) && continue
        haskey(tracks, id) && push!(violations, "duplicate track id: $(id)")
        tracks[id] = table

        label = "track[$id]"
        status = _string(table, "status", label, violations)
        current_task = _string(table, "current_task", label, violations; required=false)
        blocked_by = _string_vector(table, "blocked_by", "$label.blocked_by", violations)
        unlocks = _string_vector(table, "unlocks", "$label.unlocks", violations)
        evidence = _string_vector(table, "evidence", "$label.evidence", violations)
        _validate_status_contract(status, table, label, evidence, blocked_by, violations)
        _validate_refs(blocked_by, "$label.blocked_by", violations)
        _validate_refs(unlocks, "$label.unlocks", violations)
        required_review = _bool(table, "required_author_review", "$label.required_author_review", violations)
        review_status = _string(table, "author_review_status", label, violations)
        review_status in AUTHOR_REVIEW_STATUSES || push!(violations, "$label.author_review_status: unsupported value $(review_status)")
        !required_review && review_status != "not_required" && push!(violations, "$label.author_review_status: must be not_required when review is not required")
        promotion_required = _bool(table, "promotion_required", "$label.promotion_required", violations)
        status == "accepted" && promotion_required &&
            push!(violations, "$label: accepted -> archived is only valid when promotion_required=false")
        status == "promoted" && !promotion_required &&
            push!(violations, "$label: promoted state requires promotion_required=true")
        branch = _string(table, "current_branch", label, violations; required=false)
        isempty(branch) || occursin(BRANCH_RE, branch) || push!(violations, "$label.current_branch: invalid branch name")
        sha = _string(table, "current_sha", label, violations; required=false)
        isempty(sha) || occursin(SHA_RE, sha) || push!(violations, "$label.current_sha: expected lowercase 40-character SHA or empty")
        run_ids = _string_vector(table, "run_ids", "$label.run_ids", violations)
        for (index, run_id) in enumerate(run_ids)
            occursin(RUN_RE, run_id) || push!(violations, "$label.run_ids[$index]: invalid run id")
        end
        track_updated = _string(table, "updated_at", label, violations)
        _validate_date(track_updated, "$label.updated_at", violations)
        if status in ACTION_STATUSES && isempty(current_task)
            push!(violations, "$label.current_task: required for $(status)")
        end
    end

    primary_track in keys(tracks) || push!(violations, "ledger.primary_track does not reference a track: $(primary_track)")
    if haskey(tracks, primary_track)
        primary_status = String(get(tracks[primary_track], "status", ""))
        primary_status in ("archived", "cancelled") && push!(violations, "ledger.primary_track cannot be terminal")
    end

    for table in item_tables
        _validate_keys(table, ITEM_FIELDS, "item", violations)
        id = _string(table, "id", "item", violations; required=false)
        _validate_id(id, "item.id", violations)
        isempty(id) && continue
        haskey(items, id) && push!(violations, "duplicate item id: $(id)")
        items[id] = table

        label = "item[$id]"
        track_id = _string(table, "track_id", label, violations)
        track_id in keys(tracks) || push!(violations, "$label.track_id: unknown track $(track_id)")
        classification = _string(table, "classification", label, violations)
        classification in CLASSIFICATIONS || push!(violations, "$label.classification: unsupported value $(classification)")
        status = _string(table, "status", label, violations)
        parent = _string(table, "parent", label, violations; required=false)
        task_file = _string(table, "task_file", label, violations; required=false)
        reason = _string(table, "reason", label, violations)
        next_action = _string(table, "next_action", label, violations; required=false)
        blocked_by = _string_vector(table, "blocked_by", "$label.blocked_by", violations)
        backlog_file = _string(table, "backlog_file", label, violations; required=false)
        evidence = _string_vector(table, "evidence", "$label.evidence", violations)
        _validate_status_contract(status, table, label, evidence, blocked_by, violations)
        _validate_refs(blocked_by, "$label.blocked_by", violations)
        isempty(parent) || occursin(REF_RE, parent) || push!(violations, "$label.parent: invalid reference $(parent)")
        status in ACTION_STATUSES && isempty(next_action) && push!(violations, "$label.next_action: required for $(status)")
        status == "deferred" && (isempty(parent) || isempty(reason) || isempty(next_action) || isempty(backlog_file)) &&
            push!(violations, "$label: deferred items require parent, reason, next_action and backlog_file")
        status in ACTION_STATUSES && isempty(task_file) && push!(violations, "$label.task_file: required for $(status)")
        if !isempty(task_file)
            repository_path = _repository_path(root, task_file)
            repository_path === nothing && push!(violations, "$label.task_file: must stay inside repository")
            repository_path === nothing || isfile(repository_path) || push!(violations, "$label.task_file: file does not exist: $(task_file)")
            if haskey(task_files, task_file) && task_files[task_file] != id
                push!(violations, "task file belongs to multiple items: $(task_file)")
            else
                task_files[task_file] = id
            end
        end
        if !isempty(backlog_file)
            repository_path = _repository_path(root, backlog_file)
            repository_path === nothing && push!(violations, "$label.backlog_file: must stay inside repository")
            repository_path === nothing || isfile(repository_path) || push!(violations, "$label.backlog_file: file does not exist: $(backlog_file)")
        end
        _validate_evidence(root, evidence, "$label.evidence", violations)
    end

    for (id, table) in tracks
        current_task = String(get(table, "current_task", ""))
        isempty(current_task) && continue
        haskey(items, current_task) || begin
            push!(violations, "track[$id].current_task: unknown item $(current_task)")
            continue
        end
        item_track = String(get(items[current_task], "track_id", ""))
        item_track == id || push!(violations, "track[$id].current_task $(current_task) belongs to $(item_track)")
        item_status = String(get(items[current_task], "status", ""))
        track_status = String(get(table, "status", ""))
        item_status == track_status || push!(violations, "track[$id] status $(track_status) disagrees with current item $(current_task) status $(item_status)")
    end

    _dependency_graph(
        [(id, table) for (id, table) in tracks],
        [(id, table) for (id, table) in items],
        keys(tracks),
        keys(items),
        violations,
    )
    for (id, table) in tracks
        evidence = _string_vector(table, "evidence", "track[$id].evidence", violations)
        _validate_evidence(root, evidence, "track[$id].evidence", violations)
        required_review = Bool(get(table, "required_author_review", false))
        review_status = String(get(table, "author_review_status", ""))
        status = String(get(table, "status", ""))
        status == "promoted" && required_review && review_status != "accepted" && push!(violations, "track[$id]: promoted requires accepted author review")
    end
    return unique(violations)
end

function validate_transition(from::AbstractString, to::AbstractString; promotion_required::Bool=false)
    from in STATUSES || return false
    to in get(ALLOWED_TRANSITIONS, String(from), Set{String}()) || return false
    from == "accepted" && to == "archived" && promotion_required && return false
    return true
end

function summarize_porcelain(status_text::AbstractString)
    lines = filter(!isempty, split(replace(String(status_text), "\r\n" => "\n", "\r" => "\n"), '\n'))
    paths = String[]
    untracked = 0
    tracked = 0
    for line in lines
        if startswith(line, "?? ")
            untracked += 1
            push!(paths, line[4:end])
        else
            tracked += 1
            push!(paths, length(line) > 3 ? line[4:end] : line)
        end
    end
    return (; dirty=!isempty(lines), tracked, untracked, paths)
end

function _git_output(root::String, args...)
    try
        return readchomp(Cmd(["git", "-C", root, args...]))
    catch err
        return "<git unavailable: $(sprint(showerror, err))>"
    end
end

function preflight_report(root::AbstractString=ROOT; ledger_rel::AbstractString=DEFAULT_LEDGER_REL, track_id=nothing, io::IO=stdout)
    root = normpath(abspath(String(root)))
    violations = validate_ledger(root; ledger_rel=ledger_rel)
    isempty(violations) || return violations
    parsed = TOML.parsefile(isabspath(ledger_rel) ? String(ledger_rel) : joinpath(root, ledger_rel))
    primary_track = String(parsed["primary_track"])
    selected = track_id === nothing ? primary_track : String(track_id)
    tracks = get(parsed, "tracks", Any[])
    selected_index = findfirst(table -> String(get(table, "id", "")) == selected, tracks)
    selected_index === nothing && return ["preflight track not found: $(selected)"]
    selected_table = tracks[selected_index]
    branch = _git_output(root, "branch", "--show-current")
    head = _git_output(root, "rev-parse", "HEAD")
    porcelain = _git_output(root, "status", "--porcelain=v1")
    summary = summarize_porcelain(porcelain)
    selected_status = String(get(selected_table, "status", ""))
    selected_task = String(get(selected_table, "current_task", ""))
    println(io, "[task-ledger-preflight] track=$(selected)")
    println(io, "  primary_track=$(primary_track)")
    println(io, "  selected_track=$(selected)")
    println(io, "  status=$(selected_status)")
    println(io, "  current_task=$(selected_task)")
    println(io, "  branch=$(branch)")
    println(io, "  head=$(head)")
    println(io, "  dirty=$(summary.dirty) tracked=$(summary.tracked) untracked=$(summary.untracked)")
    if summary.dirty
        println(io, "  ATTENTION: preserve these paths before changing track state:")
        for path in summary.paths
            println(io, "    - $(path)")
        end
    end
    return String[]
end

function _track_selection_violation(root::String, ledger_rel::String, track_id)
    track_id === nothing && return nothing
    ledger_path = isabspath(ledger_rel) ? ledger_rel : joinpath(root, ledger_rel)
    parsed = TOML.parsefile(ledger_path)
    tracks = get(parsed, "tracks", Any[])
    any(table -> String(get(table, "id", "")) == String(track_id), tracks) ||
        return "track not found: $(track_id)"
    return nothing
end

function _parse_args(args)
    ledger_rel = DEFAULT_LEDGER_REL
    preflight = false
    track_id = nothing
    i = 1
    while i <= length(args)
        arg = String(args[i])
        if arg == "--ledger"
            i == length(args) && error("--ledger requires a path")
            i += 1
            ledger_rel = String(args[i])
        elseif startswith(arg, "--ledger=")
            ledger_rel = split(arg, '='; limit=2)[2]
        elseif arg == "--preflight"
            preflight = true
        elseif arg == "--track"
            i == length(args) && error("--track requires an id")
            i += 1
            track_id = String(args[i])
        elseif startswith(arg, "--track=")
            track_id = split(arg, '='; limit=2)[2]
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=. scripts/dev/check_task_ledger.jl [--ledger PATH] [--preflight] [--track ID]")
            return nothing
        else
            error("unknown option: $(arg)")
        end
        i += 1
    end
    return (; ledger_rel, preflight, track_id)
end

function main(args::Vector{String}=collect(String.(ARGS)))
    options = _parse_args(args)
    options === nothing && return 0
    violations = if options.preflight
        preflight_report(ROOT; ledger_rel=options.ledger_rel, track_id=options.track_id)
    else
        result = validate_ledger(ROOT; ledger_rel=options.ledger_rel)
        if isempty(result)
            selection_violation = _track_selection_violation(ROOT, String(options.ledger_rel), options.track_id)
            selection_violation === nothing || push!(result, selection_violation)
        end
        result
    end
    if !isempty(violations)
        println("[task-ledger] FAILED")
        for violation in violations
            println(" - " * violation)
        end
        return 1
    end
    if !options.preflight
        suffix = options.track_id === nothing ? "" : " track=$(options.track_id)"
        println("[task-ledger] OK ledger=$(replace(String(options.ledger_rel), '\\' => '/'))$(suffix)")
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end

end
