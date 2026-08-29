using Test
using TOML

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const LEDGER_CHECKER = joinpath(PROJECT_ROOT, "scripts", "dev", "check_task_ledger.jl")
if !isdefined(Main, :TaskLedgerGovernance)
    include(LEDGER_CHECKER)
end
const TL = Main.TaskLedgerGovernance

function _valid_ledger(; track_status="active", item_status=track_status, blocked_by="[]", evidence="[\"file:docs/dev/active/task.md\"]", current_sha="aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa", current_branch="codex/fixture", run_ids="[\"12345\"]", task_file="docs/dev/active/task.md", backlog_file="", required_author_review=false, author_review_status="not_required", promotion_required=false, parent="track:track-a", reason="Fixture reason", next_action="Continue the fixture task", extra="")
    return """
schema_version = 1
primary_track = "track-a"
updated_at = "2026-08-15"

[[tracks]]
id = "track-a"
title = "Fixture track"
status = "$(track_status)"
current_task = "item-a"
blocked_by = $(blocked_by)
unlocks = []
next_action = "Continue the fixture task"
required_author_review = $(required_author_review)
author_review_status = "$(author_review_status)"
promotion_required = $(promotion_required)
current_branch = "$(current_branch)"
current_sha = "$(current_sha)"
run_ids = $(run_ids)
evidence = $(evidence)
updated_at = "2026-08-15"

[[items]]
id = "item-a"
track_id = "track-a"
classification = "required_follow_up"
status = "$(item_status)"
parent = "$(parent)"
task_file = "$(task_file)"
reason = "$(reason)"
next_action = "$(next_action)"
blocked_by = $(blocked_by)
backlog_file = "$(backlog_file)"
evidence = $(evidence)
$(extra)
"""
end

function _fixture_root(; kwargs...)
    root = mktempdir()
    mkpath(joinpath(root, "docs", "dev", "active"))
    write(joinpath(root, "docs", "dev", "active", "task.md"), "# Fixture task\n")
    if haskey(kwargs, :backlog_file) && !isempty(kwargs[:backlog_file])
        mkpath(joinpath(root, "docs", "dev", "backlog"))
        write(joinpath(root, kwargs[:backlog_file]), "# Fixture backlog\n")
    end
    mkpath(joinpath(root, "config", "governance"))
    write(joinpath(root, "config", "governance", "task_tracks.toml"), _valid_ledger(; kwargs...))
    return root
end

@testset "task ledger real state" begin
    violations = TL.validate_ledger(PROJECT_ROOT)
    @test isempty(violations)
    parsed = TOML.parsefile(joinpath(PROJECT_ROOT, "config", "governance", "task_tracks.toml"))
    tracks = Dict(String(t["id"]) => t for t in parsed["tracks"])
    @test parsed["primary_track"] == "issue130-phase"
    @test tracks["issue130-phase"]["status"] == "active"
    @test tracks["issue130-phase"]["current_task"] == "issue130-phase-reference-accepted-runtime-fallback"
    @test tracks["rs-transport"]["status"] == "active"
    @test isempty(tracks["rs-transport"]["blocked_by"])
    @test tracks["plot-sop"]["status"] == "triaged"
    items = Dict(String(item["id"]) => item for item in parsed["items"])
    @test items["issue130-phase-reference-retirement"]["status"] == "accepted"
    @test items["issue130-phase-reference-retirement"]["classification"] == "required_follow_up"
    @test isempty(items["issue130-phase-reference-retirement"]["blocked_by"])
    @test occursin("full_hybrid_candidate", read(joinpath(PROJECT_ROOT, "docs", "dev", "task_tracking_governance.md"), String))
    @test !occursin("status = \"full_hybrid_candidate\"", read(joinpath(PROJECT_ROOT, "config", "governance", "task_tracks.toml"), String))
end

@testset "task ledger state transitions" begin
    for (from, destinations) in TL.ALLOWED_TRANSITIONS
        for to in destinations
            @test TL.validate_transition(from, to)
        end
        forbidden = setdiff(TL.STATUSES, union(destinations, Set([from])))
        for to in forbidden
            @test !TL.validate_transition(from, to)
        end
    end
    @test TL.validate_transition("accepted", "archived"; promotion_required=false)
    @test !TL.validate_transition("accepted", "archived"; promotion_required=true)
    @test !TL.validate_transition("unknown", "active")
    @test TL.CLASSIFICATIONS == Set(["blocker", "required_follow_up", "independent", "research"])
end

@testset "task ledger valid fixture" begin
    root = _fixture_root()
    @test isempty(TL.validate_ledger(root))
end

@testset "task ledger schema and triage classes" begin
    root = _fixture_root(backlog_file="docs/dev/backlog/research.md", extra="""
[[items]]
id = "item-blocker"
track_id = "track-a"
classification = "blocker"
status = "blocked"
parent = "track:track-a"
task_file = "docs/dev/active/blocker.md"
reason = "Blocking fixture"
next_action = "Resolve the dependency"
blocked_by = ["track:track-a"]
backlog_file = ""
evidence = []

[[items]]
id = "item-independent"
track_id = "track-a"
classification = "independent"
status = "triaged"
parent = "track:track-a"
task_file = ""
reason = "Independent fixture"
next_action = "Create an isolated branch"
blocked_by = []
backlog_file = ""
evidence = []

[[items]]
id = "item-research"
track_id = "track-a"
classification = "research"
status = "deferred"
parent = "track:track-a"
task_file = ""
reason = "Research fixture"
next_action = "Record the question in backlog"
blocked_by = []
backlog_file = "docs/dev/backlog/research.md"
evidence = []
""")
    write(joinpath(root, "docs", "dev", "active", "blocker.md"), "# Blocking fixture\n")
    parsed = TOML.parsefile(joinpath(root, "config", "governance", "task_tracks.toml"))
    @test Set(String(item["classification"]) for item in parsed["items"]) == TL.CLASSIFICATIONS
    @test isempty(TL.validate_ledger(root))

    root = _fixture_root()
    ledger_path = joinpath(root, "config", "governance", "task_tracks.toml")
    ledger = replace(read(ledger_path, String), "schema_version = 1\n" => "schema_version = 1\nunexpected = true\n")
    write(ledger_path, ledger)
    violations = TL.validate_ledger(root)
    @test any(occursin("unsupported field unexpected", item) for item in violations)

    root = _fixture_root()
    ledger_path = joinpath(root, "config", "governance", "task_tracks.toml")
    ledger = replace(read(ledger_path, String), "current_sha = \"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\"\n" => "")
    write(ledger_path, ledger)
    violations = TL.validate_ledger(root)
    @test any(occursin("current_sha: missing required field", item) for item in violations)
end

@testset "task ledger rejects missing contracts" begin
    root = _fixture_root(track_status="blocked", item_status="blocked", blocked_by="[]")
    violations = TL.validate_ledger(root)
    @test any(occursin("blocked state requires", item) for item in violations)

    root = _fixture_root(track_status="accepted", item_status="accepted", evidence="[]")
    violations = TL.validate_ledger(root)
    @test any(occursin("evidence: required for accepted", item) for item in violations)

    root = _fixture_root(track_status="accepted", item_status="accepted", promotion_required=true)
    violations = TL.validate_ledger(root)
    @test any(occursin("accepted -> archived is only valid", item) for item in violations)

    root = _fixture_root(track_status="active", item_status="active", current_sha="bad-sha")
    violations = TL.validate_ledger(root)
    @test any(occursin("current_sha: expected", item) for item in violations)

    root = _fixture_root(track_status="active", item_status="active", current_branch="bad branch", run_ids="[\"not-a-run\"]")
    violations = TL.validate_ledger(root)
    @test any(occursin("current_branch: invalid", item) for item in violations)
    @test any(occursin("run_ids[1]: invalid", item) for item in violations)

    root = _fixture_root(track_status="active", item_status="active", task_file="docs/dev/active/missing.md")
    violations = TL.validate_ledger(root)
    @test any(occursin("task_file: file does not exist", item) for item in violations)

    root = _fixture_root(track_status="promoted", item_status="promoted", evidence="[\"file:docs/dev/active/task.md\"]")
    violations = TL.validate_ledger(root)
    @test any(occursin("promoted state requires promotion:/gate:", item) for item in violations)

    root = _fixture_root(track_status="promoted", item_status="promoted", evidence="[\"promotion:phase-reference\"]")
    violations = TL.validate_ledger(root)
    @test any(occursin("promoted state requires promotion_required=true", item) for item in violations)

    root = _fixture_root(
        track_status="promoted",
        item_status="promoted",
        evidence="[\"promotion:phase-reference\"]",
        required_author_review=true,
        author_review_status="pending",
        promotion_required=true,
    )
    violations = TL.validate_ledger(root)
    @test any(occursin("promoted requires accepted author review", item) for item in violations)

    root = _fixture_root(track_status="deferred", item_status="deferred")
    violations = TL.validate_ledger(root)
    @test any(occursin("deferred items require", item) for item in violations)

    root = _fixture_root(track_status="deferred", item_status="deferred", parent="")
    violations = TL.validate_ledger(root)
    @test any(occursin("deferred items require", item) for item in violations)

    root = _fixture_root(track_status="deferred", item_status="deferred", backlog_file="")
    violations = TL.validate_ledger(root)
    @test any(occursin("deferred items require", item) for item in violations)

    root = _fixture_root(track_status="deferred", item_status="deferred", backlog_file="docs/dev/backlog/task.md")
    @test isempty(TL.validate_ledger(root))
end

@testset "task ledger rejects duplicate, dangling and cyclic references" begin
    root = _fixture_root(extra="""
[[items]]
id = "item-b"
track_id = "track-a"
classification = "research"
status = "triaged"
parent = "track:track-a"
task_file = ""
reason = "Research fixture"
next_action = "Place in backlog"
blocked_by = ["item:item-c"]
backlog_file = ""
evidence = []

[[items]]
id = "item-c"
track_id = "track-a"
classification = "blocker"
status = "blocked"
parent = "track:track-a"
task_file = "docs/dev/active/task.md"
reason = "Cyclic fixture"
next_action = "Resolve cycle"
blocked_by = ["item:item-b"]
backlog_file = ""
evidence = []
""")
    violations = TL.validate_ledger(root)
    @test any(occursin("task file belongs to multiple items", item) for item in violations)
    @test any(occursin("dependency graph contains a cycle", item) for item in violations)

    root = _fixture_root(extra="""
[[tracks]]
id = "track-b"
title = "Second fixture track"
status = "active"
current_task = "item-b"
blocked_by = []
unlocks = []
next_action = "Continue the second fixture task"
required_author_review = false
author_review_status = "not_required"
promotion_required = false
current_branch = "codex/fixture-b"
current_sha = "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"
run_ids = ["12346"]
evidence = ["file:docs/dev/active/task.md"]
updated_at = "2026-08-15"

[[items]]
id = "item-b"
track_id = "track-b"
classification = "independent"
status = "active"
parent = "track:track-b"
task_file = "docs/dev/active/task.md"
reason = "Cross-track duplicate fixture"
next_action = "Resolve ownership"
blocked_by = []
backlog_file = ""
evidence = ["file:docs/dev/active/task.md"]
""")
    violations = TL.validate_ledger(root)
    @test any(occursin("task file belongs to multiple items", item) for item in violations)

    root = _fixture_root(blocked_by="[\"track:missing\"]")
    violations = TL.validate_ledger(root)
    @test any(occursin("unknown reference track:missing", item) for item in violations)

    root = _fixture_root()
    ledger_path = joinpath(root, "config", "governance", "task_tracks.toml")
    ledger = replace(read(ledger_path, String), "unlocks = []" => "unlocks = [\"track:missing\"]")
    write(ledger_path, ledger)
    violations = TL.validate_ledger(root)
    @test any(occursin("unlocks: unknown reference track:missing", item) for item in violations)

    root = _fixture_root(task_file="../outside.md")
    violations = TL.validate_ledger(root)
    @test any(occursin("task_file: must stay inside repository", item) for item in violations)
end

@testset "task ledger preflight porcelain fixture" begin
    summary = TL.summarize_porcelain(" M tracked.md\n?? new.md\n")
    @test summary.dirty
    @test summary.tracked == 1
    @test summary.untracked == 1
    @test summary.paths == ["tracked.md", "new.md"]
    clean = TL.summarize_porcelain("")
    @test !clean.dirty
    @test isempty(clean.paths)

    dirty_git = (root, args...) -> begin
        args == ("branch", "--show-current") && return "codex/fixture"
        args == ("rev-parse", "HEAD") && return "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
        args == ("status", "--porcelain=v1") && return " M tracked.md\n?? new.md\n"
        error("unexpected git fixture command: $(args)")
    end
    output = IOBuffer()
    @test isempty(TL.preflight_report(PROJECT_ROOT; track_id="rs-transport", io=output, git_output=dirty_git))
    report = String(take!(output))
    @test occursin("primary_track=issue130-phase", report)
    @test occursin("selected_track=rs-transport", report)
    @test occursin("branch=codex/fixture", report)
    @test occursin("head=aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa", report)
    @test occursin("dirty=true", report)
    @test occursin("ATTENTION: preserve", report)

    clean_git = (root, args...) -> begin
        args == ("branch", "--show-current") && return ""
        args == ("rev-parse", "HEAD") && return "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"
        args == ("status", "--porcelain=v1") && return ""
        error("unexpected git fixture command: $(args)")
    end
    clean_output = IOBuffer()
    @test isempty(TL.preflight_report(PROJECT_ROOT; track_id="rs-transport", io=clean_output, git_output=clean_git))
    clean_report = String(take!(clean_output))
    @test occursin("dirty=false tracked=0 untracked=0", clean_report)
    @test !occursin("ATTENTION: preserve", clean_report)
end

@testset "task ledger harness routing contract" begin
    skill = read(joinpath(PROJECT_ROOT, ".agents", "skills", "codex-task-harness", "SKILL.md"), String)
    agents = read(joinpath(PROJECT_ROOT, "AGENTS.md"), String)
    @test occursin("task_tracks.toml", skill)
    @test occursin("blocker", skill)
    @test occursin("doc-implementation", skill)
    @test occursin("task_tracks.toml", agents)
end

@testset "task ledger CLI track selection" begin
    @test TL.main(["--track", "rs-transport"]) == 0
    @test TL.main(["--track", "missing-track"]) == 1
end
