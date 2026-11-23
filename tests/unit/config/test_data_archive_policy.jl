using Test

if !isdefined(Main, :collect_archive_violations)
    include(joinpath(@__DIR__, "..", "..", "..", "scripts", "dev", "check_data_output_path_guard.jl"))
end

const POLICY_RELATIVE = "config/governance/data_archive_policy.toml"
const ZIP_RELATIVE = "data/outputs/results/demo/archive_v1/demo.zip"

function git_run(root::AbstractString, args::Vararg{String})
    cmd = Cmd(["git", "-C", root, "-c", "user.name=guard-test", "-c", "user.email=guard-test@local", args...])
    return strip(read(cmd, String))
end

function init_repo()
    root = mktempdir()
    git_run(root, "init", "-q")
    mkpath(joinpath(root, "config", "governance"))
    mkpath(dirname(joinpath(root, ZIP_RELATIVE)))
    write(joinpath(root, POLICY_RELATIVE), "schema = \"data_archive_policy_v1\"\n")
    write(joinpath(root, "README.md"), "# demo\n")
    git_run(root, "add", "-A")
    git_run(root, "commit", "-q", "-m", "base")
    return root
end

sha256_of(path) = bytes2hex(SHA.sha256(read(path)))

function write_retained!(root, sha, fields...)
    entry = """
    [[retained_archives]]
    path = "$(ZIP_RELATIVE)"
    sha256 = "$(sha)"
    entry_count = 1
    rationale = "test"
    approved_at = "2026-11-23"
    """
    extra = join(fields, "\n")
    write(joinpath(root, POLICY_RELATIVE), "schema = \"data_archive_policy_v1\"\n\n" * entry * "\n" * extra)
end

function write_retired!(root, sha, blob, before)
    write(joinpath(root, POLICY_RELATIVE), """
    schema = "data_archive_policy_v1"

    [[retired_archives]]
    path = "$(ZIP_RELATIVE)"
    sha256 = "$(sha)"
    entry_count = 1
    declared_commit = "$(before)"
    uncommitted_source = false
    snapshot_mtime = "2026-10-26"
    replacement = "git cat-file blob $(blob)"
    rationale = "test"
    retired_at = "2026-11-23"
    git_blob = "$(blob)"
    retired_before_commit = "$(before)"
    """)
end

@testset "Unregistered new archive is rejected" begin
    root = init_repo()
    base = git_run(root, "rev-parse", "HEAD")
    write(joinpath(root, ZIP_RELATIVE), "payload")
    git_run(root, "add", "-A")
    git_run(root, "commit", "-q", "-m", "add archive")
    tip = git_run(root, "rev-parse", "HEAD")

    policy = load_archive_policy(root)
    violations = collect_archive_violations(root, GuardConfig(:range, base, tip), policy)
    @test any(occursin("not registered", v) for v in violations)
end

@testset "Retained archive must match its recorded sha256" begin
    root = init_repo()
    write(joinpath(root, ZIP_RELATIVE), "payload")
    git_run(root, "add", "-A")
    git_run(root, "commit", "-q", "-m", "add archive")
    base = git_run(root, "rev-parse", "HEAD~1")
    tip = git_run(root, "rev-parse", "HEAD")
    digest = sha256_of(joinpath(root, ZIP_RELATIVE))

    write_retained!(root, digest)
    git_run(root, "add", "-A")
    git_run(root, "commit", "-q", "-m", "register archive")
    tip = git_run(root, "rev-parse", "HEAD")
    violations = collect_archive_violations(root, GuardConfig(:range, base, tip), load_archive_policy(root))
    @test isempty(violations)

    write_retained!(root, repeat("0", 64))
    violations = collect_archive_violations(root, GuardConfig(:range, base, tip), load_archive_policy(root))
    @test any(occursin("sha256 mismatch", v) for v in violations)

    write_retained!(root, "not-a-hash")
    violations = collect_archive_violations(root, GuardConfig(:range, base, tip), load_archive_policy(root))
    @test any(occursin("64 lowercase hex", v) for v in violations)
end

@testset "Deleting an archive requires a tombstone" begin
    root = init_repo()
    write(joinpath(root, ZIP_RELATIVE), "payload")
    git_run(root, "add", "-A")
    git_run(root, "commit", "-q", "-m", "add archive")
    digest = sha256_of(joinpath(root, ZIP_RELATIVE))
    before = git_run(root, "rev-parse", "HEAD")
    blob = git_run(root, "rev-parse", "HEAD:$(ZIP_RELATIVE)")

    git_run(root, "rm", "-q", ZIP_RELATIVE)
    git_run(root, "commit", "-q", "-m", "drop archive")
    tip = git_run(root, "rev-parse", "HEAD")

    violations = collect_archive_violations(root, GuardConfig(:range, before, tip), load_archive_policy(root))
    @test any(occursin("requires a retired_archives tombstone", v) for v in violations)

    write_retired!(root, digest, blob, before)
    git_run(root, "add", "-A")
    git_run(root, "commit", "-q", "-m", "tombstone")
    tip = git_run(root, "rev-parse", "HEAD")
    violations = collect_archive_violations(root, GuardConfig(:range, before, tip), load_archive_policy(root))
    @test isempty(violations)
end

@testset "Retired archive must not remain in the working tree" begin
    root = init_repo()
    write(joinpath(root, ZIP_RELATIVE), "payload")
    git_run(root, "add", "-A")
    git_run(root, "commit", "-q", "-m", "add archive")
    digest = sha256_of(joinpath(root, ZIP_RELATIVE))
    before = git_run(root, "rev-parse", "HEAD")
    blob = git_run(root, "rev-parse", "HEAD:$(ZIP_RELATIVE)")
    write_retired!(root, digest, blob, before)

    cfg = GuardConfig(:working, nothing, nothing)
    violations = collect_archive_violations(root, cfg, load_archive_policy(root))
    @test any(occursin("still present in the working tree", v) for v in violations)

    git_run(root, "rm", "-q", ZIP_RELATIVE)
    violations = collect_archive_violations(root, cfg, load_archive_policy(root))
    @test isempty(violations)
end

@testset "Unregistered legacy archives are a hard failure" begin
    root = init_repo()
    write(joinpath(root, ZIP_RELATIVE), "payload")
    git_run(root, "add", "-A")
    git_run(root, "commit", "-q", "-m", "add archive")

    policy = load_archive_policy(root)
    violations = collect_archive_violations(root, GuardConfig(:working, nothing, nothing), policy)
    @test any(occursin("is not registered in", v) for v in violations)
    @test isempty(collect_archive_advisories(root, policy))

    write_retained!(root, sha256_of(joinpath(root, ZIP_RELATIVE)))
    @test isempty(collect_archive_violations(root, GuardConfig(:working, nothing, nothing), load_archive_policy(root)))
end

@testset "Retained status is restricted to the documented tiers" begin
    root = init_repo()
    write(joinpath(root, ZIP_RELATIVE), "payload")
    git_run(root, "add", "-A")
    git_run(root, "commit", "-q", "-m", "add archive")
    digest = sha256_of(joinpath(root, ZIP_RELATIVE))
    cfg = GuardConfig(:working, nothing, nothing)

    for status in ("active", "historical_superseded", "historical_evidence")
        write_retained!(root, digest, "status = \"$(status)\"")
        @test isempty(collect_archive_violations(root, cfg, load_archive_policy(root)))
    end

    write_retained!(root, digest, "status = \"current_data\"")
    violations = collect_archive_violations(root, cfg, load_archive_policy(root))
    @test any(occursin("status must be one of", v) for v in violations)
end

@testset "Tombstone git_blob must still hash to its recorded sha256" begin
    root = init_repo()
    write(joinpath(root, ZIP_RELATIVE), "payload")
    git_run(root, "add", "-A")
    git_run(root, "commit", "-q", "-m", "add archive")
    digest = sha256_of(joinpath(root, ZIP_RELATIVE))
    before = git_run(root, "rev-parse", "HEAD")
    blob = git_run(root, "rev-parse", "HEAD:$(ZIP_RELATIVE)")
    git_run(root, "rm", "-q", ZIP_RELATIVE)
    git_run(root, "commit", "-q", "-m", "drop archive")
    tip = git_run(root, "rev-parse", "HEAD")
    cfg = GuardConfig(:range, before, tip)

    write_retired!(root, digest, blob, before)
    @test isempty(collect_archive_violations(root, cfg, load_archive_policy(root)))

    write_retired!(root, repeat("a", 64), blob, before)
    violations = collect_archive_violations(root, cfg, load_archive_policy(root))
    @test any(occursin("hashes to", v) for v in violations)

    write_retired!(root, digest, repeat("f", 40), before)
    violations = collect_archive_violations(root, cfg, load_archive_policy(root))
    @test any(occursin("missing from the local object database", v) for v in violations)
end

@testset "Policy schema is validated" begin
    root = init_repo()
    write(joinpath(root, POLICY_RELATIVE), "schema = \"wrong\"\n")
    @test_throws ErrorException load_archive_policy(root)
    @test load_archive_policy(joinpath(root, "missing")) === nothing
end
