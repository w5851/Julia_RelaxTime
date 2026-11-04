include(joinpath(@__DIR__, "source", "scripts", "dev", "check_task_ledger.jl"))

function snapshot_git(_root, args...)
    names = Dict("branch --show-current" => "branch.txt",
                 "rev-parse HEAD" => "head.txt",
                 "status --porcelain=v1" => "status.txt")
    return chomp(read(joinpath(@__DIR__, "local_git", names[join(args, " ")]), String))
end

root = joinpath(@__DIR__, "source")
println("Git context below is the frozen local workspace snapshot; checks execute on W.")
for ledger in ("task_tracks_working.toml", "task_tracks.toml")
    println("Checking ", ledger)
    violations = TaskLedgerGovernance.preflight_report(root;
        ledger_rel=joinpath("config", "governance", ledger),
        track_id="transport-threshold-upgrade", git_output=snapshot_git)
    isempty(violations) || error(join(violations, "\n"))
end
println("Both working and selected-commit task ledgers passed.")
