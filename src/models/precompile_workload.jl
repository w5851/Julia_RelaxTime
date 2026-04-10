"""Representative compile workload for AD/solver hot paths.

This is a lightweight warmup entry used by tests and scan scripts to reduce
first-call latency in fresh Julia sessions.
"""
function run_precompile_workload(; profile::Symbol=:test)
    strict = get(ENV, "PRECOMPILE_STRICT", "0") in ("1", "true", "TRUE", "yes", "YES")
    run_precompile_profile(profile; strict=strict)
    return nothing
end
