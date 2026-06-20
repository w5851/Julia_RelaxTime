const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !(REPO_ROOT in LOAD_PATH)
    insert!(LOAD_PATH, min(2, length(LOAD_PATH) + 1), REPO_ROOT)
end
if !isdefined(Main, :Models)
    Base.include(Main, joinpath(REPO_ROOT, "src", "models", "Models.jl"))
end

include(joinpath(@__DIR__, "fixedmu_palc_common.jl"))
include(joinpath(@__DIR__, "fixedasymrho_palc_common.jl"))

using .FixedAsymRhoPALCSpike

FixedAsymRhoPALCSpike.main_run(ARGS; repo_root=REPO_ROOT, default_run_reference=false)
