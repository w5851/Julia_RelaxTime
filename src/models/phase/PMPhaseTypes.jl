const PM_BRANCH_STATUSES = (:accepted, :nonconverged, :branch_jump, :invalid_thermo)
const PM_SEED_SOURCES = (:seed0, :previous_same_branch, :manual_override)
const PM_ENDPOINT_CAUSES = (:physical_loss_candidate, :nonconvergence, :branch_jump, :out_of_grid)
const PM_COMPARISON_STATUSES = (:both, :pm_only, :maxwell_only, :neither)

Base.@kwdef struct PMSeedPair
    hadron_seed0::Vector{Float64}
    quark_seed0::Vector{Float64}
    continuity_mode::Symbol = :branch_local
    fallback_mode::Symbol = :none
end
