# Legacy compatibility shim for removed path `src/models/pnjl/core/Integrals.jl`
# Forward to the canonical implementation under `src/models/pnjl_physics/core/Integrals.jl`.

include(normpath(joinpath(@__DIR__, "..", "..", "pnjl_physics", "core", "Integrals.jl")))
