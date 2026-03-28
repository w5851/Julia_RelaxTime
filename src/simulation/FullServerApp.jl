module FullServerApp

using HTTP
using JSON3
using LinearAlgebra
using UUIDs: uuid4

include(joinpath(@__DIR__, "MomentumMapping.jl"))
using .MomentumMapping

include(joinpath(@__DIR__, "..", "models", "Models.jl"))
using .Models

const _CONSTANTS_PATH = normpath(joinpath(@__DIR__, "..", "constants", "Constants_PNJL.jl"))
if !isdefined(Main, :Constants_PNJL)
    Base.include(Main, _CONSTANTS_PATH)
end
using Main.Constants_PNJL: ħc_MeV_fm

include(joinpath(@__DIR__, "fullserver", "shared.jl"))
include(joinpath(@__DIR__, "fullserver", "compute_handlers.jl"))
include(joinpath(@__DIR__, "fullserver", "pnjl_handlers.jl"))
include(joinpath(@__DIR__, "fullserver", "pnjl_scan_jobs.jl"))
include(joinpath(@__DIR__, "fullserver", "http_utils.jl"))
include(joinpath(@__DIR__, "fullserver", "routing.jl"))

end
