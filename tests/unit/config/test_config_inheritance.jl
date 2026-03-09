using Test

const _CONFIG_LOADER_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "config", "ConfigLoader.jl"))
if !isdefined(Main, :ConfigLoader)
    Base.include(Main, _CONFIG_LOADER_PATH)
end

using Main.ConfigLoader: load_config

const PHYSICS_CONFIG_DIR = normpath(joinpath(@__DIR__, "..", "..", "..", "config", "physics"))
const NJL_MODEL_CONFIG_DIR = normpath(joinpath(@__DIR__, "..", "..", "..", "config", "models", "njl"))
const PNJL_MODEL_CONFIG_DIR = normpath(joinpath(@__DIR__, "..", "..", "..", "config", "models", "pnjl"))
const NJL2_MODEL_CONFIG_DIR = normpath(joinpath(@__DIR__, "..", "..", "..", "config", "models", "njl2"))

const DEFAULT_PHYSICS_CONFIG = Dict{String, Any}(
    "physical" => Dict(
        "hbarc" => 197.3269804,
        "alpha_em" => 0.0072973525664,
    ),
)

@testset "Config inheritance" begin
    physics_data = load_config(PHYSICS_CONFIG_DIR, DEFAULT_PHYSICS_CONFIG; profile="default")
    physical_cfg = get(physics_data.config, "physical", Dict{String, Any}())
    shared_model_cfg = get(physics_data.config, "model_shared", Dict{String, Any}())

    @test Float64(get(physical_cfg, "hbarc", 0.0)) == 197.3269804
    @test Float64(get(physical_cfg, "alpha_em", 0.0)) == 0.0072973525664

    inherited_model_cfg = isempty(shared_model_cfg) ? Dict{String, Any}[] : [Dict("model" => shared_model_cfg)]

    njl_data = load_config(NJL_MODEL_CONFIG_DIR, Dict{String, Any}(); profile="default", inherited_configs=inherited_model_cfg)
    njl_model_cfg = get(njl_data.config, "model", Dict{String, Any}())
    @test Int(get(njl_model_cfg, "N_color", 0)) == 3
    @test Float64(get(njl_model_cfg, "Lambda_MeV", 0.0)) == 602.3
    @test Float64(get(njl_model_cfg, "m_s0_MeV", 0.0)) == 140.7

    pnjl_data = load_config(PNJL_MODEL_CONFIG_DIR, Dict{String, Any}(); profile="default", inherited_configs=inherited_model_cfg)
    pnjl_model_cfg = get(pnjl_data.config, "model", Dict{String, Any}())
    pnjl_polyakov_cfg = get(pnjl_data.config, "polyakov", Dict{String, Any}())
    @test Int(get(pnjl_model_cfg, "N_color", 0)) == 3
    @test Float64(get(pnjl_polyakov_cfg, "T0_MeV", 0.0)) == 210.0

    njl2_data = load_config(NJL2_MODEL_CONFIG_DIR, Dict{String, Any}(); profile="default")
    njl2_model_cfg = get(njl2_data.config, "model", Dict{String, Any}())
    @test Int(get(njl2_model_cfg, "N_color", 0)) == 3
    @test Int(get(njl2_model_cfg, "N_flavor", 0)) == 2
end
