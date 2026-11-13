using Test

@testset "RelaxTime load order preserves module identity" begin
    root = normpath(joinpath(@__DIR__, "..", "..", ".."))
    # Fresh processes are required: an existing Main.RelaxTime would hide the
    # duplicate module/type/cache problem this test exercises.
    probe = raw"""
    using Test
    order, root = ARGS
    source(path) = joinpath(root, "src", path)
    legacy_paths = (
        (:ParameterTypes, "ParameterTypes.jl"),
        (:Constants_PNJL, "Constants_PNJL.jl"),
        (:PNJLQuarkDistributions, "QuarkDistribution.jl"),
        (:PNJLQuarkDistributions_Aniso, "QuarkDistribution_Aniso.jl"),
        (:ParameterAdapters, "utils/ParameterAdapters.jl"),
        (:ValidationUtils, "utils/ValidationUtils.jl"),
        (:GaussLegendre, "integration/GaussLegendre.jl"),
        (:OneLoopIntegrals, "relaxtime/OneLoopIntegrals.jl"),
        (:OneLoopIntegralsCorrection, "relaxtime/OneLoopIntegralsAniso.jl"),
        (:PolarizationAniso, "relaxtime/PolarizationAniso.jl"),
        (:PolarizationCache, "relaxtime/PolarizationCache.jl"),
        (:MesonInteractionKernel, "relaxtime/MesonInteractionKernel.jl"),
        (:KinematicChecks, "relaxtime/KinematicChecks.jl"),
    )
    if order == "legacy-first"
        for (name, path) in legacy_paths
            isdefined(Main, name) || Base.include(Main, source(path))
        end
    elseif order == "models-first"
        Base.include(Main, source("models/Models.jl"))
    end
    existing = Dict(name => getfield(Main, name) for (name, _) in legacy_paths if isdefined(Main, name))
    isdefined(Main, :RelaxTime) || Base.include(Main, source("relaxtime/RelaxTime.jl"))
    rt = Main.RelaxTime

    @testset "$order" begin
        for (name, mod) in existing
            @test getfield(Main, name) === mod
        end
        for name in (:ParameterAdapters, :ValidationUtils, :GaussLegendre,
                     :ParticleSymbols, :PhaseSpaceSampling, :KinematicChecks,
                     :OneLoopIntegrals, :OneLoopIntegralsCorrection,
                     :PolarizationAniso, :PolarizationCache, :MesonInteractionKernel,
                     :TotalPropagator, :ScatteringAmplitude, :TotalCrossSection,
                     :AverageScatteringRate, :RelaxationTime, :TransportCoefficients)
            @test getfield(Main, name) === getfield(rt, name)
        end
        @test Main.MesonInteractionKernel.FullKMTInteraction === rt.MesonInteractionKernel.FullKMTInteraction
        @test Main.AverageScatteringRate.ThermalWeightedSampling.WeightedSigmaTable ===
            rt.AverageScatteringRate.ThermalWeightedSampling.WeightedSigmaTable
        @test Main.KinematicChecks.check_kinematic_threshold === rt.DifferentialCrossSection.check_kinematic_threshold
        @test Main.GaussLegendre.standard_nodes_weights(7)[1] === rt.GaussLegendre.standard_nodes_weights(7)[1]
        for (name, path) in legacy_paths
            isdefined(Main, name) || Base.include(Main, source(path))
        end

        args = (:P, 1.0, 0.0, 1.5, 1.5, 0.1, 0.1, 0.7, 0.5, 0.5, 0.0, -2.0, -2.0, 0)
        Main.PolarizationCache.reset_cache!()
        try
            a = Main.PolarizationCache.polarization_aniso_cached(args...)
            b = rt.PolarizationCache.polarization_aniso_cached(args...)
            @test all(isfinite, a)
            @test a == b
            @test rt.PolarizationCache.get_cache_stats().cache_hits == 1
            @test rt.PolarizationCache.get_cache_stats().cache_size == 1
            rt.TotalPropagator.reset_cache!()
            @test Main.PolarizationCache.get_cache_stats().cache_size == 0
        finally
            Main.PolarizationCache.reset_cache!()
            rt.PolarizationCache.reset_cache!()
        end

        isdefined(Main, :Models) || Base.include(Main, source("models/Models.jl"))
        workflow = Main.Models.transport_workflow_module()
        @test workflow.RelaxationTime === rt.RelaxationTime
        @test workflow.AFieldBuilder === rt.AFieldBuilder
    end
    """
    for order in ("native-first", "legacy-first", "models-first")
        cmd = `$(Base.julia_cmd()) --project=$root --startup-file=no -e $probe $order $root`
        @test success(pipeline(cmd; stdout=stdout, stderr=stderr))
    end
end
