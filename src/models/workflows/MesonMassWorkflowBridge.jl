if !isdefined(Main, :ModelsMesonMassWorkflowBridge)
    @eval module ModelsMesonMassWorkflowBridge

    const _MESON_WORKFLOW_PATH = normpath(joinpath(@__DIR__, "MesonMassWorkflow.jl"))

    function module_ref()
        if !isdefined(Main, :MesonMassWorkflow)
            Base.include(Main, _MESON_WORKFLOW_PATH)
        end
        return Main.MesonMassWorkflow
    end

    end
end
