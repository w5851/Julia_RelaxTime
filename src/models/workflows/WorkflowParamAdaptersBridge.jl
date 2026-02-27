if !isdefined(Main, :ModelsWorkflowParamAdaptersBridge)
    @eval module ModelsWorkflowParamAdaptersBridge

    const _WORKFLOW_PARAM_ADAPTERS_PATH = normpath(joinpath(@__DIR__, "WorkflowParamAdapters.jl"))

    function module_ref()
        if !isdefined(Main, :WorkflowParamAdapters)
            Base.include(Main, _WORKFLOW_PARAM_ADAPTERS_PATH)
        end
        return Main.WorkflowParamAdapters
    end

    end
end
