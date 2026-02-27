if !isdefined(Main, :ModelsTransportWorkflowBridge)
    @eval module ModelsTransportWorkflowBridge

    const _TRANSPORT_WORKFLOW_PATH = normpath(joinpath(@__DIR__, "TransportWorkflow.jl"))

    function module_ref()
        if !isdefined(Main, :TransportWorkflow)
            Base.include(Main, _TRANSPORT_WORKFLOW_PATH)
        end
        return Main.TransportWorkflow
    end

    end
end
