const WORKFLOW_STAGE_SKELETON = (
    :prepare_inputs,
    :solve_core,
    :postprocess,
    :export_artifacts,
    :emit_diagnostics,
    :emit_repro_manifest,
)

@inline workflow_pipeline_stage_ids() = WORKFLOW_STAGE_SKELETON
