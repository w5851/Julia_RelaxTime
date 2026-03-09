const VALIDATION_PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const VALIDATION_DATA_ROOT = joinpath(VALIDATION_PROJECT_ROOT, "tests", "validation", "data")

validation_data_path(parts...) = joinpath(VALIDATION_DATA_ROOT, parts...)
validation_raw_long_path(parts...) = validation_data_path("raw_long", parts...)
validation_targets_path(parts...) = validation_data_path("targets", parts...)
validation_provenance_path(parts...) = validation_data_path("provenance", parts...)