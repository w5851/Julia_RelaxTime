"""
Test runner for parameter struct migration tests.

Run all struct migration tests:
    julia --project=. -e 'include("tests/unit/struct_migration/runtests.jl")'

Run with verbose output:
    julia --project=. -e 'ENV["STRUCT_MIGRATION_VERBOSE"]="1"; include("tests/unit/struct_migration/runtests.jl")'
"""

using Test

const STRUCT_MIGRATION_DIR = @__DIR__

# Test files to run (in order)
const TEST_FILES = [
    "test_parameter_types.jl",
    "test_particle_symbols_struct.jl",
    "test_particle_symbols_edge_cases.jl",
    "test_differential_cross_section_property.jl",
    "test_differential_cross_section_edge_cases.jl",
    "test_scattering_amplitude_property.jl",
    "test_scattering_amplitude_edge_cases.jl",
    "test_total_propagator_property.jl",
    "test_total_propagator_edge_cases.jl",
    "test_total_cross_section_property.jl",
    "test_total_cross_section_edge_cases.jl",
    "test_average_scattering_rate_property.jl",
    "test_average_scattering_rate_edge_cases.jl",
    "test_relaxation_time_property.jl",
    "test_relaxation_time_edge_cases.jl",
    # Integration and validation tests
    "test_integration_full_chain.jl",
    "test_conversion_roundtrip_property.jl",
    "test_ensure_quark_params_has_A_property.jl",
]

@testset "Parameter Struct Migration" begin
    for test_file in TEST_FILES
        test_path = joinpath(STRUCT_MIGRATION_DIR, test_file)
        if isfile(test_path)
            @testset "$test_file" begin
                include(test_path)
            end
        else
            @warn "Test file not found: $test_path"
        end
    end
end

println("\n" * "="^70)
println("Parameter Struct Migration Tests Complete!")
println("="^70)
