using TestItems

@testitem "Unit Smoke" tags=[:unit, :smoke] default_imports=false begin
    using Test
    ENV["UNIT_PROFILE"] = "smoke"
    include(joinpath(@__DIR__, "..", "tests", "unit", "runtests.jl"))
end

@testitem "Integration Smoke" tags=[:integration, :smoke] default_imports=false begin
    using Test
    ENV["INTEGRATION_PROFILE"] = "smoke"
    include(joinpath(@__DIR__, "..", "tests", "integration", "runtests.jl"))
end

@testitem "Regression Smoke" tags=[:regression, :smoke] default_imports=false begin
    using Test
    ENV["REGRESSION_PROFILE"] = "smoke"
    include(joinpath(@__DIR__, "..", "tests", "regression", "runtests.jl"))
end

@testitem "Validation" tags=[:validation] default_imports=false begin
    using Test
    include(joinpath(@__DIR__, "..", "tests", "validation", "runtests.jl"))
end