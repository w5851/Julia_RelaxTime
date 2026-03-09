using TestItems

@testitem "Unit Smoke" tags=[:unit, :smoke] default_imports=false begin
    using Test
    ENV["UNIT_PROFILE"] = "smoke"
    include(joinpath(@__DIR__, "unit.jl"))
end

@testitem "Integration Smoke" tags=[:integration, :smoke] default_imports=false begin
    using Test
    ENV["INTEGRATION_PROFILE"] = "smoke"
    include(joinpath(@__DIR__, "integration.jl"))
end

@testitem "Regression Smoke" tags=[:regression, :smoke] default_imports=false begin
    using Test
    ENV["REGRESSION_PROFILE"] = "smoke"
    include(joinpath(@__DIR__, "regression.jl"))
end

@testitem "Validation" tags=[:validation] default_imports=false begin
    using Test
    include(joinpath(@__DIR__, "validation.jl"))
end