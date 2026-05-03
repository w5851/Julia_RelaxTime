using Test

include(joinpath(@__DIR__, "..", "..", "..", "scripts", "dev", "gen_deps.jl"))

@testset "gen_deps parser filters doc examples" begin
    sample = """
    \"\"\"
    Example:
        include("FakeFromDoc.jl")
    \"\"\"
    # include("FakeFromComment.jl")
    include("RealDependency.jl") # trailing comment include("IgnoredTrailingComment.jl")
    using .LocalModule
    """

    sanitized = sanitize_source_text(sample)

    @test !occursin("FakeFromDoc.jl", sanitized)
    @test !occursin("FakeFromComment.jl", sanitized)
    @test occursin("include(\"RealDependency.jl\")", sanitized)
    @test occursin("using .LocalModule", sanitized)
end
