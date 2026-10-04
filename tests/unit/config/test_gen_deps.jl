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

@testset "gen_deps resolves paths without executing source" begin
    mktempdir() do root
        src = joinpath(root, "src")
        mkpath(joinpath(src, "models"))
        for file in ("literal.jl", "variable.jl", "nested.jl", "base.jl")
            write(joinpath(src, file), "")
        end
        write(joinpath(src, "models", "entry.jl"), raw"""
        error("source must never execute")
        """ * "\n" * raw"""
        const ROOT = normpath(joinpath(@__DIR__, ".."))
        include("../literal.jl")
        path = joinpath(ROOT, "variable.jl")
        Base.include(Main, path)
        if enabled
            nested_path = joinpath(ROOT, "nested.jl")
            include(nested_path)
        end
        include(Base.joinpath(dirname(dirname(@__FILE__)), "base.jl"))
        # include("comment.jl")
        example = "include('string.jl')"
        quote
            include("quoted.jl")
        end
        """)
        unresolved = String[]
        adj = collect_dependency_graph(root; unresolved)
        @test adj["src/models/entry.jl"] == Set("src/" .* ["literal.jl", "variable.jl", "nested.jl", "base.jl"])
        @test isempty(unresolved)
    end
end

@testset "gen_deps reports dynamic paths and respects scope" begin
    mktempdir() do root
        mkpath(joinpath(root, "src"))
        write(joinpath(root, "src", "known.jl"), "")
        write(joinpath(root, "src", "entry.jl"), raw"""
        path = joinpath(@__DIR__, "known.jl")
        function load_file(path)
            include(path)
        end
        module Inner
            include(path)
        end
        include(path)
        if enabled
            path = "other.jl"
        end
        include(path)
        include(abspath("depends_on_cwd.jl"))
        include("missing.jl")
        """)
        unresolved = String[]
        adj = collect_dependency_graph(root; unresolved)
        @test adj["src/entry.jl"] == Set(["src/known.jl", "src/missing.jl"])
        @test count(x -> occursin("include(path)", x), unresolved) == 3
        @test any(x -> occursin("depends_on_cwd", x), unresolved)
        @test any(x -> occursin("missing include target src/missing.jl", x), unresolved)
        write(joinpath(root, "src", "broken.jl"), "function broken(")
        empty!(unresolved)
        collect_dependency_graph(root; unresolved)
        @test any(x -> occursin("parse error", x), unresolved)
    end
end
