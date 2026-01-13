using Test

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const UNIT_DIR = joinpath(ROOT, "tests", "unit")
const INTEG_DIR = joinpath(UNIT_DIR, "integration")

function should_include(path::String)
    file = lowercase(basename(path))
    return endswith(file, ".jl") && startswith(file, "test_")
end

files = sort(filter(should_include, readdir(INTEG_DIR; join=true)))

@testset "Integration triage" begin
    for f in files
        @testset "$(basename(f))" begin
            try
                include(f)
            catch err
                # Make the error visible, but keep going.
                @test false
                @info "Errored file" file=f exception=(err, catch_backtrace())
            end
        end
    end
end
