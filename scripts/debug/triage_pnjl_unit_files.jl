using Test

unit_dir = joinpath(@__DIR__, "..", "..", "tests", "unit", "pnjl") |> normpath

# Keep this consistent with tests/unit/runtests.jl (DEFAULT_SKIP)
const DEFAULT_SKIP = Set([
    "test_core_integrals.jl",
    "test_implicit_jacobian.jl",
    "test_scans.jl",
])

files = sort(filter(f ->
    endswith(f, ".jl") &&
    startswith(basename(f), "test_") &&
    !(lowercase(basename(f)) in DEFAULT_SKIP),
    readdir(unit_dir; join=true)
))

bad = String[]

for f in files
    try
        ts = @testset "$(basename(f))" begin
            include(f)
        end

        c = Test.get_test_counts(ts)
        if c.fails > 0 || c.errors > 0 || c.broken > 0
            push!(bad, "$(basename(f))  fail=$(c.fails)  err=$(c.errors)  broken=$(c.broken)")
        end
    catch e
        push!(bad, "$(basename(f))  THREW $(typeof(e))")
        @error "Exception while including test file" file=f exception=(e, catch_backtrace())
    end
end

println("PNJL problematic files:")
for line in bad
    println("  ", line)
end
println("Total bad: $(length(bad)) / $(length(files))")
