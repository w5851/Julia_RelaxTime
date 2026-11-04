using Test

module MottDiagnosticEntrypointFixtures

const SCRIPT_DIR = normpath(joinpath(@__DIR__, "..", "..", "..", "scripts", "analysis", "relaxtime"))

function namespace()
    scope = Module(gensym(:MottEntrypoint))
    Core.eval(scope, :(include(path) = Base.include(@__MODULE__, path)))
    return scope
end

function dependency(path, marker)
    mkpath(dirname(path))
    # Stop at the actual include, before any physical setup or calculation.
    write(path, "const fixture_loaded = :$(marker)\nerror(\"diagnostic fixture stop\")\n")
end

function load_entrypoint(scope, path)
    try
        Base.include(scope, path)
        return nothing
    catch err
        return err
    end
end

end

@testset "Mott diagnostics resolve repository and frozen library names" begin
    fixtures = MottDiagnosticEntrypointFixtures
    entries = (
        ("audit_sigma_mott_cutoff.jl", "run_sigma_mott_spectrum.jl", "run_sigma_mott_spectrum_v2.jl"),
        ("run_sigma_mott_precision_reference.jl", "run_sigma_mott_spectrum.jl", "run_sigma_mott_spectrum_v2.jl"),
        ("run_sigma_mott_thermal_audit.jl", "run_sigma_mott_spectrum.jl", "run_sigma_mott_spectrum_v2.jl"),
        ("check_sigma_mott_asymptotics.jl", "run_sigma_mott_asymptotics.jl", "run_sigma_mott_asymptotics_v3.jl"),
        ("check_sigma_mott_background_asymptotics.jl", "run_sigma_mott_asymptotics.jl", "run_sigma_mott_asymptotics_v3.jl"),
    )
    for (entry, canonical, archived) in entries
        @testset "$entry" begin
            mktempdir() do folder
                script = joinpath(folder, entry)
                cp(joinpath(fixtures.SCRIPT_DIR, entry), script)
                scope = fixtures.namespace()
                @test fixtures.load_entrypoint(scope, script) isa LoadError
                @test !isdefined(scope, :fixture_loaded)

                fixtures.dependency(joinpath(folder, canonical), :canonical)
                scope = fixtures.namespace()
                err = fixtures.load_entrypoint(scope, script)
                @test err isa LoadError
                @test occursin("diagnostic fixture stop", sprint(showerror, err))
                @test Base.invokelatest(getfield, scope, :fixture_loaded) === :canonical

                fixtures.dependency(joinpath(folder, archived), :archived)
                scope = fixtures.namespace()
                err = fixtures.load_entrypoint(scope, script)
                @test err isa LoadError
                @test occursin("diagnostic fixture stop", sprint(showerror, err))
                @test Base.invokelatest(getfield, scope, :fixture_loaded) === :archived
            end
        end
    end
end

@testset "Mott midpoint screening honors the explicit library override" begin
    fixtures = MottDiagnosticEntrypointFixtures
    mktempdir() do folder
        script = joinpath(folder, "screen_transport_other_mott_midpoints.jl")
        cp(joinpath(fixtures.SCRIPT_DIR, basename(script)), script)
        fixtures.dependency(joinpath(folder, "run_transport_other_mott_refinement.jl"), :canonical)
        custom = joinpath(folder, "selected", "refinement.jl")
        fixtures.dependency(custom, :custom)
        for (value, expected) in ((nothing, :canonical), (joinpath("selected", "refinement.jl"), :custom), (custom, :custom))
            withenv("MOTT_REFINEMENT_LIBRARY" => value) do
                scope = fixtures.namespace()
                err = fixtures.load_entrypoint(scope, script)
                @test err isa LoadError
                @test occursin("diagnostic fixture stop", sprint(showerror, err))
                @test Base.invokelatest(getfield, scope, :fixture_loaded) === expected
            end
        end
    end
end
