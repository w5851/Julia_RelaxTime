using Test

if !isdefined(Main, :TransportBaselineExport)
    include(joinpath(@__DIR__, "..", "..", "..", "scripts", "dev", "export_transport_fixedpoint_baseline.jl"))
end
const TBE = Main.TransportBaselineExport

_fixture_transport(pt) = (equilibrium=(converged=true,), transport=(eta=pt.T, sigma=pt.mu, zeta=pt.xi))

@testset "transport export requires a new explicit destination" begin
    loaded_before = isdefined(Main, :Models)
    @test TBE.main(["--help"]; io=devnull) == 0
    @test TBE.main(String[]; io=devnull) == 1
    @test TBE.main(["--output"]; io=devnull) == 1
    mktempdir() do root
        output = joinpath(root, "candidate.csv")
        write(output, "existing baseline")
        @test TBE.main(["--output", output]; io=devnull) == 1
        @test read(output, String) == "existing baseline"
        @test isdefined(Main, :Models) == loaded_before
    end
end

@testset "transport candidate publication and failure cleanup" begin
    mktempdir() do root
        output = joinpath(root, "candidate.csv")
        @test TBE.main(["--output", output]; io=devnull, compute_point=_fixture_transport) == 0
        lines = readlines(output)
        @test length(lines) == 14
        @test first(lines) == "T,mu,xi,eta,sigma,zeta"
        first_values = parse.(Float64, split(lines[2], ','))
        @test first_values == [0.5, 0.0, 0.0, 0.5, 0.0, 0.0]
        @test sort(readdir(root)) == ["candidate.csv"]
    end
    for failure in (:exception, :nonfinite, :unconverged)
        mktempdir() do root
            output = joinpath(root, "candidate.csv")
            calls = Ref(0)
            compute = pt -> begin
                calls[] += 1
                calls[] < 3 && return _fixture_transport(pt)
                failure == :exception && error("synthetic solver failure")
                failure == :nonfinite && return (equilibrium=(converged=true,), transport=(eta=NaN, sigma=1.0, zeta=1.0))
                return (equilibrium=(converged=false,), transport=(eta=1.0, sigma=1.0, zeta=1.0))
            end
            report = IOBuffer()
            @test TBE.main(["--output", output]; io=report, compute_point=compute) == 1
            @test occursin("T = 0.66, mu = 0.5, xi = 0.0", String(take!(report)))
            @test calls[] == 3
            @test !ispath(output)
            @test isempty(readdir(root))
        end
    end
    mktempdir() do root
        output = joinpath(root, "candidate.csv")
        compute = pt -> begin
            pt == last(TBE.POINTS) && write(output, "concurrent writer")
            return _fixture_transport(pt)
        end
        @test TBE.main(["--output", output]; io=devnull, compute_point=compute) == 1
        @test read(output, String) == "concurrent writer"
        @test readdir(root) == ["candidate.csv"]
    end
end
