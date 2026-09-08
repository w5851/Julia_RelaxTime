using Test
const ROOT=normpath(joinpath(@__DIR__,"..","..",".."))
isdefined(Main,:Models) || Base.include(Main,joinpath(ROOT,"src","models","Models.jl"))
include(joinpath(ROOT,"scripts","relaxtime","run_charged_gbu_freezeout_scan.jl"))
@testset "Opt-in Models/CLI wiring leaves legacy path unchanged" begin
    @test isdefined(Main.Models,:run_charged_gbu_freezeout_scan)
    @test isdefined(Main.Models,:run_freezeout_meson_density_scan)
    cli=ChargedGBUFreezeoutCLI
    @test cli.parse_args(["--help"])===nothing
    opts=cli.parse_args(["--sqrts-list","3,7.7,200","--workers","1","--no-plot","--resume"])
    @test opts[:energies]==[3.,7.7,200.]
    @test opts[:workers]==1 && opts[:resume] && !opts[:make_plot]
    @test_throws ErrorException cli.parse_args(["--unknown"])
    @test_throws ErrorException cli.parse_args(["--output"])
    # Dispatch and validation are exercised without starting an equilibrium solve.
    @test_throws ArgumentError Main.Models.run_charged_gbu_freezeout_scan(output="unused",energies=[NaN])
    w=Main.ChargedGBUResearchWorkflow
    mktempdir() do dir
        path=joinpath(dir,"checkpoint.json")
        w.writejson(path,(density=NaN,passed=false))
        @test w.readjson(path).density===nothing
        open(io->write(io,"{}"),path,"w")
        @test_throws ErrorException w.readjson(path)
    end
    legacy=read(joinpath(ROOT,"scripts","relaxtime","run_freezeout_meson_density_scan.jl"),String)
    @test occursin(":regime => :stable",legacy)
end
