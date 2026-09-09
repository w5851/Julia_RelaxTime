using Test,TOML
isdefined(Main,:Models) || Base.include(Main,joinpath(@__DIR__,"..","..","..","src","models","Models.jl"))
include(joinpath(@__DIR__,"..","..","..","src","models","workflow_apps","ChargedGBUResearchWorkflow.jl"))
const W=ChargedGBUResearchWorkflow
@testset "Infinite GBU explicit production contract" begin
    c=TOML.parsefile(W.DEFAULT_CONFIG)
    @test W.validate_config(c)===c
    @test endswith(replace(W.default_figure_output("data/outputs/results/relaxtime/meson_density/charged_gbu_infinite/freezeout_test"),'\\'=>'/'),
        "data/outputs/figures/relaxtime/meson_density/charged_gbu_infinite/freezeout_test")
    @test W.energy_grid([3.,200.,7.7])==[200.,7.7,3.]
    for e in (Float64[],[0.],[NaN],[3.,3.]);@test_throws ArgumentError W.energy_grid(e);end
    for (key,value) in (("thermal_target","finite"),("meson_feedback",true),("production_default",false),("default_tier","formal"))
        bad=deepcopy(c);bad[key]=value;@test_throws ArgumentError W.validate_config(bad)
    end
    bad=deepcopy(c);bad["gates"]["q_relative"]=.1;@test_throws ArgumentError W.validate_config(bad)
    pt=(T_MeV=150.,muB_MeV=300.)
    r=Dict(ch=>(passed=true,density=ch[1]=='K' ? .2 : 1.) for ch in ("pi_plus","pi_minus","K_plus","K_minus"))
    row=W.ratio_row(7.7,pt,r)
    @test row.plus_passed && row.minus_passed
    @test row.Kplus_over_pi_plus==.2
    r["pi_plus"]=(passed=false,density=1.)
    @test isnan(W.ratio_row(7.7,pt,r).Kplus_over_pi_plus)
    @test !W.ratio_row(7.7,pt,Dict()).plus_passed
    a=(density=1.,bound=2.,landau=.1,pair=-1.1,omega_tail=1e-12)
    tails=[(density=10.0^(-j-6),omega_tail=1e-14) for j in 1:4]
    @test W.reduce_bands(a,a,tails,c["gates"]).passed
    @test !W.reduce_bands(merge(a,(bound=2.1,pair=-1.2)),a,tails,c["gates"]).passed
    @test !W.reduce_bands(a,a,reverse(tails),c["gates"]).passed
    @test W.jsonsafe((density=NaN,))["density"]===nothing
    mktempdir() do dir
        W.writejson(joinpath(dir,"result.json"),(density=1.,))
        W.writejson(joinpath(dir,"manifest.json"),(version=1,))
        hashes=W.output_hashes(dir)
        @test haskey(hashes,"result.json") && haskey(hashes,"result.json.sha256")
        @test !haskey(hashes,"manifest.json")
        @test !haskey(hashes,"manifest.json.sha256")
        W.writejson(joinpath(dir,"manifest.json"),(version=2,))
        @test W.output_hashes(dir)==hashes
    end
end
