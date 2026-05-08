using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const OUTDIR = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "meson_mass", "path_scan", "test_outputs")
const FREEZEOUT_OUTFILE = joinpath(OUTDIR, "freezeout_meson_mass_path_scan_smoke.csv")
const ISENTROPIC_OUTFILE = joinpath(OUTDIR, "isentropic_meson_mass_path_scan_smoke.csv")

@testset "meson mass path scan smoke" begin
    mkpath(OUTDIR)
    isfile(FREEZEOUT_OUTFILE) && rm(FREEZEOUT_OUTFILE)
    isfile(ISENTROPIC_OUTFILE) && rm(ISENTROPIC_OUTFILE)

    freezeout = Models.run_freezeout_meson_mass_scan(
        sqrt_s_NN_values=[7.7];
        output_path=FREEZEOUT_OUTFILE,
        overwrite=true,
        p_num=8,
        t_num=4,
        max_iter=20,
        mesons=(:pi, :K),
    )
    @test freezeout.workflow_entry == "Models.run_freezeout_meson_mass_scan"
    @test isfile(FREEZEOUT_OUTFILE)
    freezeout_text = read(FREEZEOUT_OUTFILE, String)
    @test occursin("path_family,path_profile,path_segment,path_point_index,path_order_key,path_label", freezeout_text)
    @test occursin("freezeout", freezeout_text)

    isentropic = Models.run_isentropic_meson_mass_scan(
        T_MeV_values=[160.0];
        sigma_target=30.0,
        output_path=ISENTROPIC_OUTFILE,
        overwrite=true,
        p_num=8,
        t_num=4,
        max_iter=20,
        mesons=(:pi,),
    )
    @test isentropic.workflow_entry == "Models.run_isentropic_meson_mass_scan"
    @test isfile(ISENTROPIC_OUTFILE)
    isentropic_text = read(ISENTROPIC_OUTFILE, String)
    @test occursin("sigma_target", isentropic_text)
    @test occursin("isentropic", isentropic_text)
end
