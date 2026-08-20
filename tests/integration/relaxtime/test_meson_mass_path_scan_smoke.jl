using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "meson mass path scan smoke" begin
    mktempdir() do outdir
        freezeout_outfile = joinpath(outdir, "freezeout_meson_mass_path_scan_smoke.csv")
        isentropic_outfile = joinpath(outdir, "isentropic_meson_mass_path_scan_smoke.csv")

        freezeout = Models.run_freezeout_meson_mass_scan(
            sqrt_s_NN_values=[7.7];
            output_path=freezeout_outfile,
            overwrite=true,
            p_num=8,
            t_num=4,
            max_iter=20,
            mesons=(:pi, :K),
        )
        @test freezeout.workflow_entry == "Models.run_freezeout_meson_mass_scan"
        @test isfile(freezeout_outfile)
        freezeout_text = read(freezeout_outfile, String)
        @test occursin("path_family,path_profile,path_segment,path_point_index,path_order_key,path_label", freezeout_text)
        @test occursin("root_sign_flipped_K", freezeout_text)
        @test occursin("freezeout", freezeout_text)

        isentropic = Models.run_isentropic_meson_mass_scan(
            T_MeV_values=[160.0];
            sigma_target=30.0,
            output_path=isentropic_outfile,
            overwrite=true,
            p_num=8,
            t_num=4,
            max_iter=20,
            mesons=(:pi,),
        )
        @test isentropic.workflow_entry == "Models.run_isentropic_meson_mass_scan"
        @test isfile(isentropic_outfile)
        isentropic_text = read(isentropic_outfile, String)
        @test occursin("sigma_target", isentropic_text)
        @test occursin("isentropic", isentropic_text)
        @test !occursin("seed_guess is required", isentropic_text)
    end
end
