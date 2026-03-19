using Test

const _OFFLINE_PATCH_UTILS_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "scripts", "relaxtime", "offline_patch_utils.jl"))
Base.include(Main, _OFFLINE_PATCH_UTILS_PATH)

using Main.OfflinePatchUtils: read_flagged_points, select_patch_points, PatchPoint

@testset "offline patch utils" begin
    @testset "read flagged points and deduplicate" begin
        tmp = tempname() * ".csv"
        open(tmp, "w") do io
            println(io, "# schema: scan_csv_v1")
            println(io, "T_MeV,muB_MeV,xi,quality_flag,quality_reason,quality_metric")
            println(io, "190,0,-0.1,true,tau_u_ubar_ratio_high,7.5")
            println(io, "190,0,-0.1,true,tau_u_ubar_ratio_high,8.0")
            println(io, "200,0,-0.2,false,ok,1.0")
            println(io, "200,0,-0.2,true,eta_over_s_nonfinite,NaN")
        end

        pts = read_flagged_points(tmp)
        @test length(pts) == 2
        @test pts[1].T_MeV == 190.0
        @test pts[1].quality_metric == 8.0
        @test pts[2].T_MeV == 200.0
        @test isnan(pts[2].quality_metric)
    end

    @testset "select by reason and top-k" begin
        pts = PatchPoint[
            PatchPoint(190.0, 0.0, -0.10, "tau_u_ubar_ratio_high", 8.0),
            PatchPoint(190.0, 0.0, -0.30, "tau_u_ubar_ratio_high", 6.2),
            PatchPoint(200.0, 0.0, -0.20, "eta_over_s_nonfinite", NaN),
        ]

        sel1 = select_patch_points(pts; reason_filter="tau_u_ubar_ratio_high")
        @test length(sel1) == 2

        sel2 = select_patch_points(pts; max_points=1)
        @test length(sel2) == 1
        @test sel2[1].quality_metric == 8.0
    end
end
