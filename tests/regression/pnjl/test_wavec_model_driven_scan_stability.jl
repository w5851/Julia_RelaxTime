using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

function _read_data_line(path::AbstractString)
    lines = readlines(path)
    length(lines) == 2 || error("expected exactly one data row in $(path)")
    return lines[2]
end

@testset "Wave-C model-driven scan stability" begin
    @testset "models backend supports non-PNJL model_kind" begin
        mktempdir() do outdir
            stats = Main.Models.run_tmu_scan(
                T_values=[150.0],
                mu_values=[0.0],
                xi_values=[0.0],
                output_path=joinpath(outdir, "models_rpnjl.csv"),
                overwrite=true,
                resume=false,
                use_phase_aware=false,
                solver_backend=:models,
                model_kind=:RPNJL,
                p_num=8,
                t_num=4,
            )
            @test stats.total == 1
            @test stats.success == 1
        end
    end

    @testset "Tmu model-driven outputs stay deterministic" begin
        mktempdir() do outdir
            out_pnjl_a = joinpath(outdir, "tmu_pnjl_a.csv")
            out_pnjl_b = joinpath(outdir, "tmu_pnjl_b.csv")
            out_rpnjl = joinpath(outdir, "tmu_rpnjl.csv")

            stats_pnjl_a = Main.Models.run_tmu_scan(
                T_values=[150.0],
                mu_values=[0.0],
                xi_values=[0.0],
                output_path=out_pnjl_a,
                overwrite=true,
                resume=false,
                use_phase_aware=false,
                solver_backend=:models,
                model_kind=:PNJL,
                p_num=8,
                t_num=4,
            )
            stats_pnjl_b = Main.Models.run_tmu_scan(
                T_values=[150.0],
                mu_values=[0.0],
                xi_values=[0.0],
                output_path=out_pnjl_b,
                overwrite=true,
                resume=false,
                use_phase_aware=false,
                solver_backend=:models,
                model_kind=:PNJL,
                p_num=8,
                t_num=4,
            )
            stats_rpnjl = Main.Models.run_tmu_scan(
                T_values=[150.0],
                mu_values=[0.0],
                xi_values=[0.0],
                output_path=out_rpnjl,
                overwrite=true,
                resume=false,
                use_phase_aware=false,
                solver_backend=:models,
                model_kind=:RPNJL,
                p_num=8,
                t_num=4,
            )

            @test stats_pnjl_a.total == 1
            @test stats_pnjl_b.total == 1
            @test stats_rpnjl.total == 1

            line_pnjl_a = _read_data_line(out_pnjl_a)
            line_pnjl_b = _read_data_line(out_pnjl_b)
            line_rpnjl = _read_data_line(out_rpnjl)
            @test line_pnjl_a == line_pnjl_b
            @test occursin("governance.selection=", line_pnjl_a)
            @test occursin("governance.selection=", line_rpnjl)

            @test length(split(line_pnjl_a, ',')) == length(split(line_rpnjl, ','))
        end
    end
end
