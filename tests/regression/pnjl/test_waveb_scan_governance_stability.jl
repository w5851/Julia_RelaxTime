using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

function _read_last_data_line(path::AbstractString)
    lines = readlines(path)
    length(lines) >= 2 || error("expected data rows in $(path)")
    return lines[end]
end

@testset "Wave-B scan governance stability" begin
    @testset "Tmu scan emits governance selection tag" begin
        mktempdir() do outdir
            out = joinpath(outdir, "tmu_waveb_governance.csv")
            stats = Main.Models.run_tmu_scan(
                T_values=[150.0],
                mu_values=[0.0],
                xi_values=[0.0],
                output_path=out,
                overwrite=true,
                resume=false,
                use_phase_aware=false,
                solver_backend=:models,
                p_num=12,
                t_num=4,
                iterations=60,
            )
            @test stats.total == 1
            line = _read_last_data_line(out)
            @test occursin("governance.selection=", line)
        end
    end

    @testset "Trho scan emits governance selection tag" begin
        mktempdir() do outdir
            out = joinpath(outdir, "trho_waveb_governance.csv")
            stats = Main.Models.run_trho_scan(
                T_values=[150.0],
                rho_values=[0.2],
                xi_values=[0.0],
                output_path=out,
                overwrite=true,
                resume=false,
                reverse_rho=false,
                seed_policy=:candidates,
                solver_backend=:models,
                p_num=12,
                t_num=4,
                iterations=60,
            )
            @test stats.total == 1
            line = _read_last_data_line(out)
            @test occursin("governance.selection=", line)
        end
    end

    @testset "fallback reason tags are normalized and ordered" begin
        msg = Main.Models.ScanCommon.join_messages([
            "seed[quark] failed (iterations=20, residual=0.100000, converged=false): bad",
            "hybrid weighted-block fallback rescued",
            "governance.selection=pressure_max_under_constraints;seed=hadron",
        ])

        @test msg == "governance.selection=pressure_max_under_constraints;seed=hadron | fallback.reason=seed_failed | fallback.reason=weighted_block_rescued"
    end
end
