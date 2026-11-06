using Test
using CSV
using JSON3

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const EVALUATOR = joinpath(ROOT, "scripts", "analysis", "pnjl_c2_cep_limited_feasibility.jl")
const JOB = joinpath(ROOT, "scripts", "analysis", "pnjl_c2_cep_limited_feasibility_job.jl")

include(JOB)
module CEPReplayEvaluatorContract
include(Main.EVALUATOR)
end

@testset "C2 CEP limited feasibility contracts" begin
    @testset "target xi subset remains explicit" begin
        @test CEPReplayEvaluatorContract._parse_target_xi("0.125,0.39375,0.5") ==
            (0.125, 0.39375, 0.5)
        @test CEPReplayEvaluatorContract._parse_target_xi("") == Main.TARGET_XI
        @test_throws ErrorException CEPReplayEvaluatorContract._parse_target_xi("0.126")
        @test_throws ErrorException CEPReplayEvaluatorContract._parse_target_xi("0.125,0.125")
    end

    @testset "temperature keys use Julia 1.12-compatible rounding" begin
        @test CEPReplayEvaluatorContract._temperature_key(113.1328125) == 113.1328125
        @test CEPReplayEvaluatorContract._temperature_key(113.1328125004) == 113.1328125
    end

    @testset "CEP fine-pool materialization is solver-free and complete" begin
        xi = -0.45
        temperature = 100.0
        sha = "4c9703c3be45b76608ab57d375082e29418bfd05"
        function write_fixture(root, rows)
            oracle = joinpath(root, "oracle")
            eval_dir = joinpath(oracle, "production_eval")
            mkpath(eval_dir)
            CSV.write(joinpath(eval_dir, "prod_eval_T100p000000_memoized.csv"), rows)
            CSV.write(joinpath(oracle, "trho_scan.csv"), rows[1:min(length(rows), 3)])
            oracle
        end
        rows = [(xi=xi, T_MeV=temperature, rho=i * Main.RHO_FINE_STEP,
            mu_avg_MeV=300.0 + i, residual_norm=0.0, iterations=2, converged=true)
            for i in 0:round(Int, Main.RHO_MAX / Main.RHO_FINE_STEP)]
        mktempdir() do tmp
            oracle = write_fixture(tmp, rows)
            materialized = Main._materialize_fine_pool(oracle, xi, temperature, sha)
            @test materialized.rows == length(Main._rho_grid())
            @test materialized.recovered_rows == length(Main._rho_grid())
            @test materialized.aggregate_rows == 3
            @test isfile(materialized.path)
            @test isfile(materialized.provenance_path)
            provenance = JSON3.read(read(materialized.provenance_path, String))
            @test provenance.solver_called == false
            @test provenance.method == "production_eval_materialization_v1"
            @test length(Main._curve_rows(materialized.path, xi, temperature, sha)) ==
                length(Main._rho_grid())
        end

        @testset "CSV coordinate round-off is normalized" begin
            rounded_rows = [(xi=xi, T_MeV=(i == 0 ? temperature - 5e-7 : temperature),
                rho=i * Main.RHO_FINE_STEP, mu_avg_MeV=300.0 + i,
                residual_norm=0.0, iterations=2, converged=true)
                for i in 0:round(Int, Main.RHO_MAX / Main.RHO_FINE_STEP)]
            mktempdir() do tmp
                oracle = write_fixture(tmp, rounded_rows)
                materialized = Main._materialize_fine_pool(oracle, xi, temperature, sha)
                first_row = first(collect(CSV.File(materialized.path)))
                @test Float64(first_row.T_MeV) == temperature
                @test length(Main._curve_rows(materialized.path, xi, temperature, sha)) ==
                    length(Main._rho_grid())
            end
        end

        @testset "missing and duplicate keys are rejected" begin
            mktempdir() do tmp
                @test_throws ErrorException Main._materialize_fine_pool(
                    write_fixture(tmp, rows[1:end-1]), xi, temperature, sha)
            end
            mktempdir() do tmp
                duplicate_rows = vcat(rows, [first(rows)])
                @test_throws ErrorException Main._materialize_fine_pool(
                    write_fixture(tmp, duplicate_rows), xi, temperature, sha)
            end
        end
    end
end
