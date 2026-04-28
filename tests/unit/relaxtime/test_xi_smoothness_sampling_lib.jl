using Test

const _SAMPLING_LIB_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "scripts", "relaxtime", "xi_smoothness_sampling_lib.jl"))

function _write_anchor_csv(path::String, rows::Vector{String}; header::String)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, header)
        for row in rows
            println(io, row)
        end
    end
end

function _capture_argument_error_message(f::Function)
    try
        f()
        return nothing
    catch err
        err isa ArgumentError || rethrow(err)
        return sprint(showerror, err)
    end
end

@testset "xi smoothness sampling lib" begin
    @test isfile(_SAMPLING_LIB_PATH)

    if isfile(_SAMPLING_LIB_PATH)
        Base.include(Main, _SAMPLING_LIB_PATH)
        using Main.XiSmoothnessSampling: sample_params

        tmp = mktempdir()
        boundary_csv = joinpath(tmp, "boundary.csv")
        crossover_csv = joinpath(tmp, "crossover.csv")

        _write_anchor_csv(
            boundary_csv,
            [
                "0.0,80.0,330.0,0.0,0.0",
                "0.0,120.0,300.0,0.0,0.0",
                "0.2,120.0,320.0,0.0,0.0",
                "0.0,200.0,120.0,0.0,0.0",
            ];
            header="xi,T_MeV,mu_transition_MeV,rho_hadron,rho_quark",
        )
        _write_anchor_csv(
            crossover_csv,
            [
                "0.0,0.0,210.0,180.0,0.0,0.0",
                "0.0,120.0,180.0,165.0,0.0,0.0",
                "0.4,120.0,170.0,155.0,0.0,0.0",
                "0.0,260.0,130.0,130.0,0.0,0.0",
            ];
            header="xi,mu_MeV,T_crossover_chiral_MeV,T_crossover_deconf_MeV,rho_chiral,rho_deconf",
        )

        cfg = (
            total = 24,
            seed = 20260427,
            random_count = 12,
            near_count = 12,
            T_range = (50.0, 270.0),
            muq_range = (0.0, 360.0),
            boundary_csv = boundary_csv,
            crossover_csv = crossover_csv,
        )

        rows = sample_params(
            cfg.total;
            seed=cfg.seed,
            random_count=cfg.random_count,
            near_count=cfg.near_count,
            T_range=cfg.T_range,
            muq_range=cfg.muq_range,
            boundary_csv=cfg.boundary_csv,
            crossover_csv=cfg.crossover_csv,
        )

        @test length(rows) == 24
        @test count(r -> r.source == "random_uniform", rows) == 12
        @test count(r -> r.source == "near_phase_line", rows) == 12
        @test all(r -> 50.0 <= r.T_MeV <= 270.0, rows)
        @test all(r -> 0.0 <= r.muq_MeV <= 360.0, rows)
        @test all(r -> isapprox(r.muB_MeV, 3.0 * r.muq_MeV; atol=1e-12, rtol=0.0), rows)

        expected_ids = ["S" * lpad(string(i), 3, '0') for i in 1:24]
        @test [r.sample_id for r in rows] == expected_ids

        rows_again = sample_params(
            cfg.total;
            seed=cfg.seed,
            random_count=cfg.random_count,
            near_count=cfg.near_count,
            T_range=cfg.T_range,
            muq_range=cfg.muq_range,
            boundary_csv=cfg.boundary_csv,
            crossover_csv=cfg.crossover_csv,
        )
        @test length(rows) == length(rows_again)
        for (r1, r2) in zip(rows, rows_again)
            for fn in fieldnames(typeof(r1))
                @test isequal(getfield(r1, fn), getfield(r2, fn))
            end
        end

        @testset "near phase anchors keep mu as mu_q (no extra division)" begin
            near_rows = sample_params(
                6;
                seed=cfg.seed,
                random_count=0,
                near_count=6,
                T_range=cfg.T_range,
                muq_range=cfg.muq_range,
                boundary_csv=cfg.boundary_csv,
                crossover_csv=cfg.crossover_csv,
            )
            observed_anchor_muq = sort(unique(r.anchor_muq_MeV for r in near_rows if r.source == "near_phase_line"))
            @test observed_anchor_muq == [0.0, 120.0, 260.0, 300.0, 330.0]
        end

        @testset "bad anchor row raises ArgumentError with file and line" begin
            bad_boundary_csv = joinpath(tmp, "bad_boundary.csv")
            _write_anchor_csv(
                bad_boundary_csv,
                [
                    "# comment before data row",
                    "0.0,80.0,330.0,0.0,0.0",
                    "0.0,not_a_number,300.0,0.0,0.0",
                ];
                header="xi,T_MeV,mu_transition_MeV,rho_hadron,rho_quark",
            )

            msg = _capture_argument_error_message() do
                sample_params(
                    1;
                    seed=cfg.seed,
                    random_count=0,
                    near_count=1,
                    T_range=cfg.T_range,
                    muq_range=cfg.muq_range,
                    boundary_csv=bad_boundary_csv,
                    crossover_csv=cfg.crossover_csv,
                )
            end
            @test msg !== nothing
            @test occursin(basename(bad_boundary_csv), msg)
            @test occursin(":4:", msg)
            @test occursin("T_MeV", msg)
        end

        @testset "missing required column raises clear error" begin
            bad_crossover_csv = joinpath(tmp, "bad_crossover_missing_col.csv")
            _write_anchor_csv(
                bad_crossover_csv,
                ["0.0,0.0,180.0,0.0,0.0"]; 
                header="xi,mu_MeV,T_crossover_deconf_MeV,rho_chiral,rho_deconf",
            )

            msg = _capture_argument_error_message() do
                sample_params(
                    1;
                    seed=cfg.seed,
                    random_count=0,
                    near_count=1,
                    T_range=cfg.T_range,
                    muq_range=cfg.muq_range,
                    boundary_csv=cfg.boundary_csv,
                    crossover_csv=bad_crossover_csv,
                )
            end
            @test msg !== nothing
            @test occursin("missing required column", msg)
            @test occursin("T_crossover_chiral_MeV", msg)
        end
    end
end
