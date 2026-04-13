using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SCRIPT_PATH = joinpath(REPO_ROOT, "scripts", "relaxtime", "run_mott_phase_plot_modes.jl")

@testset "mott phase plot modes smoke" begin
    @test isfile(SCRIPT_PATH)

    tmp = mktempdir()
    in_csv = joinpath(tmp, "derived.csv")
    out_dir = joinpath(tmp, "figs")

    open(in_csv, "w") do io
        write(io, "# format: scan_csv_v1\n")
        write(io, "T_MeV,muB_MeV,xi,M_pi,M_K,Gamma_pi,Gamma_K,M_u_plus_M_d,M_u_plus_M_s\n")
        for xi in (-0.3, -0.15, 0.0, 0.15, 0.3)
            for T in (150.0, 160.0)
                write(io, string(T, ',', 0.0, ',', xi, ',', 0.7, ',', 0.8, ',', 0.01, ',', 0.02, ',', 1.0, ',', 1.2, "\n"))
            end
        end
    end

    python_cmd = Sys.which("python") !== nothing ? "python" : (Sys.which("python3") !== nothing ? "python3" : nothing)
    has_matplotlib = python_cmd !== nothing && success(`$(python_cmd) -c "import matplotlib"`)
    if !has_matplotlib
        @test true
        return
    end

    run(`julia --project=. $SCRIPT_PATH --in $in_csv --out-dir $out_dir`)

    @test isfile(joinpath(out_dir, "mode_a", "mott_mode_a__xi-0p3.png"))
    @test isfile(joinpath(out_dir, "mode_a", "mott_mode_a__xi-0p15.png"))
    @test isfile(joinpath(out_dir, "mode_a", "mott_mode_a__xi0.png"))
    @test isfile(joinpath(out_dir, "mode_a", "mott_mode_a__xi0p15.png"))
    @test isfile(joinpath(out_dir, "mode_a", "mott_mode_a__xi0p3.png"))

    observables = ["M_pi", "M_K", "Gamma_pi", "Gamma_K", "M_u_plus_M_d", "M_u_plus_M_s"]
    for obs in observables
        @test isfile(joinpath(out_dir, "mode_b", "mott_mode_b__$(obs).png"))
    end
end
