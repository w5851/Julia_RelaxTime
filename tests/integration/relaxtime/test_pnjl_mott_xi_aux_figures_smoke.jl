using Test
using Base64

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SCRIPT_PATH = joinpath(REPO_ROOT, "scripts", "analysis", "gen_pnjl_mott_xi_aux_figures.py")

function _write_png(path::String)
    # 1x1 transparent PNG
    b64 = "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO7Z5JkAAAAASUVORK5CYII="
    open(path, "w") do io
        write(io, base64decode(b64))
    end
end

@testset "pnjl mott xi aux figures smoke" begin
    @test isfile(SCRIPT_PATH)

    python_cmd = Sys.which("python") !== nothing ? "python" : (Sys.which("python3") !== nothing ? "python3" : nothing)
    has_matplotlib = python_cmd !== nothing && success(`$(python_cmd) -c "import matplotlib"`)
    if !has_matplotlib
        @test true
        return
    end

    tmp = mktempdir()
    derived_csv = joinpath(tmp, "derived.csv")
    scan_csv = joinpath(tmp, "scan.csv")
    gap_csv = joinpath(tmp, "gap.csv")
    mode_ab_dir = joinpath(tmp, "mode_ab")
    fig_dir = joinpath(tmp, "figures")
    doc_path = joinpath(tmp, "xi_dependence_analysis.md")
    mkpath(mode_ab_dir)

    _write_png(joinpath(mode_ab_dir, "mott_mode_ab__M_K__xi3.png"))
    _write_png(joinpath(mode_ab_dir, "mott_mode_ab__M_pi__xi3.png"))

    open(derived_csv, "w") do io
        write(io, "T_MeV,muB_MeV,xi,M_pi,M_K,Gamma_pi,Gamma_K,M_u_plus_M_d,M_u_plus_M_s\n")
        for xi in (-0.3, 0.0, 0.3)
            T0_pi = 198.0 + 20.0 * xi
            T0_k = 196.0 + 20.0 * xi
            for T in (190.0, 195.0, 200.0, 205.0, 210.0, 215.0, 220.0)
                mthr_pi = 1.0 + 0.002 * (T - 190.0) + 0.1 * xi
                mthr_k = 1.8 + 0.003 * (T - 190.0) + 0.08 * xi
                mpi = mthr_pi + 0.02 * (T - T0_pi)
                mk = mthr_k + 0.018 * (T - T0_k)
                gpi = max(0.0, 0.02 * (T - T0_pi))
                gk = max(0.0, 0.018 * (T - T0_k))
                write(io, string(T, ',', 0.0, ',', xi, ',', mpi, ',', mk, ',', gpi, ',', gk, ',', mthr_pi, ',', mthr_k, "\n"))
            end
        end
    end

    open(scan_csv, "w") do io
        write(io, "T_MeV,muB_MeV,xi,M_pi,M_K,Gamma_pi,Gamma_K,m_u,m_s\n")
        for xi in (-0.3, 0.0, 0.3)
            T0_pi = 198.0 + 20.0 * xi
            T0_k = 196.0 + 20.0 * xi
            for T in (190.0, 195.0, 200.0, 205.0, 210.0, 215.0, 220.0)
                mu = 0.3 + 0.004 * (T - 190.0) + 0.2 * xi
                ms = 1.9 + 0.002 * (T - 190.0) + 0.1 * xi
                mthr_pi = 2.0 * mu
                mthr_k = mu + ms
                mpi = mthr_pi + 0.02 * (T - T0_pi)
                mk = mthr_k + 0.018 * (T - T0_k)
                gpi = max(0.0, 0.02 * (T - T0_pi))
                gk = max(0.0, 0.018 * (T - T0_k))
                write(io, string(T, ',', 0.0, ',', xi, ',', mpi, ',', mk, ',', gpi, ',', gk, ',', mu, ',', ms, "\n"))
            end
        end
    end

    open(gap_csv, "w") do io
        write(io, "T_MeV,muB_MeV,xi,Phi,Phibar,m_u,m_s\n")
        for xi in (-0.3, 0.0, 0.3)
            for T in (190.0, 195.0, 200.0, 205.0, 210.0, 215.0, 220.0)
                phi = 0.62 + 0.004 * (T - 200.0) - 0.05 * xi
                phibar = phi
                mu = 0.35 + 0.003 * (T - 200.0) + 0.1 * xi
                ms = 2.0 + 0.001 * (T - 200.0) + 0.05 * xi
                write(io, string(T, ',', 0.0, ',', xi, ',', phi, ',', phibar, ',', mu, ',', ms, "\n"))
            end
        end
    end

    open(doc_path, "w") do io
        write(io, "# PNJL/Mott xi\n\n## 5. 图像读取要点\n")
    end

    run(`$(python_cmd) $SCRIPT_PATH --derived-csv $derived_csv --scan-csv $scan_csv --gap-csv $gap_csv --mode-ab-dir $mode_ab_dir --fig-dir $fig_dir --doc $doc_path`)

    @test isfile(joinpath(fig_dir, "fig1_tmott_vs_xi_fit.png"))
    @test isfile(joinpath(fig_dir, "fig4_orderparam_direct_indirect.png"))
    @test isfile(joinpath(fig_dir, "fig2_gamma_delta_dualaxis.png"))

    @test isfile(joinpath(mode_ab_dir, "mott_mode_ab__M_K__xi3_annotated.png"))
    @test isfile(joinpath(mode_ab_dir, "mott_mode_ab__M_pi__xi3_annotated.png"))

    doc_text = read(doc_path, String)
    @test occursin("<!-- AUX_FIGURES:BEGIN -->", doc_text)
    @test occursin("Delta = M_thr - M_mes", doc_text)

    @test isfile(joinpath(fig_dir, "fig5_taylor_decomposition.png"))
end
