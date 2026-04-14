const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(@__DIR__, "t190_imag_path_evidence.jl"))

const OUTDIR_T200 = raw"D:\Desktop\Temp\relaxtime_t200_window"
const OUT_DETAIL_T200 = joinpath(OUTDIR_T200, "t200_imag_path_evidence_detail.csv")
const OUT_SUMMARY_T200 = joinpath(OUTDIR_T200, "t200_imag_path_evidence_summary.csv")

function main_t200()
    ensure_outdir(OUTDIR_T200)
    T_MeV = 200.0
    muB_MeV = 0.0
    ds = 1.0e-8

    xi_simple = [-0.34, -0.32, -0.30, -0.22, -0.20, 0.34, 0.36, 0.38, 0.40]
    xi_mixed = [-0.34, -0.32, -0.30, -0.22, -0.20, 0.34, 0.36, 0.38, 0.40]

    open(OUT_DETAIL_T200, "w") do io
        println(io, "branch,process,xi,delta_s,k0_s,k_norm,lambda_branch,term,has_pole,E0,Emin,Emax,term_im,den_re,den_im,pi_im,b0_re_from_pi,b0_im_from_pi,b0_im_used")
        for xi in xi_simple
            write_simple_rows!(io, T_MeV, muB_MeV, xi, ds)
        end
        for xi in xi_mixed
            write_mixed_rows!(io, T_MeV, muB_MeV, xi, ds)
        end
    end

    build_summary(OUT_DETAIL_T200, OUT_SUMMARY_T200)
    println("Wrote detail:  $OUT_DETAIL_T200")
    println("Wrote summary: $OUT_SUMMARY_T200")
end

main_t200()
