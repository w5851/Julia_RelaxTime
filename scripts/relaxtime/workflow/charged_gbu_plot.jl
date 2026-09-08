using Plots
function render_ratio(rows,output)
    x=[r.sqrt_s_NN_GeV for r in rows]
    p=plot(;xlabel="sqrt(s_NN) [GeV]",ylabel="K/pi partial-yield ratio",xscale=:log10,
        title="Infinite-thermal GBU | quark-only BQS",size=(1100,680),legend=:topright,
        grid=true,bottom_margin=8Plots.mm)
    plot!(p,x,[r.plus_passed ? r.Kplus_over_pi_plus : NaN for r in rows];marker=:circle,color=:black,label="K+/pi+",linewidth=2)
    plot!(p,x,[r.minus_passed ? r.Kminus_over_pi_minus : NaN for r in rows];marker=:diamond,color=:red3,label="K-/pi-",linewidth=2)
    failed=[r.sqrt_s_NN_GeV for r in rows if !(r.plus_passed && r.minus_passed)]
    isempty(failed) || scatter!(p,failed,zeros(length(failed));marker=:xcross,color=:orange,label="failed energy (not zero yield)")
    plot!(p;plot_title="rhoQ/rhoB=0.4; rhoS=0; no meson feedback; lines guide the eye")
    savefig(p,joinpath(output,"freezeout_ratios.png"));savefig(p,joinpath(output,"freezeout_ratios.pdf"))
    return p
end
