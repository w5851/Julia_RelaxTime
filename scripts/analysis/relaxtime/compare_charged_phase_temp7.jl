"""Read-only curve comparison across explicitly different backgrounds/conventions.
The overlay is display-only. No cross-background error metric or density is computed.
"""
module ChargedPhaseTemp7Comparison
using CSV, JSON3, SHA
ENV["GKSwstype"] = get(ENV,"GKSwstype","100")
using Plots
const ROOT = normpath(joinpath(@__DIR__,"..","..",".."))

function main()
    temp = get(ENV,"CHARGED_TEMP7_INPUT","")
    isempty(temp) && error("CHARGED_TEMP7_INPUT must explicitly name the historical audit directory")
    input = joinpath(ROOT,"data","outputs","results","relaxtime","analysis",
        "charged_rpa_phase_backend","negative_density_phase_fig2_like")
    output = get(ENV,"CHARGED_TEMP7_OUTPUT",joinpath(dirname(input),"temp7_controlled_overlay_20260905"))
    isdir(output) && !isempty(readdir(output)) && error("refusing to overwrite $(output)")
    old_csv = joinpath(temp,"output","fig2_paper2020_fold0pi_eta0_pv_b0","fig2_phase_shift_curves.csv")
    new_csv = joinpath(input,"charged_phase_profile_detail.csv")
    old = collect(CSV.File(old_csv))
    new = collect(CSV.File(new_csv))
    manifest = JSON3.read(read(joinpath(input,"plot_manifest.json"),String))
    all(r->r.T_MeV==90.0,old) || error("unexpected temp7 temperature")
    panels = []
    metrics = []
    fields = ((:phase_principal,:raw_phase,"principal"),
              (:phase_unwrapped,:unwrapped_phase,"unwrapped (different algorithms)"),
              (:phase_display,:phase_display_fold_0_pi,"fold_0_pi (display only)"))
    for (old_field,new_field,title) in fields
        p = plot(xlabel="external omega [GeV]",ylabel="phase / pi",title=title,
            legend=:bottomright,xlims=(0,1.6),linewidth=1.8,gridalpha=0.15,
            tickfontsize=9,guidefontsize=10,titlefontsize=11,legendfontsize=8)
        for (channel,pair,color,style) in (("K_plus","(s,u)",:firebrick,:solid),("K_minus","(u,s)",:darkorange,:dash))
            r = filter(r->r.meson==channel && r.muq_MeV==350.0,old)
            plot!(p,[x.omega_MeV/1000 for x in r],[getproperty(x,old_field)/π for x in r];
                label="temp7 $(channel) $(pair), q=0",color=color,linestyle=style)
        end
        for (q,color,style) in ((0.0,:royalblue,:solid),(1.0,:forestgreen,:dash))
            r = filter(r->r.channel=="K_plus" && r.variant=="pv_cut" && r.q_inv_fm==q,new)
            plot!(p,[x.omega_MeV/1000 for x in r],[getproperty(x,new_field)/π for x in r];
                label="BQS K_plus (u,s), q=$(q) fm^-1",color=color,linestyle=style)
        end
        push!(panels,p)
    end
    for channel in ("K_plus","K_minus"), muq in (300.0,350.0)
        r = filter(r->r.meson==channel && r.muq_MeV==muq,old)
        push!(metrics,(channel=channel,muq_MeV=muq,rows=length(r),
            max_principal_unwrapped_difference=maximum(abs(x.phase_principal-x.phase_unwrapped) for x in r)))
    end
    fig = plot(panels...;layout=(1,3),size=(1800,590),margin=5Plots.mm,
        plot_title="DIAGNOSTIC topology only: temp7 T=90, mu(u,d,s)=(350,350,70) MeV; BQS T=170, muB=240 MeV\nDifferent flavor labels / PV cut prescriptions; no fitted shifts, no density branch selected",
        plot_titlefontsize=11)
    mkpath(output)
    savefig(fig,joinpath(output,"temp7_bqs_phase_overlay.png"))
    record = Dict("status"=>"diagnostic_only_not_production","solver_called"=>false,
        "temp7_csv"=>old_csv,"temp7_sha256"=>bytes2hex(sha256(read(old_csv))),
        "current_csv"=>new_csv,"current_sha256"=>bytes2hex(sha256(read(new_csv))),
        "temp7_background"=>Dict("T_MeV"=>90,"mu_u_MeV"=>350,"mu_d_MeV"=>350,"mu_s_MeV"=>70,"q_inv_fm"=>0),
        "current_background"=>manifest.background,"current_mu_inv_fm"=>manifest.chemical_potentials_inv_fm,
        "flavor_conventions"=>"temp7 paper_2020: K+=(s,u), K-=(u,s); current K+=(u,s)",
        "phase_conventions"=>"temp7 atan(-Im inverse_from_legacy_B0,Re), forward unwrap; current -arg(inverse_PV_ret), reverse unwrap",
        "density_use"=>false,"cross_background_parity_claim"=>false,"temp7_metrics"=>metrics,
        "source_sha256"=>bytes2hex(sha256(read(@__FILE__))))
    open(joinpath(output,"manifest.json"),"w") do io
        JSON3.write(io,record)
    end
    println("[temp7-comparison] $(output)")
end
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
end
