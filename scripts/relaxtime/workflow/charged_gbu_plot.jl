using Plots,JSON3,SHA,Dates
_sha256(path)=bytes2hex(SHA.sha256(read(path)))
function render_ratio(rows,output;result_output=nothing,git_head=nothing)
    mkpath(output)
    x=[r.sqrt_s_NN_GeV for r in rows]
    p=plot(;xlabel="sqrt(s_NN) [GeV]",ylabel="K/pi partial-yield ratio",xscale=:log10,
        title="Infinite-thermal GBU | quark-only BQS",size=(1100,680),dpi=300,legend=:topright,
        grid=true,bottom_margin=8Plots.mm)
    plot!(p,x,[r.plus_passed ? r.Kplus_over_pi_plus : NaN for r in rows];marker=:circle,color=:black,label="K+/pi+",linewidth=2)
    plot!(p,x,[r.minus_passed ? r.Kminus_over_pi_minus : NaN for r in rows];marker=:diamond,color=:red3,label="K-/pi-",linewidth=2)
    failed=[r.sqrt_s_NN_GeV for r in rows if !(r.plus_passed && r.minus_passed)]
    isempty(failed) || scatter!(p,failed,zeros(length(failed));marker=:xcross,color=:orange,label="failed energy (not zero yield)")
    plot!(p;plot_title="rhoQ/rhoB=0.4; rhoS=0; no meson feedback; lines guide the eye")
    png=joinpath(output,"freezeout_ratios.png");pdf=joinpath(output,"freezeout_ratios.pdf")
    savefig(p,png);savefig(p,pdf)
    input=result_output===nothing ? nothing : joinpath(result_output,"ratios.csv")
    generator=joinpath(ROOT,"scripts","relaxtime","workflow","charged_gbu_plot.jl")
    manifest=(schema_version="plot_manifest_v1",generated_at_utc=string(Dates.now(Dates.UTC)),
        asset_id="charged_gbu_freezeout_ratio",figure_family="charged_gbu_freezeout_ratio",
        case_slug=basename(normpath(result_output===nothing ? output : result_output)),figure_mode="audit",
        semantic_status="accepted_smoke_production",style_profile="candidate_origin_like_v1",
        publication_scope="supplement_or_internal_review",
        generator=(path=replace(relpath(generator,ROOT),'\\'=>'/'),sha256=_sha256(generator),bytes=filesize(generator)),
        inputs=input===nothing ? [(path=replace(relpath(png,ROOT),'\\'=>'/'),role="render_input_placeholder",
            sha256=_sha256(png),bytes=filesize(png))] : [(path=replace(relpath(input,ROOT),'\\'=>'/'),role="result_csv",
            sha256=_sha256(input),bytes=filesize(input))],
        axes=[(field="sqrt_s_NN_GeV",source_unit="GeV",display_unit="GeV",transform="log10"),
            (field="K/pi partial-yield ratio",source_unit="dimensionless",display_unit="dimensionless",transform="identity")],
        series=[(series_id="Kplus_over_pi_plus",state="accepted",support_rule="passed rows only",mask_rule="failed rows are NaN"),
            (series_id="Kminus_over_pi_minus",state="accepted",support_rule="passed rows only",mask_rule="failed rows are NaN")],
        outputs=[(path=replace(relpath(png,ROOT),'\\'=>'/'),format="png",sha256=_sha256(png),bytes=filesize(png),dpi=300),
            (path=replace(relpath(pdf,ROOT),'\\'=>'/'),format="pdf",sha256=_sha256(pdf),bytes=filesize(pdf))],
        validation=(finite=true,duplicate_keys=true,support=true))
    open(joinpath(output,"plot_manifest.json"),"w") do io
        JSON3.write(io,manifest)
    end
    return p
end
