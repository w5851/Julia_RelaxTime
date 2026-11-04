#!/usr/bin/env julia
# SOURCE_ROOT INPUT_DIR OUTPUT_DIR; quote-remote only.
const MOTT_REFINEMENT_LIBRARY_PATH = joinpath(@__DIR__,
    get(ENV, "MOTT_REFINEMENT_LIBRARY", "run_transport_other_mott_refinement.jl"))
include(MOTT_REFINEMENT_LIBRARY_PATH)
profile=SP_OTHER.load_profile(joinpath(INPUT_DIR,"base.toml"))
thresholds,backgrounds=NamedTuple[],NamedTuple[]

function screen_point(row,state,segment,kind)
    for meson in MESONS
        cuts=meson in MIXED ? (:uu,:ss) : (:pair,)
        for cut in cuts
            value=endpoint_other(state,meson,cut;n_nodes=128)
            push!(thresholds,merge(point_fields(row),(continuation_segment=segment,kind=kind,
                meson=String(meson),cut=String(cut),inverse_re=real(value),inverse_im=imag(value),
                inverse_abs=abs(value),finite=isfinite(real(value))&&isfinite(imag(value)))))
        end
    end
end

for source in SOURCE_OTHER
    state=frozen_state(source,profile,32,16)
    screen_point(source,state,String(source.continuation_segment),"frozen_v15")
end
for segment in unique(String(r.continuation_segment) for r in SOURCE_OTHER)
    scan=sort(filter(r->String(r.continuation_segment)==segment,SOURCE_OTHER);by=r->Float64(r.xi))
    target=(curve_id=segment,window=String(first(scan).window),T_MeV=Float64(first(scan).T_MeV),muB_MeV=Float64(first(scan).muB_MeV))
    seed=original_seed_other(SP_OTHER.profile_model(profile),first(scan))
    started=time()
    for (index,(a,b)) in enumerate(zip(scan,scan[2:end]))
        Float64(b.xi)-Float64(a.xi)>.006 || continue
        xi=(Float64(a.xi)+Float64(b.xi))/2
        bg=background_other(target,xi,profile,"base";seed=seed)
        seed=bg.equilibrium.x_state
        row=merge(bg.row,(point_id=segment*"_mid"*lpad(string(index),3,'0'),))
        push!(backgrounds,merge(row,(continuation_segment=segment,point_left=String(a.point_id),point_right=String(b.point_id))))
        screen_point(row,bg.state,segment,"solved_midpoint")
    end
    CSV.write(joinpath(OUTPUT_DIR,"backgrounds.csv"),backgrounds)
    CSV.write(joinpath(OUTPUT_DIR,"thresholds.csv"),thresholds)
    println("MIDPOINT_SEGMENT ",segment," new_total=",length(backgrounds)," seconds=",round(time()-started;digits=1));flush(stdout)
end
write_json(joinpath(OUTPUT_DIR,"summary.json"),(status="completed_diagnostic",original_points=length(SOURCE_OTHER),
    solved_midpoints=length(backgrounds),threshold_rows=length(thresholds),
    finite_threshold_rows=count(r->r.finite,thresholds),background_residual_max=maximum(r.background_residual for r in backgrounds),
    calculation_sha="2e5c411c716d5cbc524c3fc6cfb335fe05b39444",production_kernel_changed=false,
    input_rows_sha256=digest(joinpath(INPUT_DIR,"source_rows.csv")),script_sha256=digest(@__FILE__),
    scope="New stable backgrounds at every continuous interval wider than 0.006; no interpolation of background or transport."))
