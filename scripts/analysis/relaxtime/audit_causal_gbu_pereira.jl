"""Pereira formula parity and separate PNJL/thermal-tail diagnostics.

No equilibrium solve, density scan, branch choice, regulator switch, or
production default change. Use GBU_PEREIRA_OUTPUT for a NEW output directory.
GBU_PEREIRA_PDF optionally binds the user-supplied paper by SHA256.
"""
module CausalGBUPereiraAudit
include("causal_gbu_research_utils.jl")
include("causal_gbu_pereira_reference.jl")
const R=CausalGBUResearch
const P=CausalGBUPereiraReference
using CSV, JSON3

function reference(bg,c,q,zs,phi,bar,thermal,n)
    f1,f2=R.charged_rpa_spec(c).pair
    a,u,b,v=bg.m[f1],bg.mu[f1],bg.m[f2],bg.mu[f2]
    vac=P.cylindrical_reference(q,a,u,b,v,bg.T,bg.vacuum,zs;
        Phi=phi,PhiBar=bar,component=:vacuum,nz=n,ny=n)
    th=P.cylindrical_reference(q,a,u,b,v,bg.T,thermal,zs;
        Phi=phi,PhiBar=bar,component=:thermal,nz=n,ny=n)
    return (pi=vac.pi_p.+th.pi_p,b0=vac.b0.+th.b0,
        contact=vac.contact+th.contact,vacuum=vac,thermal=th)
end

function cut(bg,c,q,l,phi,bar,thermal,n)
    f1,f2=R.charged_rpa_spec(c).pair
    args=(l,q,bg.m[f1],bg.mu[f1],bg.m[f2],bg.mu[f2],bg.T)
    v=P.reference_cut(args...,bg.vacuum;Phi=phi,PhiBar=bar,component=:vacuum,nodes=n)
    t=P.reference_cut(args...,thermal;Phi=phi,PhiBar=bar,component=:thermal,nodes=n)
    return (pair=v.pair+t.pair,landau=v.landau+t.landau,imaginary=v.imaginary+t.imaginary)
end

function main(;qs=(0.,0.8,1.4,2.7,3.2,6.2),nodes=384)
    nodes>=16 && iseven(nodes) || throw(ArgumentError("even nodes >=16 required"))
    all(q->isfinite(q) && q>=0,qs) && !isempty(qs) || throw(ArgumentError("invalid q probes"))
    base=joinpath(R.ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    bg=R.frozen_background(joinpath(base,"negative_density_phase_fig2_like"))
    output=get(ENV,"GBU_PEREIRA_OUTPUT",joinpath(base,"pereira_formula_checks_20260906"))
    paper=get(ENV,"GBU_PEREIRA_PDF","")
    !isempty(paper) && !isfile(paper) && error("paper not found")
    hashes=R.start_output(output)
    paper_hash=isempty(paper) ? "" : R.hashfile(paper)
    variants=[("fermi_hard",1.,1.,bg.vacuum),("pnjl_hard",bg.Phi,bg.PhiBar,bg.vacuum)]
    append!(variants,[("pnjl_tail$(Int(L))",bg.Phi,bg.PhiBar,L) for L in (8.,12.,16.,24.)])
    zs=ComplexF64[-3.5+0.7im,-0.5+0.6im,0.5+0.6im,3.5+0.7im,10.0+1im]
    loops,cuts,summary,effects,weak,failures=NamedTuple[],NamedTuple[],NamedTuple[],NamedTuple[],NamedTuple[],NamedTuple[]
    tolerance=1e-6
    cut_tolerance=1e-8
    save(name,rows)=isempty(rows) ? nothing : CSV.write(joinpath(output,name*".csv"),rows)
    for c in R.CHANNELS,q in qs
        first_flavor,second_flavor=R.charged_rpa_spec(c).pair
        a,u,b,v=bg.m[first_flavor],bg.mu[first_flavor],bg.m[second_flavor],bg.mu[second_flavor]
        variant_results=Dict{String,Any}()
        for (name,phi,bar,Lth) in variants
            try
                ref=reference(bg,c,q,zs,phi,bar,Lth,nodes)
                coarse=reference(bg,c,q,zs,phi,bar,Lth,div(nodes,2))
                g=R.build_spectral_bubble(q,a,u,b,v,bg.T;Phi=phi,PhiBar=bar,
                    vacuum_cutoff_inv_fm=bg.vacuum,thermal_cutoff_inv_fm=Lth,
                    momentum_nodes=nodes,angle_nodes=nodes)
                small=R.build_spectral_bubble(q,a,u,b,v,bg.T;Phi=phi,PhiBar=bar,
                    vacuum_cutoff_inv_fm=bg.vacuum,thermal_cutoff_inv_fm=Lth,
                    momentum_nodes=div(nodes,2),angle_nodes=div(nodes,2))
                max_pi,max_b0,max_current_node=0.,0.,0.
                for (i,z) in enumerate(zs)
                    current=R.spectral_bubble(g,real(z)-u+v;eta_inv_fm=imag(z))
                    old=R.spectral_bubble(small,real(z)-u+v;eta_inv_fm=imag(z))
                    ep,eb=abs(current.value-ref.pi[i]),abs(current.B0-ref.b0[i])
                    en=max(abs(current.value-old.value),abs(current.B0-old.B0))
                    max_pi=max(max_pi,ep);max_b0=max(max_b0,eb);max_current_node=max(max_current_node,en)
                    push!(loops,(channel=String(c),variant=name,q_inv_fm=q,
                        lambda_real_inv_fm=real(z),eta_inv_fm=imag(z),
                        reference_pi_real=real(ref.pi[i]),reference_pi_imag=imag(ref.pi[i]),
                        current_pi_real=real(current.value),current_pi_imag=imag(current.value),
                        pi_error=ep,B0_error=eb,current_node_change=en,production_authorized=false))
                end
                ls=unique([u-v,0.,-0.5hypot(q,a-b),0.5hypot(q,a-b),
                           -hypot(q,a+b)-0.3,hypot(q,a+b)+0.3,-8.,8.])
                max_cut=0.
                for l in ls
                    ref_cut=cut(bg,c,q,l,phi,bar,Lth,128)
                    cv=Main.RelaxTime.OneLoopIntegrals.B0_spectral_cut(l,q,a,u,b,v,bg.T;
                        Φ=phi,Φbar=bar,pmax_inv_fm=bg.vacuum,component=:vacuum,energy_nodes=128)
                    ct=Main.RelaxTime.OneLoopIntegrals.B0_spectral_cut(l,q,a,u,b,v,bg.T;
                        Φ=phi,Φbar=bar,pmax_inv_fm=Lth,component=:thermal,energy_nodes=128)
                    ep,el=abs(ref_cut.pair-cv.pair-ct.pair),abs(ref_cut.landau-cv.landau-ct.landau)
                    max_cut=max(max_cut,ep,el)
                    push!(cuts,(channel=String(c),variant=name,q_inv_fm=q,lambda_inv_fm=l,
                        k0_inv_fm=l-u+v,reference_pair=ref_cut.pair,reference_landau=ref_cut.landau,
                        pair_error=ep,landau_error=el,production_authorized=false))
                end
                contact_error=abs(ref.contact-g.contact_inv_fm2)
                exact_contact=P.vacuum_A(q,a,bg.vacuum)+P.vacuum_A(q,b,bg.vacuum)
                analytic_error=abs(ref.vacuum.contact-exact_contact)
                reference_node=max(maximum(abs.(ref.pi.-coarse.pi)),maximum(abs.(ref.b0.-coarse.b0)),
                                   abs(ref.contact-coarse.contact))
                passed=max(max_pi,max_b0,max_current_node,contact_error,analytic_error,reference_node)<tolerance &&
                       max_cut<cut_tolerance
                push!(summary,(channel=String(c),variant=name,q_inv_fm=q,Phi=phi,PhiBar=bar,
                    thermal_cutoff_inv_fm=Lth,contact_inv_fm2=ref.contact,
                    vacuum_contact_inv_fm2=ref.vacuum.contact,thermal_contact_inv_fm2=ref.thermal.contact,
                    pi_error=max_pi,B0_error=max_b0,contact_error=contact_error,
                    analytic_vacuum_contact_error=analytic_error,reference_node_change=reference_node,
                    current_node_change=max_current_node,cut_error=max_cut,passed=passed,
                    equilibrium_recomputed=false,density_computed=false,production_authorized=false))
                variant_results[name]=ref
                println("[pereira] $(c) q=$(q) $(name) pass=$(passed) loop=$(max_pi) cut=$(max_cut) node=$(reference_node)")
            catch err
                err isa InterruptException && rethrow()
                push!(failures,(channel=String(c),variant=name,q_inv_fm=q,reason=sprint(showerror,err)))
                println(stderr,"[pereira-failed] $(c) q=$(q) $(name): $(sprint(showerror,err))")
            end
            save("complex_checks",loops);save("cut_checks",cuts);save("summary",summary);save("failures",failures)
            flush(stdout)
        end
        for (left,right,label) in (("fermi_hard","pnjl_hard","occupation_only"),
                ("pnjl_hard","pnjl_tail24","thermal_extension"),
                ("pnjl_tail8","pnjl_tail12","thermal_tail_8_to_12"),
                ("pnjl_tail12","pnjl_tail16","thermal_tail_12_to_16"),
                ("pnjl_tail16","pnjl_tail24","thermal_tail_16_to_24"))
            haskey(variant_results,left) && haskey(variant_results,right) || continue
            x,y=variant_results[left],variant_results[right]
            for (i,z) in enumerate(zs)
                d=y.pi[i]-x.pi[i]
                push!(effects,(channel=String(c),q_inv_fm=q,comparison=label,
                    lambda_real_inv_fm=real(z),eta_inv_fm=imag(z),pi_change_real=real(d),
                    pi_change_imag=imag(d),pi_change_abs=abs(d),contact_change=y.contact-x.contact,
                    density_effect_quantified=false,production_authorized=false))
            end
        end
        save("extension_effects",effects)
    end
    # Poisson-weighted real-axis integral vs the complex loop. This is a
    # distributional B0 test, NOT a GBU density eta-convergence certification.
    c,q=:K_plus,1.4
    f1,f2=R.charged_rpa_spec(c).pair
    a,b=bg.m[f1],bg.m[f2];width=0.8;center=2.0;L=24.
    radius=hypot(a,L)+hypot(b,L)
    knots=sort!(unique!(vcat(collect(range(-radius,radius;length=257)),
        [-hypot(q,a+b),-hypot(q,a-b),0.,hypot(q,a-b),hypot(q,a+b),
         -hypot(a,bg.vacuum)-hypot(b,bg.vacuum),hypot(a,bg.vacuum)+hypot(b,bg.vacuum)])))
    integral=0.
    for i in 1:length(knots)-1
        xs,ws=R.gauleg(knots[i],knots[i+1],16)
        for (x,w) in zip(xs,ws)
            integral+=w*width/(pi*((x-center)^2+width^2))*cut(bg,c,q,x,bg.Phi,bg.PhiBar,L,64).imaginary
        end
    end
    zseq=complex.(fill(center,7),width.+[0.,0.2,0.1,0.05,0.025,0.0125,0.00625])
    values=reference(bg,c,q,zseq,bg.Phi,bg.PhiBar,L,nodes).b0
    for i in eachindex(zseq)
        push!(weak,(channel=String(c),q_inv_fm=q,eta_inv_fm=imag(zseq[i])-width,
            test_width_inv_fm=width,pv_weighted_integral=integral,
            eta_weighted_integral=imag(values[i]),difference_from_pv=abs(imag(values[i])-integral),
            representation_error=abs(imag(values[1])-integral),production_authorized=false))
    end
    save("weak_retarded_check",weak)
    weak_ok=abs(imag(values[1])-integral)<1e-4 &&
            all(diff([abs(imag(values[i])-imag(values[1])) for i in 2:length(values)]).<0)
    pass=isempty(failures) && length(summary)==length(qs)*length(variants)*length(R.CHANNELS) &&
         all(r.passed for r in summary) && weak_ok
    R.finish_output(output,hashes,Dict("status"=>pass ? "formula_parity_passed" : "formula_parity_failed",
        "all_checks_passed"=>pass,"background"=>bg,"solver_called"=>false,"density_computed"=>false,
        "nodes"=>nodes,"coarse_nodes"=>div(nodes,2),"qs"=>collect(qs),"variants"=>variants,
        "comparison_tolerance"=>tolerance,"cut_tolerance"=>cut_tolerance,"weak_passed"=>weak_ok,
        "failed_evaluations"=>length(failures),"failed_cases"=>count(r->!r.passed,summary),
        "paper_doi"=>"10.1103/PhysRevC.109.025206","paper_path"=>paper,"paper_sha256"=>paper_hash,
        "paper_equations"=>["26","55","56","60","C18","C19","D2-D10"],
        "limitations"=>["Fixed masses and chemical potentials; fermi_hard is not a newly solved NJL equilibrium",
            "PNJL and thermal extension are project variants, not the paper's original all-hard NJL model",
            "Finite UHP probes and weighted B0 parity do not certify RPA stability or GBU density convergence",
            "No regulator selection, production promotion, or full freezeout result is authorized"],
        "production_authorized"=>false))
    println("[pereira] output=$(output) all_passed=$(pass)")
    pass || error("Pereira audit failed; evidence retained in $(output)")
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
