"""Automatic pre-production closure audit; no density or background scan.

Compare identical domains in two loop coordinates, spectral/PV/UHP values,
contact moments, signed gap counts and vacuum/thermal regulator dependence.
All failed cases remain in the output. This does not select a regulator.
"""
module CausalGBUFullClosure
include("causal_gbu_coordinate_oracles.jl")
include("causal_gbu_gap_completeness.jl")
const O=CausalGBUCoordinateOracles
const G=CausalGBUGapCompleteness
const C=O.C
const R=O.R
using CSV,JSON3
using Main.RelaxTime.OneLoopIntegrals: A, B0_retarded, B0_spectral_cut

function settings_from_env()
    np=parse(Int,get(ENV,"GBU_FULL_LOOP_NODES","512"))
    mesh=parse(Int,get(ENV,"GBU_FULL_MESH","512"))
    np>=16 && mesh>=16 && iseven(np) && iseven(mesh) || throw(ArgumentError("even nodes>=16 required"))
    qraw=get(ENV,"GBU_FULL_Q","density_shells")
    q=qraw=="density_shells" ? sort!(unique(vcat([0.,0.01,0.02,1.,3.,6.],first(R.gauleg(0.,8.,24))))) :
        sort!(unique(parse.(Float64,split(qraw,','))))
    !isempty(q) && first(q)==0 && all(x->isfinite(x) && x>=0,q) || throw(ArgumentError("q list must include 0"))
    return (np=np,nx=div(np,2),mesh=mesh,ne=128,nr=128,thermal=24.,qs=q)
end

function direct_values(bg,ch,q,s,regulator,zs)
    a,b=R.charged_rpa_spec(ch).pair
    args=(q,bg.m[a],bg.mu[a],bg.m[b],bg.mu[b],bg.T)
    if regulator===:two_line
        grid=R.build_spectral_bubble(args...;Phi=bg.Phi,PhiBar=bg.PhiBar,
            vacuum_cutoff_inv_fm=bg.vacuum,thermal_cutoff_inv_fm=s.thermal,
            momentum_nodes=s.np,angle_nodes=s.nx)
        vals=[R.spectral_bubble(grid,real(z)-bg.mu[a]+bg.mu[b];eta_inv_fm=imag(z)) for z in zs]
        return (values=getproperty.(vals,:value),b0=getproperty.(vals,:B0),contact=grid.contact_inv_fm2,
            identity=maximum(v.contact_identity_residual for v in vals),grid=grid)
    end
    grid=C.centered_atoms(args...,bg.Phi,bg.PhiBar,bg.vacuum,s.thermal;np=s.np,nx=s.nx)
    vals=[C.atom_value(grid,z) for z in zs]
    return (values=getproperty.(vals,:value),b0=getproperty.(vals,:b0),contact=grid.contact,
        identity=maximum(v.contact_residual for v in vals),grid=grid)
end

function cut_value(bg,ch,q,l,s,regulator;reverse=false)
    a,b=R.charged_rpa_spec(ch).pair
    reverse && ((a,b)=(b,a))
    total=0.
    for (component,L) in ((:vacuum,bg.vacuum),(:thermal,s.thermal))
        total+=regulator===:two_line ?
            B0_spectral_cut(l,q,bg.m[a],bg.mu[a],bg.m[b],bg.mu[b],bg.T;
                Φ=bg.Phi,Φbar=bg.PhiBar,pmax_inv_fm=L,energy_nodes=s.ne,component=component).imaginary :
            C.centered_cut(l,q,bg.m[a],bg.mu[a],bg.m[b],bg.mu[b],bg.T,bg.Phi,bg.PhiBar,L;
                ne=s.ne,component=component)
    end
    return C.NC/(8pi^2)*(l^2-q^2-(bg.m[a]-bg.m[b])^2)*total
end

function run_case(bg,ch,q,regulator,s,q0vac)
    a,b=R.charged_rpa_spec(ch).pair
    m1,m2,mu1,mu2=bg.m[a],bg.m[b],bg.mu[a],bg.mu[b]
    shift=mu1-mu2
    th=hypot(q,m1+m2)
    z0=2.0+0.6im
    zboost=sqrt(z0^2+q^2)
    zs=[shift+0.4im,1.0+0.6im,-1.0+0.6im,th+0.5+0.4im,-th-0.5+0.4im,12.0+0.8im,zboost]
    atom=direct_values(bg,ch,q,s,regulator,zs)
    args=(q,m1,mu1,m2,mu2,bg.T,bg.Phi,bg.PhiBar,bg.vacuum,s.thermal,zs,regulator)
    radial=O.radial_loop(args...;np=s.np,nr=s.nx)
    coarse=O.radial_loop(args...;np=div(s.np,2),nr=div(s.nx,2))
    loop_difference=maximum(abs.(radial.values.-atom.values))
    node_difference=maximum(abs.(radial.values.-coarse.values))
    contact_difference=abs(radial.contact-atom.contact)
    moment=max(radial.vacuum.moment_identity_residual,radial.thermal.moment_identity_residual)
    identity=max(atom.identity,radial.vacuum.contact_identity_residual,radial.thermal.contact_identity_residual)
    complexrows,rootrows,meshrows=NamedTuple[],NamedTuple[],NamedTuple[]
    cutrows=NamedTuple[]
    signedcounts=Int[]
    roots_by_mesh=Vector{Vector{Float64}}()
    representation_differences=Float64[]
    for mesh in (div(s.mesh,2),s.mesh)
        p=if regulator===:two_line
            R.build_bubble_dispersion(atom.grid;segment_nodes=mesh,energy_nodes=s.ne)
        else
            C.centered_profile(bg,ch,q;mesh=mesh,ne=s.ne).profile
        end
        spectral=[R.cauchy_transform(p,z) for z in zs]
        repr=maximum(abs.(spectral.-radial.values))
        push!(representation_differences,repr)
        for (i,z) in enumerate(zs)
            push!(complexrows,(channel=String(ch),q_inv_fm=q,regulator=String(regulator),mesh=mesh,
                probe=i,lambda_real_inv_fm=real(z),eta_inv_fm=imag(z),
                pi_radial_real=real(radial.values[i]),pi_radial_imag=imag(radial.values[i]),
                pi_original_real=real(atom.values[i]),pi_original_imag=imag(atom.values[i]),
                pi_spectral_real=real(spectral[i]),pi_spectral_imag=imag(spectral[i]),
                coordinate_difference=abs(radial.values[i]-atom.values[i]),
                radial_node_difference=abs(radial.values[i]-coarse.values[i]),
                spectral_difference=abs(spectral[i]-radial.values[i]),
                pi_vacuum_real=real(radial.vacuum.values[i]),pi_vacuum_imag=imag(radial.vacuum.values[i]),
                pi_thermal_real=real(radial.thermal.values[i]),pi_thermal_imag=imag(radial.thermal.values[i]),
                production_authorized=false))
        end
        cuts=C.support_intervals(q,m1,m2,bg.vacuum,s.thermal,regulator)
        geometry=G.support_geometry_check(G.signed_cells(p),cuts)
        roots=G.count_all_real_gaps(p,bg.coupling[ch],q;root_nodes=s.nr)
        for r in roots.rows
            push!(rootrows,merge((channel=String(ch),q_inv_fm=q,regulator=String(regulator),mesh=mesh,
                shift_inv_fm=shift),r,(production_authorized=false,)))
        end
        ls=sort!([parse(Float64,l) for r in roots.rows for l in split(r.roots_lambda_inv_fm,';') if !isempty(l)])
        push!(roots_by_mesh,ls)
        push!(signedcounts,roots.count)
        inverse0=1-4bg.coupling[ch]*R.cauchy_transform(p,shift)
        static_ok=real(inverse0)>0 && abs(imag(inverse0))<1e-8
        max_imag_error,max_reciprocity=0.,0.
        # Non-knot probes on both signs, within Landau and pair support as well as gaps.
        probe_lambdas=unique(vcat([shift,0.,-0.3,0.3,-th*0.8,th*0.8,-th-0.4,th+0.4],
            [left+(right-left)*f for (left,right) in cuts for f in (0.17,0.43,0.79)]))
        for l in probe_lambdas
            directcut=cut_value(bg,ch,q,l,s,regulator)
            rev=cut_value(bg,ch,q,-l,s,regulator;reverse=true)
            sp=R.cauchy_transform(p,l)
            max_imag_error=max(max_imag_error,abs(imag(sp)-directcut))
            max_reciprocity=max(max_reciprocity,abs(directcut+rev))
            push!(cutrows,(channel=String(ch),q_inv_fm=q,regulator=String(regulator),mesh=mesh,
                lambda_inv_fm=l,k0_inv_fm=l-shift,pv_real=real(sp),interpolated_imag=imag(sp),
                direct_cut_imag=directcut,reversed_negative_cut_imag=rev,
                interpolation_error=abs(imag(sp)-directcut),reciprocity_error=abs(directcut+rev),
                production_authorized=false))
        end
        push!(meshrows,(channel=String(ch),q_inv_fm=q,regulator=String(regulator),mesh=mesh,
            signed_gap_count=length(roots.rows),signed_real_root_count=roots.count,
            positive_k0_roots=count(>(shift),ls),negative_k0_roots=count(<(shift),ls),
            all_interpolant_gaps_passed=roots.passed,geometry_passed=geometry.passed,
            support_excess_width_inv_fm=geometry.max_excess_width_inv_fm,
            support_missing_width_inv_fm=geometry.max_missing_width_inv_fm,
            support_roundoff_allowance_inv_fm=geometry.roundoff_allowance_inv_fm,
            tail_radius_lambda_inv_fm=roots.tail.radius_inv_fm,tail_inverse_deviation_bound=roots.tail.inverse_deviation_bound,
            tail_passed=roots.tail.passed,static_passed=static_ok,
            max_spectral_difference=repr,max_direct_cut_interpolation_error=max_imag_error,
            max_cut_reciprocity_error=max_reciprocity,
            scope=roots.scope,full_physical_spectrum_certified=false,production_authorized=false))
    end
    meshstable=signedcounts[1]==signedcounts[2]
    rootdrift=meshstable ? maximum(abs.(roots_by_mesh[1].-roots_by_mesh[2]);init=0.) : Inf
    lastmesh=last(meshrows)
    coordinate_ok=loop_difference<1e-6 && node_difference<1e-6 && contact_difference<1e-8
    algebra_ok=identity<1e-10 && moment<1e-10
    # Keep the established representation, loop and algebra targets unchanged.
    representation_ok=last(representation_differences)<1e-4
    gaps_ok=all(r.all_interpolant_gaps_passed && r.geometry_passed && r.tail_passed for r in meshrows) && meshstable
    reciprocity_ok=lastmesh.max_cut_reciprocity_error<1e-10
    cut_ok=lastmesh.max_direct_cut_interpolation_error<1e-4
    summary=(channel=String(ch),q_inv_fm=q,regulator=String(regulator),
        max_coordinate_difference=loop_difference,max_radial_node_difference=node_difference,
        contact_coordinate_difference=contact_difference,max_contact_identity_residual=identity,
        max_contact_moment_residual=moment,contact_inv_fm2=radial.contact,
        vacuum_contact_inv_fm2=radial.vacuum.contact,thermal_contact_inv_fm2=radial.thermal.contact,
        coarse_spectral_difference=first(representation_differences),fine_spectral_difference=last(representation_differences),
        mesh_root_counts_stable=meshstable,max_root_mesh_shift_inv_fm=rootdrift,
        signed_real_root_count=lastmesh.signed_real_root_count,
        positive_k0_roots=lastmesh.positive_k0_roots,negative_k0_roots=lastmesh.negative_k0_roots,
        vacuum_boost_violation_real=real(last(radial.vacuum.values)-q0vac),
        vacuum_boost_violation_imag=imag(last(radial.vacuum.values)-q0vac),
        coordinate_passed=coordinate_ok,algebra_passed=algebra_ok,representation_passed=representation_ok,
        all_real_gap_checks_passed=gaps_ok,reciprocity_passed=reciprocity_ok,
        cut_interpolation_passed=cut_ok,static_passed=lastmesh.static_passed,
        passed=coordinate_ok && algebra_ok && representation_ok && gaps_ok && reciprocity_ok && cut_ok && lastmesh.static_passed,
        density_computed=false,production_authorized=false)
    return (summary=summary,complexrows=complexrows,rootrows=rootrows,meshrows=meshrows,cutrows=cutrows)
end

function q0_audit(bg,ch,s)
    a,b=R.charged_rpa_spec(ch).pair
    zs=[0.3+0.4im,2.0+0.6im,5.0+0.8im]
    ps,ws=R.gauleg(0.,s.thermal,s.np)
    expectedA=A(bg.m[a],bg.mu[a],bg.T,bg.Phi,bg.PhiBar,ps,ws)+A(bg.m[b],bg.mu[b],bg.T,bg.Phi,bg.PhiBar,ps,ws)
    args=(0.,bg.m[a],bg.mu[a],bg.m[b],bg.mu[b],bg.T,bg.Phi,bg.PhiBar,bg.vacuum,s.thermal,zs)
    two=O.radial_loop(args...,:two_line;np=s.np,nr=s.nx)
    cen=O.radial_loop(args...,:centered;np=s.np,nr=s.nx)
    # Old q=0 B0 integrates both pieces to model Lambda. Match THAT domain only.
    hs=merge(s,(thermal=bg.vacuum,))
    hard=direct_values(bg,ch,0.,hs,:two_line,zs)
    old=[B0_retarded(real(z),0.,bg.m[a],bg.mu[a],bg.m[b],bg.mu[b],bg.T;
        Φ=bg.Phi,Φbar=bg.PhiBar,eta_inv_fm=imag(z),energy_nodes=s.np) for z in zs]
    old_error=maximum(abs.(old.-hard.b0))
    contacterror=abs(two.contact-expectedA)
    routeerror=maximum(abs.(two.values.-cen.values))
    row=(channel=String(ch),q_inv_fm=0.,q0_regulator_value_difference=routeerror,
        q0_contact_A_difference=contacterror,q0_all_hard_old_B0_difference=old_error,
        old_B0_comparison_scope="q0_same_all_hard_domain_only",
        passed=routeerror<1e-12 && contacterror<1e-10 && old_error<1e-8,production_authorized=false)
    return row,two.vacuum.values[2]
end

function main()
    s=settings_from_env()
    base=joinpath(R.ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    bg=R.frozen_background(joinpath(base,"negative_density_phase_fig2_like"))
    output=get(ENV,"GBU_FULL_OUTPUT",joinpath(base,"full_signed_regulator_closure_20260906"))
    hashes=R.start_output(output)
    summaries,complexrows,rootrows,meshrows,cutrows,q0rows,failures=(NamedTuple[] for _ in 1:7)
    save(name,rows)=isempty(rows) ? nothing : CSV.write(joinpath(output,name*".csv"),rows)
    vacuum_q0=Dict{Symbol,ComplexF64}()
    for ch in R.CHANNELS
        zero,vac=q0_audit(bg,ch,s)
        push!(q0rows,zero)
        vacuum_q0[ch]=vac
        save("q0_checks",q0rows)
    end
    # The serial pilot measured independent expensive cases. Four worker threads
    # share read-only inputs; only this task writes the checkpoint CSVs.
    # gauleg's standard-node cache is protected by its existing ReentrantLock.
    jobs=[(ch,q,regulator) for ch in R.CHANNELS for q in s.qs for regulator in (:two_line,:centered)]
    workers=min(4,Threads.nthreads())
    for firstjob in 1:workers:length(jobs)
        batch=jobs[firstjob:min(firstjob+workers-1,length(jobs))]
        tasks=map(batch) do (ch,q,regulator)
            Threads.@spawn begin
                started=time_ns()
                try
                    (result=run_case(bg,ch,q,regulator,s,vacuum_q0[ch]),error="",elapsed_s=(time_ns()-started)/1e9)
                catch err
                    err isa InterruptException && rethrow()
                    (result=nothing,error=sprint(showerror,err),elapsed_s=(time_ns()-started)/1e9)
                end
            end
        end
        for ((ch,q,regulator),task) in zip(batch,tasks)
            item=fetch(task)
            if item.result!==nothing
                result=item.result
                push!(summaries,merge(result.summary,(elapsed_s=item.elapsed_s,)))
                append!(complexrows,result.complexrows);append!(rootrows,result.rootrows)
                append!(meshrows,result.meshrows);append!(cutrows,result.cutrows)
                println("[gbu-full] $(ch) $(regulator) q=$(q) passed=$(result.summary.passed) roots=$(result.summary.signed_real_root_count) coordinate=$(result.summary.max_coordinate_difference)")
            else
                push!(failures,(channel=String(ch),q_inv_fm=q,regulator=String(regulator),reason=item.error))
                println(stderr,"[gbu-full-failed] $(ch) $(regulator) q=$(q): $(item.error)")
            end
            save("summary",summaries);save("complex_probes",complexrows);save("signed_gap_counts",rootrows)
            save("spectral_mesh_checks",meshrows);save("real_axis_cut_checks",cutrows);save("failures",failures)
            flush(stdout)
        end
    end
    R.finish_output(output,hashes,Dict("status"=>"preproduction_closure_audit_completed","background"=>bg,
        "settings"=>s,"worker_threads"=>workers,"expected_cases"=>8length(s.qs),"completed_cases"=>length(summaries),
        "failed_evaluations"=>length(failures),"passed_cases"=>count(r->r.passed,summaries),
        "q0_passed"=>all(r.passed for r in q0rows),"all_checks_passed"=>isempty(failures) &&
            length(summaries)==8length(s.qs) && all(r.passed for r in summaries) && all(r.passed for r in q0rows),
        "solver_called"=>false,"density_computed"=>false,"regulator_selected"=>false,
        "root_scope"=>"all_signed_real_gaps_of_each_compact_interpolant_with_endpoint_and_tail_checks",
        "full_physical_spectrum_certified"=>false,"limitations"=>[
            "Real-gap completeness is for each interpolant, not a proof for the exact continuum kernel",
            "Cut-embedded zeros, UHP stability and exact Mott endpoints are separate",
            "Vacuum decomposition keeps frozen in-medium masses, not a new equilibrium solution",
            "No finite-q equality between distinct regulators or with q0 extrapolation is imposed",
            "No density grid extension or production promotion is performed"]))
    println("[gbu-full-done] $(output)")
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
