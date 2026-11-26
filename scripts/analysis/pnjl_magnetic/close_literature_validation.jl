#!/usr/bin/env julia
# Final targeted checks for the literature cross-validation. Reuses the
# independently implemented formulas; no full phase diagram is generated.
using JSON, Dates, LinearAlgebra
include(joinpath(@__DIR__, "crosscheck_literature.jl"))

function independent_roots(points, muB, T, eB, params, controls)
    seeds = Vector{Float64}[]
    for (a, b) in zip(points[1:end-1], points[2:end])
        lo, hi = a["muB_MeV"], b["muB_MeV"]
        lo == hi && continue
        min(lo, hi) <= muB <= max(lo, hi) || continue
        w = (muB - lo)/(hi - lo)
        push!(seeds, (1-w).*Float64.(a["x_state"]) + w.*Float64.(b["x_state"]))
    end
    roots = Any[]
    failures = Any[]
    mu = fill(muB/(3H), 3)
    for seed in seeds
        r = solve_fd(x -> MR.evaluate(x,T,mu,eB,params,controls).residual, seed;
            ftol=2e-12, h=2e-6)
        if !r.converged
            push!(failures, (;seed, residual=r.residual))
            continue
        end
        any(norm(r.x-q.x_state) < 1e-7 for q in roots) && continue
        ev = MR.evaluate(r.x,T,mu,eB,params,controls)
        push!(roots, (;x_state=r.x, residual=r.residual, omega=ev.omega,
            rho=ev.rho, masses_MeV=ev.masses.*H))
    end
    (;roots, failures)
end

function bidirectional_match(a, b)
    distances_a = [isempty(b) ? Inf : minimum(norm(x-y) for y in b) for x in a]
    distances_b = [isempty(a) ? Inf : minimum(norm(x-y) for x in a) for y in b]
    nearest = [isempty(b) ? 0 : argmin([norm(x-y) for y in b]) for x in a]
    (;count_a=length(a), count_b=length(b), distances_a, distances_b,
        bijective=length(a)==length(b)==length(unique(nearest)),
        max_distance=max(maximum(distances_a;init=0.),maximum(distances_b;init=0.)))
end

function branch_checks(folder)
    prod = JSON.parsefile(joinpath(folder,"corrected_palc_p128.json"))
    lower = JSON.parsefile(joinpath(folder,"corrected_palc_p64.json"))
    ref = JSON.parsefile(joinpath(folder,"resolved_mfir_o48_fine_v2.json"))
    rows = Any[]
    for c in prod["cases"]
        ebG = c["eB_GeV2"]; eb = ebG*(1000/H)^2; T = c["T_MeV"]/H
        params = MR.Parameters()
        controls = MR.Controls(order=96, n_max=128)
        rc = only(filter(q -> q["eB_GeV2"]==ebG,ref["cases"]))
        lc = filter(q -> q["eB_GeV2"]==ebG,lower["cases"])
        points = Any[]
        for pt in c["points"]
            x = Float64.(pt["x_state"]); mu = fill(pt["muB_MeV"]/(3H),3)
            ev = MR.evaluate(x,T,mu,eb,params,controls)
            # Raw PALC vertices can have tiny negative Polyakov fields at
            # the vacuum endpoint. Keep that observation instead of relabelling
            # the raw states as accepted public-solver candidates.
            clamped = copy(x); clamped[4:5] = clamp.(x[4:5],0.,1.)
            cev = MR.evaluate(clamped,T,mu,eb,params,controls)
            push!(points, (;index=pt["index"], muB_MeV=pt["muB_MeV"],
                independent_gap=norm(ev.condensate_gap),
                independent_stationarity=norm(ev.stationarity),
                density_delta=norm(ev.rho-Float64.(pt["density_inv_fm3"])),
                omega_delta=ev.omega-pt["omega"],
                min_mass_jacobian_singular_value=ev.mass_jacobian_sigma_min,
                polyakov_boundary_violation=max(0.,-minimum(x[4:5]),maximum(x[4:5])-1),
                boundary_projection_residual=norm(cev.residual),
                boundary_projection_density_delta=norm(cev.rho-ev.rho)))
        end
        targets = Any[]
        for t in c["targets"]
            muB = t["muB_MeV"]
            ir = independent_roots(rc["points"],muB,T,eb,params,controls)
            a = [Float64.(r["x_state"]) for r in t["palc_roots"]]
            b = [r.x_state for r in ir.roots]
            low_match = if isempty(lc)
                nothing
            else
                lt = only(filter(q -> q["muB_MeV"]==muB,only(lc)["targets"]))
                bidirectional_match(a,[Float64.(r["x_state"]) for r in lt["palc_roots"]])
            end
            push!(targets, (;muB_MeV=muB, reference_match=bidirectional_match(a,b),
                p64_match=low_match, independent_roots=ir.roots, failures=ir.failures))
        end
        println("branch eB=",ebG," vertices=",length(points),
            " max independent gap=",maximum(p.independent_gap for p in points),
            " target mismatch=",maximum(t.reference_match.max_distance for t in targets));flush(stdout)
        push!(rows,(;eB_GeV2=ebG,T_MeV=c["T_MeV"],
            production_order=c["p_num"],reference_order=controls.order,
            reference_n_max=controls.n_max,points,targets))
    end
    rows
end

function cold_gradient(phi, mu, eB, params)
    M = MR.masses(phi,eB,params)
    loop = cold_eval(M,mu,eB,params;with_omega=false).condensates
    MR.mass_jacobian(phi,eB,params)'*(loop-phi)
end

function cold_stability_checks(folder)
    source = JSON.parsefile(joinpath(folder,"literature_cold.json"))
    refined = JSON.parsefile(joinpath(folder,"literature_cold_refined.json"))
    jobs = Any[]
    for c in source["cases"], co in c["target"]["coexistences"]
        co["converged"] || continue
        any(j.id==c["id"] && abs(j.co["muq_MeV"]-co["muq_MeV"])<1e-6 for j in jobs) && continue
        push!(jobs,(;id=c["id"], set=c["set"], eB_GeV2=c["eB_GeV2"],co))
    end
    push!(jobs,(;id="Grunfeld2014_set1_B001",set=1,eB_GeV2=.01,co=first(refined["crossings"])))
    rows = Any[]
    for j in jobs
        p = literature_params(set=j.set); eb=j.eB_GeV2*(1000/H)^2
        mu = j.co["muq_MeV"]/H
        phases = j.co["phases"]
        seeds = [Float64.(q["masses_MeV"])/H for q in phases]
        for c in source["cases"]
            c["id"]==j.id || continue
            append!(seeds,[Float64.(r["masses"]) for r in c["target"]["catalogue"]["roots"]])
        end
        append!(seeds,[(1-w).*seeds[1]+w.*seeds[2] for w in (.25,.5,.75)])
        eigenchecks = Any[]
        for phase in phases
            phi = Float64.(phase["condensates"])
            checks = [begin
                hess = jacobian_fd(z->cold_gradient(z,mu,eb,p),phi;h)
                (;h, eigenvalues=eigvals(Symmetric((hess+hess')/2)),
                    symmetry_error=norm(hess-hess'))
            end for h in (2e-5,5e-6,1e-6)]
            push!(eigenchecks,(;masses_MeV=phase["masses_MeV"],checks))
        end
        sidechecks = Any[]
        for offset in (-.01,0.,.01)
            mus = mu+offset/H
            a=solve_fd(M->cold_eval(M,mus,eb,p;with_omega=false).gap,seeds[1];h=2e-7,ftol=1e-12)
            b=solve_fd(M->cold_eval(M,mus,eb,p;with_omega=false).gap,seeds[2];h=2e-7,ftol=1e-12)
            ea=cold_eval(a.x,mus,eb,p);ebval=cold_eval(b.x,mus,eb,p)
            alternatives=Any[];failed=0
            for seed in seeds
                r=solve_fd(M->cold_eval(M,mus,eb,p;with_omega=false).gap,seed;h=2e-7,ftol=1e-12)
                if !r.converged || minimum(r.x)<=0
                    failed+=1;continue
                end
                any(norm(r.x-q.masses)<1e-7 for q in alternatives) && continue
                ev=cold_eval(r.x,mus,eb,p)
                push!(alternatives,(;masses=r.x,omega=ev.omega,residual=r.residual))
            end
            push!(sidechecks,(;offset_muq_MeV=offset,phase_a_converged=a.converged,
                phase_b_converged=b.converged,omega_a_minus_b=ea.omega-ebval.omega,
                rhoB_a=sum(ea.rho)/3,rhoB_b=sum(ebval.rho)/3,
                phase_separation_MeV=norm(a.x-b.x)*H,
                known_candidate_count=length(alternatives),failed_seed_attempts=failed,
                selected_pair_above_lowest_known_omega=min(ea.omega,ebval.omega)-minimum(q.omega for q in alternatives)))
        end
        println("cold ",j.id," muq=",mu*H," Hessian minima=",
            [minimum(e.checks[end].eigenvalues) for e in eigenchecks],
            " pressure crossing=",[s.omega_a_minus_b for s in sidechecks]);flush(stdout)
        push!(rows,(;id=j.id,eB_GeV2=j.eB_GeV2,muq_MeV=mu*H,eigenchecks,sidechecks))
    end
    rows
end

function cep_shape_checks(folder)
    source=JSON.parsefile(joinpath(folder,"literature_cep.json"))
    # One light CEP and the lowest-temperature strange CEP, selected before
    # this shape check. Finite neighborhoods supplement the cusp equations.
    selected = [first(filter(c->c["target"]["sector"]=="light",source["cases"])),
        argmin(c->last(c["checks"])["T_MeV"],source["cases"])]
    rows=Any[]
    for c in selected
        t=c["target"];r=last(c["checks"]);p=literature_params(imc=t["imc"])
        eb=t["eB"]*(1000/H)^2;ctl=MR.Controls(order=64,n_max=240)
        k=t["sector"]=="light" ? 1 : 3
        z=vcat(Float64.(r["x_state"]),r["muq"]);s=z[k]
        checks=Any[]
        for dT in (-.5,.5)
            values=Any[]
            for ds in range(-.20,.20;length=17)
                q=cep_derivatives(s+ds,(r["T_MeV"]+dT)/H,eb,p,ctl,k,z;delta=.002)
                push!(values,(;condensate_offset=ds,dmu_dphi=q.dmu,
                    muB_MeV=3H*q.center.x[6],residual=q.center.residual))
            end
            sign_changes=count(i->signbit(values[i].dmu_dphi)!=signbit(values[i+1].dmu_dphi),1:length(values)-1)
            push!(checks,(;temperature_offset_MeV=dT,sign_changes,values))
        end
        println("CEP shape ",t["id"]," sign changes below/above=",[q.sign_changes for q in checks]);flush(stdout)
        push!(rows,(;id=t["id"],checks))
    end
    rows
end

function main_close(folder,output)
    result=Dict{String,Any}("generated_at"=>string(now()),"julia_version"=>string(VERSION),
        "method"=>"independent condensate gap on complete saved branches; bidirectional root matching; T=0 Hessians and pressure exchange; local CEP shape")
    save()=open(io->JSON.print(io,result,2),output,"w")
    result["branches"]=branch_checks(folder);save()
    result["cold_stability"]=cold_stability_checks(folder);save()
    result["cep_shapes"]=cep_shape_checks(folder);save()
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==2 || error("usage: close_literature_validation.jl INPUT_FOLDER OUTPUT.json")
    main_close(ARGS[1],ARGS[2])
end
