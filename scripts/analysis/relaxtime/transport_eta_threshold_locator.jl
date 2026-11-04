# Diagnostic threshold evaluation on the frozen background. This does not
# replace the historical B0 kernel used by the transport or pole solvers.
const PAIR_ETA = Main.RelaxTime.PairLoopQuadrature

function eta_pair_inverse(state, mass_offset; n_nodes=128)
    qp, tp, kc = state.qp, state.tp, state.kc
    s = (2qp.m.u + mass_offset)^2
    polarizations = map((:u, :s)) do flavor
        m, mu, av = getproperty(qp.m, flavor), getproperty(qp.μ, flavor), getproperty(qp.A, flavor)
        plan = PAIR_ETA.PairLoopPlan(m, mu, tp; n_nodes=n_nodes)
        b0 = PAIR_ETA.pair_b0_pv(plan, s)
        -Main.Constants_PNJL.N_color/(8pi^2) * (2av-s*b0)
    end
    elems = EC.mixing_matrix_elements(polarizations..., kc, :P)
    MP._mixed_inverse(:P, :light, elems)
end

function locate_all()
    profile = SP.load_profile(joinpath(INPUT_DIR, "base.toml"))
    histories, location_rows, comparisons, pole_candidates, pole_rows, background_controls =
        (NamedTuple[] for _ in 1:6)
    locations = Any[]
    tolerance = Float64(INPUT_CONFIG.xi_bisection_tolerance)
    function checkpoint()
        for (file, rows) in (("location_iterations",histories), ("location_brackets",location_rows),
                ("kernel_comparison",comparisons), ("pole_candidates",pole_candidates),
                ("poles",pole_rows), ("background_controls",background_controls))
            isempty(rows) || CSV.write(joinpath(OUTPUT_DIR,file*".csv"),rows)
        end
    end
    for target in INPUT_CONFIG.targets
        target_locations = Dict{String,Any}()
        original = only(r for r in SOURCE_ROWS if String(r.point_id)==String(target.point_left))
        seed = original_seed(SP.profile_model(profile), original)
        for level_name in ("base", "A64", "background_ref")
            level = numerical_level(profile, level_name)
            cache = Dict{Float64,Any}()
            getbg(x) = get!(cache,Float64(x)) do
                background(target,x,profile,level;seed=seed)
            end
            function solve_location(evaluator,offset,n_nodes)
                println("LOCATION_START ",target.curve_id," ",level_name," ",evaluator,
                    " offset=",offset," nodes=",n_nodes); flush(stdout)
                function inverse_at(x)
                    bg=getbg(x);st=bg.state
                    value = evaluator=="endpoint_subtracted" ? eta_pair_inverse(st,offset;n_nodes=n_nodes) :
                        MP.meson_mass_equation(:eta,2st.qp.m.u+offset,0.0,0.0,st.qp,st.tp,st.kc)
                    push!(histories,(curve_id=String(target.curve_id),level=level_name,evaluator=evaluator,
                        mass_offset_inv_fm=offset,pair_nodes=n_nodes,xi=Float64(x),inverse_re=real(value),
                        inverse_im=imag(value),light_mass_inv_fm=st.qp.m.u,strange_mass_inv_fm=st.qp.m.s,
                        background_residual=bg.equilibrium.residual_norm))
                    if !(isfinite(real(value)) && isfinite(imag(value)) && abs(imag(value))<1e-8)
                        checkpoint()
                        error("invalid locator: curve=$(target.curve_id) level=$level_name evaluator=$evaluator offset=$offset xi=$x value=$value")
                    end
                    real(value)
                end
                lo,hi=Float64(target.xi_left)-0.01,Float64(target.xi_right)+0.01
                flo,fhi=inverse_at(lo),inverse_at(hi)
                flo*fhi<0 || error("eta threshold lost its original bracket")
                for _ in 1:50
                    hi-lo<=tolerance && break
                    mid=(lo+hi)/2;fm=inverse_at(mid)
                    if flo*fm<=0
                        hi=mid;fhi=fm
                    else
                        lo=mid;flo=fm
                    end
                end
                result=(curve_id=String(target.curve_id),window=String(target.window),
                    T_MeV=Float64(target.T_MeV),muB_MeV=Float64(target.muB_MeV),level=level_name,
                    evaluator=evaluator,mass_offset_inv_fm=offset,pair_nodes=n_nodes,xi_left=lo,
                    xi_right=hi,xi_midpoint=(lo+hi)/2,inverse_left=flo,inverse_right=fhi)
                push!(location_rows,result);checkpoint()
                println("LOCATION_DONE ",target.curve_id," ",level_name," ",evaluator,
                    " offset=",offset," nodes=",n_nodes," xi=",result.xi_midpoint);flush(stdout)
                result
            end
            exact=solve_location("endpoint_subtracted",0.0,128)
            target_locations[level_name]=(xi=exact.xi_midpoint,xi_left=exact.xi_left,xi_right=exact.xi_right,
                bisection_tolerance=tolerance,pair_nodes=128,method="equal-mass endpoint-subtracted threshold inverse")
            if level_name=="base"
                refined=solve_location("endpoint_subtracted",0.0,256)
                target_locations["pair_ref"]=(xi=refined.xi_midpoint,pair_nodes=256,
                    difference=abs(refined.xi_midpoint-exact.xi_midpoint))
                # Finite-offset original-kernel roots are controls, not an
                # assumed sqrt(offset) extrapolation through its endpoint tolerance.
                for eps_mass in Float64.(INPUT_CONFIG.mass_offsets_inv_fm)
                    native=solve_location("original_kernel",-eps_mass,0)
                    target_locations["original_offset_"*string(eps_mass)]=(
                        xi=native.xi_midpoint,mass_offset_inv_fm=-eps_mass,
                        difference=native.xi_midpoint-exact.xi_midpoint)
                end
                bg=getbg(exact.xi_midpoint);st=bg.state
                for offset in (-1e-4,-1e-6,-1e-7,-1e-8,-1e-10,-1e-12,0.0,1e-12,1e-10,1e-8,1e-6,1e-4)
                    old=MP.meson_mass_equation(:eta,2st.qp.m.u+offset,0.0,0.0,st.qp,st.tp,st.kc)
                    stable=eta_pair_inverse(st,offset;n_nodes=256)
                    push!(comparisons,(curve_id=String(target.curve_id),xi=exact.xi_midpoint,
                        mass_offset_inv_fm=offset,original_re=real(old),original_im=imag(old),
                        original_finite=isfinite(real(old))&&isfinite(imag(old)),
                        endpoint_subtracted_re=real(stable),endpoint_subtracted_im=imag(stable)))
                end
                for source_id in (String(target.point_left),String(target.point_right))
                    ref=only(r for r in SOURCE_ROWS if String(r.point_id)==source_id)
                    bg=getbg(Float64(ref.xi))
                    error_m=maximum(abs(bg.equilibrium.masses[k]/Float64(getproperty(ref,(:m_u,:m_d,:m_s)[k]))-1) for k in 1:3)
                    push!(background_controls,(curve_id=String(target.curve_id),point_id=source_id,
                        xi=Float64(ref.xi),mass_max_relative_difference=error_m,
                        Phi_difference=bg.row.Phi-Float64(ref.Phi),Phibar_difference=bg.row.Phibar-Float64(ref.Phibar)))
                    error_m<1e-7 || error("formal endpoint background was not reproduced")
                end
                checks=[getbg(exact.xi_midpoint+d) for d in (-0.001,-0.0001,0.0001,0.001)]
                for direction in ("forward","reverse")
                    previous=nothing
                    for bg in (direction=="forward" ? checks : reverse(checks))
                        row=merge(bg.row,(point_id=String(target.curve_id)*"_"*string(bg.row.xi),))
                        cs,sm,previous=probe_roots(row,bg.state,:eta,direction,previous)
                        append!(pole_candidates,cs);push!(pole_rows,sm)
                    end
                end
                checkpoint()
            end
        end
        push!(locations,(curve_id=String(target.curve_id),window=String(target.window),
            T_MeV=Float64(target.T_MeV),muB_MeV=Float64(target.muB_MeV),levels=target_locations))
        write_json(joinpath(OUTPUT_DIR,"locations.json"),locations)
    end
    points=NamedTuple[]
    for (target,location) in zip(INPUT_CONFIG.targets,locations)
        center=location.levels["base"].xi
        curve=sort(filter(r->String(r.window)==String(target.window),SOURCE_ROWS);by=r->Float64(r.xi))
        at=findfirst(r->String(r.point_id)==String(target.point_left),curve)
        values=[(xi=center+Float64(d),kind="mott_centered",original_point_id="") for d in INPUT_CONFIG.xi_offsets]
        append!(values,[(xi=Float64(curve[i].xi),kind="original_control",original_point_id=String(curve[i].point_id)) for i in at:at+2])
        sort!(values;by=v->abs(v.xi-center))
        for (i,v) in enumerate(values)
            push!(points,(point_id=String(target.curve_id)*"_p"*lpad(string(i),2,'0'),
                curve_id=String(target.curve_id),window=String(target.window),T_MeV=Float64(target.T_MeV),
                muB_MeV=Float64(target.muB_MeV),xi=v.xi,xi_mott=center,delta_xi=v.xi-center,
                kind=v.kind,original_point_id=v.original_point_id))
        end
    end
    CSV.write(joinpath(OUTPUT_DIR,"transport_plan.csv"),points)
    write_json(joinpath(OUTPUT_DIR,"summary.json"),(status="completed_diagnostic",source_root=SOURCE_ROOT,
        calculation_sha=String(INPUT_CONFIG.calculation_sha),input_config_sha256=digest(joinpath(INPUT_DIR,"config.json")),
        analysis_script_sha256=digest(@__FILE__),driver_sha256=digest(joinpath(@__DIR__,"run_transport_eta_refinement.jl")),
        helper_sha256=digest(joinpath(@__DIR__,"run_transport_mott_probe.jl")),
        root_method="equal-mass endpoint subtraction; original finite-offset controls; stable multiseed backgrounds",
        production_kernel_changed=false,locations=locations,transport_plan_points=length(points),
        output_sha256=Dict(f=>digest(joinpath(OUTPUT_DIR,f)) for f in readdir(OUTPUT_DIR))))
end
