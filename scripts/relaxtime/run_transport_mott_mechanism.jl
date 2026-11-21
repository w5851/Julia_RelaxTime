# Energy-resolved dominant channel and species-weight decomposition.
include(joinpath(@__DIR__,"run_frozen_transport_scan.jl"))
include(joinpath(@__DIR__,"run_frozen_meson_spectra.jl"))

function mechanism_main(args)
    length(args)==3 || error("INPUT OUTPUT PROFILE")
    input,outdir,profile_name=args
    runtime=SP_FROZEN.build_runtime(SP_FROZEN.load_profile(profile_name))
    rows=collect(CSV.File(input))
    mkpath(outdir)
    spectrum,areas,poles=NamedTuple[],NamedTuple[],NamedTuple[]
    as=RT_FROZEN.AverageScatteringRate
    ts=as.ThermalWeightedSampling
    hc=Main.Constants_PNJL.ħc_MeV_fm
    for (index,row) in enumerate(rows)
        id=frozen_row_id(row,index)
        qp,tp,n,old_tau=frozen_background(row)
        state=SP_FROZEN.fixed_point_state(runtime.profile,qp,tp)
        for meson in (:pi,:K)
            pole=spectral_pole(meson,state,nothing)
            threshold=RT_FROZEN.MottTransition.mott_threshold_mass(meson,qp)
            push!(poles,(point_id=id,T_MeV=row.T_MeV,xi=row.xi,meson=String(meson),
                mass_MeV=pole.mass*hc,width_MeV=pole.width*hc,threshold_MeV=threshold*hc,
                residual=pole.residual,status=pole.status))
        end
        problem=as._prepare_thermal_cache_problem(:dubar_to_dubar,state.quark_params,state.thermo_params,
            state.K_coeffs;runtime.tau_kwargs...)
        cache=as._build_thermal_weighted_cache(:dubar_to_dubar,problem)
        measure=ts.thermal_measure(problem.plan,problem.context.bounds;
            relative_velocity=:medium_moller,max_events=runtime.profile.sampling.max_events)
        threshold=2qp.m.u
        for excess in range(0.01,500.0;length=401)
            s=(threshold+excess/hc)^2
            s<last(problem.context.bounds) || continue
            sigma=as.get_sigma(cache,s;direct=true)
            push!(spectrum,(point_id=id,T_MeV=row.T_MeV,xi=row.xi,
                excess_energy_MeV=excess,sqrt_s_MeV=sqrt(s)*hc,sigma_fm2=sigma))
        end
        bounds=(0.0,25.0,75.0,150.0,Inf)
        contributions=zeros(4)
        for (s,w) in zip(measure.nodes,measure.weights)
            energy=(sqrt(s)-threshold)*hc
            bin=searchsortedlast(collect(bounds),energy)
            1<=bin<=4 || error("unexpected collision energy")
            contributions[bin]+=w*cache.table(s)
        end
        isapprox(sum(contributions),cache.rate;rtol=1e-10) || error("energy partition mismatch")
        for j in 1:4
            push!(areas,(point_id=id,T_MeV=row.T_MeV,xi=row.xi,
                lower_excess_MeV=bounds[j],upper_excess_MeV=bounds[j+1],
                rate_contribution=contributions[j],total_rate=cache.rate,
                inverse_ubar_contribution=n.u*contributions[j]))
        end
        frozen_atomic_csv(joinpath(outdir,"spectrum.csv"),spectrum)
        frozen_atomic_csv(joinpath(outdir,"thermal_areas.csv"),areas)
        frozen_atomic_csv(joinpath(outdir,"poles.csv"),poles)
        RT_FROZEN.PolarizationCache.reset_cache!()
        GC.gc()
        println("MECHANISM_DONE ",id," ",cache.rate);flush(stdout)
    end
    write(joinpath(outdir,"manifest.json"),JSON3.write((input_sha256=bytes2hex(sha256(read(input))),
        profile=SP_FROZEN.effective_config(runtime.profile),rows=length(rows),
        energy_partition_MeV=(0,25,75,150),source_sha=get(ENV,"PNJL_SOURCE_TREE_SHA256","unrecorded"))))
end

if abspath(PROGRAM_FILE)==@__FILE__
    mechanism_main(ARGS)
end
