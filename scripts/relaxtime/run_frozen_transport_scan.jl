# Recompute transport on published mean-field rows without a phase scan.
using CSV, JSON3, SHA, Dates, LinearAlgebra, StaticArrays, NLsolve
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
if !isdefined(Main, :Models)
    include(joinpath(ROOT, "src", "models", "Models.jl"))
end
const SP_FROZEN = Models.TransportSamplingProfiles
const RT_FROZEN = Main.RelaxTime
const SPECIES_FROZEN = (:u,:d,:s,:ubar,:dbar,:sbar)
BLAS.set_num_threads(1)

function frozen_options(args)
    allowed = Set(("input","outdir","sampling-profile","point-list","shard-index","shard-count","store-tables"))
    opts = Dict{String,String}()
    iseven(length(args)) || error("arguments require --name value pairs")
    for i in 1:2:length(args)
        startswith(args[i], "--") || error("expected --option")
        key = args[i][3:end]
        key in allowed || error("unknown option: $key")
        haskey(opts,key) && error("duplicate option: $key")
        opts[key] = args[i+1]
    end
    all(k->haskey(opts,k),("input","outdir","sampling-profile")) ||
        error("required: --input CSV --outdir DIR --sampling-profile PROFILE")
    n, i = parse(Int,get(opts,"shard-count","1")), parse(Int,get(opts,"shard-index","1"))
    1 <= i <= n || error("shard-index must be in 1:shard-count")
    return opts
end

function frozen_atomic_csv(path, rows)
    isempty(rows) && return
    CSV.write(path*".writing", rows)
    mv(path*".writing", path;force=true)
end

function frozen_row_id(row, index)
    if hasproperty(row, :probe_id)
        return String(row.probe_id)
    end
    if hasproperty(row, :point_id) && !ismissing(row.point_id)
        return String(row.point_id)
    end
    return "p"*lpad(index,5,'0')
end

function frozen_background(row)
    qp = (m=(u=Float64(row.m_u),d=Float64(row.m_d),s=Float64(row.m_s)),
        μ=(u=Float64(row.muq_fm),d=Float64(row.muq_fm),s=Float64(row.muq_fm)))
    tp = (T=Float64(row.T_fm),Φ=Float64(row.Phi),Φbar=Float64(row.Phibar),ξ=Float64(row.xi))
    n = NamedTuple{SPECIES_FROZEN}(Tuple(Float64(getproperty(row,Symbol("n_",s))) for s in SPECIES_FROZEN))
    tau = NamedTuple{SPECIES_FROZEN}(Tuple(Float64(getproperty(row,Symbol("tau_",s))) for s in SPECIES_FROZEN))
    return qp,tp,n,tau
end

function frozen_bulk(runtime, row, qp, tp)
    target = SVector(qp.m.u,qp.m.d,qp.m.s)
    # Only invert the algebraic mass-condensate map. No gap/phase solve occurs.
    f(x) = Vector(Models.calculate_mass_vec(runtime.model,SVector{3}(x))-target)
    root = nlsolve(f,[-0.9,-0.9,-1.9];ftol=1e-12,xtol=1e-12)
    maximum(abs,f(root.zero)) < 1e-10 || error("mass-condensate inversion failed")
    state = SVector(root.zero...,tp.Φ,tp.Φbar)
    return Models.ThermoDerivatives.bulk_viscosity_coefficients(tp.T,qp.μ.u;
        model=runtime.model,xi=tp.ξ,p_num=runtime.profile.background.p_nodes,
        t_num=runtime.profile.background.cos_nodes,base_state=state,base_masses=target,
        base_mu_vec=SVector(qp.μ.u,qp.μ.d,qp.μ.s),base_state_source=:workflow_equilibrium,
        allow_local_polish=false)
end

function frozen_terminal_kernels(qp,tp,config,bc)
    tc = RT_FROZEN.TransportCoefficients
    J = NamedTuple{SPECIES_FROZEN}(Tuple(tc.electric_conductivity(qp,tp;
        tau=NamedTuple{SPECIES_FROZEN}(Tuple(s===f ? 1.0 : 0.0 for s in SPECIES_FROZEN)),config)/tp.T
        for f in SPECIES_FROZEN))
    Z = NamedTuple{SPECIES_FROZEN}(Tuple(tc.bulk_viscosity_isentropic(qp,tp;
        bulk_coeffs_isentropic=bc,tau=NamedTuple{SPECIES_FROZEN}(Tuple(s===f ? 1.0 : 0.0 for s in SPECIES_FROZEN)),
        config) for f in SPECIES_FROZEN))
    return J,Z
end

function frozen_emit(runtime,row,id,dir;store_tables=false)
    qp,tp,n,old_tau = frozen_background(row)
    config = runtime.transport_config
    bc = frozen_bulk(runtime,row,qp,tp)
    tc = RT_FROZEN.TransportCoefficients
    original_config = SP_FROZEN.build_runtime(SP_FROZEN.load_profile("thermal_weighted_v1")).transport_config
    old_eta = tc.shear_viscosity(qp,tp;tau=old_tau,config=original_config)
    old_sigma = tc.electric_conductivity(qp,tp;tau=old_tau,config=original_config)
    old_zeta = tc.bulk_viscosity_isentropic(qp,tp;bulk_coeffs_isentropic=bc,tau=old_tau,config=original_config)
    replay = (eta=old_eta/Float64(row.eta)-1,sigma=old_sigma/Float64(row.sigma)-1,
        zeta=old_zeta/Float64(row.zeta)-1)
    maximum(abs,values(replay)) < 1e-6 || error("fixed-background terminal replay drift: $replay")
    J,Z = frozen_terminal_kernels(qp,tp,config,bc)
    records = NamedTuple[]
    previous = Dict{Symbol,Float64}()
    rate_path = joinpath(dir,"rates.csv")
    if isfile(rate_path)
        for r in CSV.File(rate_path)
            push!(records,NamedTuple(r))
            previous[Symbol(r.channel)] = Float64(r.rate)
        end
    end
    function on_channel(r,cache)
        push!(records,r)
        frozen_atomic_csv(rate_path,records)
        if store_tables
            table = cache.table
            frozen_atomic_csv(joinpath(dir,r.channel*"_sigma_table.csv"),
                [(s_inv_fm2=table.nodes[j],sigma_fm2=table.sigma[j]) for j in eachindex(table.nodes)])
        end
        println("RATE ",id," ",r.channel," ",r.rate," ",r.seconds)
        flush(stdout)
    end
    result = SP_FROZEN.solve_frozen_transport(runtime,qp,tp,n,Float64(row.s_fm3inv);
        bulk_coeffs=bc,on_channel,existing_rates=previous)
    tau,inv_tau,tr = result.tau,result.tau_inv,result.transport
    check_sigma = sum(tau[s]*J[s] for s in SPECIES_FROZEN)
    check_zeta = sum(tau[s]*Z[s] for s in SPECIES_FROZEN)
    isapprox(check_sigma,tr.sigma/tp.T;rtol=1e-11) || error("conductivity kernel decomposition failed")
    isapprox(check_zeta,tr.zeta;rtol=1e-11) || error("bulk kernel decomposition failed")
    data = Dict(Symbol(k)=>getproperty(row,Symbol(k)) for k in propertynames(row))
    for s in SPECIES_FROZEN
        data[Symbol("tau_",s)] = tau[s]
        data[Symbol("tauinv_",s)] = inv_tau[s]
        data[Symbol("J_",s)] = J[s]
        data[Symbol("Z_",s)] = Z[s]
    end
    data[:eta],data[:sigma],data[:zeta] = tr.eta,tr.sigma,tr.zeta
    data[:eta_over_s] = tr.eta/result.entropy
    data[:zeta_over_s] = tr.zeta/result.entropy
    data[:sigma_over_T] = tr.sigma/tp.T
    data[:sigma_over_T_over_eta_over_s] = data[:sigma_over_T]/data[:eta_over_s]
    data[:zeta_over_s_over_eta_over_s] = data[:zeta_over_s]/data[:eta_over_s]
    data[:sampling_profile] = runtime.profile.id
    data[:a_prescription] = String(runtime.profile.physics.a_method)
    data[:b0_direction] = String(runtime.profile.physics.b0_direction)
    data[:relative_velocity] = String(runtime.profile.sampling.relative_velocity)
    state = result.loop_state
    data[:soft_pion_residual] = 1 + Main.Constants_PNJL.N_color/(2pi^2)*state.K_coeffs.K123_plus*
        (state.quark_params.A.u+state.quark_params.A.d) - Main.Constants_PNJL.m_ud0_inv_fm/qp.m.u
    data[:su_sd_fraction] = 2n.u*result.rates[:us_to_us]/inv_tau.s
    data[:dubar_inverse_contribution] = n.u*result.rates[:dubar_to_dubar]
    data[:source_tau_u] = Float64(row.tau_u)
    data[:source_tau_ubar] = Float64(row.tau_ubar)
    data[:source_eta_over_s] = Float64(row.eta_over_s)
    data[:source_sigma_over_T] = Float64(row.sigma_over_T)
    data[:bulk_replay_error] = replay.zeta
    data[:terminal_replay_error] = max(abs(replay.eta),abs(replay.sigma))
    # Reevaluate diagnostics that depend on the new relaxation times.
    data[:quality_metric] = max(tau.u/tau.ubar,tau.ubar/tau.u)
    data[:quality_flag] = data[:quality_metric] > 10.0
    data[:quality_reason] = data[:quality_flag] ? "tau_u_ubar_ratio_high" : "ok"
    names = Tuple(sort!(collect(keys(data));by=String))
    output = NamedTuple{names}(Tuple(data[k] for k in names))
    frozen_atomic_csv(joinpath(dir,"transport.csv"),[output])
    println("POINT_DONE ",id," ",data[:eta_over_s]," ",data[:zeta_over_s]," ",data[:sigma_over_T])
    flush(stdout)
    return output
end

function frozen_main(args)
    opts = frozen_options(args)
    input,outdir = abspath(opts["input"]),abspath(opts["outdir"])
    profile = SP_FROZEN.load_profile(opts["sampling-profile"])
    runtime = SP_FROZEN.build_runtime(profile)
    rows = collect(CSV.File(input))
    n,i = parse(Int,get(opts,"shard-count","1")),parse(Int,get(opts,"shard-index","1"))
    wanted = haskey(opts,"point-list") ? Set(split(opts["point-list"],',')) : nothing
    selected = [(index,row,frozen_row_id(row,index)) for (index,row) in enumerate(rows)
        if (index-1)%n==i-1 && (wanted===nothing || frozen_row_id(row,index) in wanted)]
    isempty(selected) && error("no rows selected")
    all(x->occursin(r"^[A-Za-z0-9_.-]+$",x[3]),selected) || error("unsafe point ID")
    mkpath(outdir)
    source_files = String[]
    for subtree in ("src","config/physics","config/models/pnjl")
        for (dir,_,files) in walkdir(joinpath(ROOT,subtree)), name in files
            endswith(name,".jl") || endswith(name,".toml") || continue
            push!(source_files,joinpath(dir,name))
        end
    end
    push!(source_files,@__FILE__)
    source_digest = bytes2hex(sha256(join([replace(relpath(file,ROOT),'\\'=>'/')*":"*
        bytes2hex(sha256(read(file))) for file in sort!(source_files)],"\n")))
    identity = (schema="frozen_transport_v1",input_sha256=bytes2hex(sha256(read(input))),
        profile_sha256=profile.source_sha256,profile_id=profile.id,
        source_tree_sha256=get(ENV,"PNJL_SOURCE_TREE_SHA256","unrecorded"),
        numerical_source_sha256=source_digest,
        points=[x[3] for x in selected],shard_index=i,shard_count=n)
    contract = joinpath(outdir,"contract.json")
    if isfile(contract)
        JSON3.write(JSON3.read(read(contract,String))) == JSON3.write(identity) || error("resume identity mismatch")
    else
        write(contract,JSON3.write(identity))
    end
    write(joinpath(outdir,"effective_config.json"),JSON3.write(SP_FROZEN.effective_config(profile)))
    println("HOST ",get(ENV,"COMPUTERNAME",get(ENV,"HOSTNAME","?"))," JULIA ",VERSION,
        " PROFILE ",profile.id," POINTS ",length(selected))
    failures = String[]
    results = NamedTuple[]
    for (_,row,id) in selected
        dir = joinpath(outdir,id)
        mkpath(dir)
        try
            if isfile(joinpath(dir,"transport.csv"))
                push!(results,NamedTuple(only(CSV.File(joinpath(dir,"transport.csv")))))
            else
                measured = @timed frozen_emit(runtime,row,id,dir;
                    store_tables=get(opts,"store-tables","false")=="true")
                push!(results,measured.value)
                write(joinpath(dir,"timing.json"),JSON3.write((seconds=measured.time,
                    allocated_bytes=measured.bytes,gc_seconds=measured.gctime)))
            end
            frozen_atomic_csv(joinpath(outdir,"transport.csv"),results)
        catch e
            e isa InterruptException && rethrow()
            msg = sprint(showerror,e,catch_backtrace())
            write(joinpath(dir,"error.txt"),msg)
            push!(failures,id)
            println(stderr,"POINT_FAILED ",id," ",msg)
            flush(stderr)
        end
    end
    write(joinpath(outdir,"summary.json"),JSON3.write((completed=length(results),failed=failures)))
    isempty(failures) || error("failed frozen points: $(join(failures,','))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    frozen_main(ARGS)
end
