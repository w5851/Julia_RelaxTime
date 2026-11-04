#!/usr/bin/env julia
# Reinsert only the four investigated channel rates into the frozen formal rows.
# Execute on quote-remote: SOURCE_ROOT INPUT_DIR OUTPUT_DIR RUN_DIR_NAMES [POINT_IDS CHANNELS OUTER_CHANNELS]
using CSV, JSON3, SHA, Dates, LinearAlgebra

length(ARGS) in (4,7) || error("expected SOURCE_ROOT INPUT_DIR OUTPUT_DIR RUN_DIR_NAMES [POINT_IDS CHANNELS OUTER_CHANNELS]")
const ROOT=abspath(ARGS[1])
const INPUT=abspath(ARGS[2])
const OUTPUT=abspath(ARGS[3])
const CASE=dirname(OUTPUT)
isdir(OUTPUT) && !isempty(readdir(OUTPUT)) && error("refusing to overwrite retained analysis")
mkpath(OUTPUT)
include(joinpath(ROOT,"src","models","Models.jl"))
const AS=Main.RelaxTime.AverageScatteringRate
const RT=Main.RelaxTime.RelaxationTime
const TC=Main.RelaxTime.TransportCoefficients
const SPECIES=(:u,:d,:s,:ubar,:dbar,:sbar)
const POINTS=length(ARGS)==7 ? Tuple(split(ARGS[5],',')) : ("p004","p005","p006","p014","p015","p016")
const CHANNELS=length(ARGS)==7 ? Tuple(Symbol.(split(ARGS[6],','))) : (:udbar_to_udbar,:dubar_to_dubar,:uubar_to_uubar,:uubar_to_ddbar)
const OUTER_CHANNELS=length(ARGS)==7 ? Tuple(Symbol.(split(ARGS[7],','))) : (:udbar_to_udbar,:dubar_to_dubar)
length(unique(POINTS))==length(POINTS)>0 || error("points must be nonempty and unique")
length(unique(CHANNELS))==length(CHANNELS)>0 || error("channels must be nonempty and unique")
all(c->c in RT.REQUIRED_RATE_KEYS_FOR_TAU,CHANNELS) || error("unsupported selected channel")
all(c->c in CHANNELS,OUTER_CHANNELS) || error("outer channels must belong to the s study")
BLAS.set_num_threads(1)

digest(path)=bytes2hex(sha256(read(path)))
function write_json(path,data)
    open(path,"w") do io
        JSON3.pretty(io,data)
        println(io)
    end
end

function main()
    profile=Models.TransportSamplingProfiles.load_profile("thermal_weighted_v1")
    cfg=profile.terminal_integral
    p,w=AS._cached_semi_infinite_momentum_grid(cfg.p_nodes,cfg.scale_inv_fm)
    c,cw=AS._cached_interval_grid(-1.0,1.0,cfg.cos_nodes)
    integration=TC.TransportIntegrationConfig(p_grid=p,p_w=w,cos_grid=c,cos_w=cw)
    source=collect(CSV.File(joinpath(INPUT,"source_rows.csv")))
    source_rates=collect(CSV.File(joinpath(INPUT,"source_rates.csv")))
    runs=split(ARGS[4],',')
    all(r->occursin(r"^[a-zA-Z0-9_-]+$",r),runs) || error("invalid retained run names")
    all_channels=reduce(vcat,[collect(CSV.File(joinpath(CASE,r,"channels.csv"))) for r in runs])
    corrections=filter(r->hasproperty(r,:reference_rate),all_channels)
    outer=filter(r->hasproperty(r,:joint1_rate),all_channels)
    length(corrections)+length(outer)==length(all_channels) || error("unknown channel result schema")
    length(corrections)==length(POINTS)*length(CHANNELS) || error("incomplete channel study")
    index=Dict{Tuple{String,Symbol},Any}()
    for row in corrections
        key=(String(row.point_id),Symbol(row.channel))
        haskey(index,key) && error("duplicate channel result $key")
        index[key]=row
    end
    all(haskey(index,(p,c)) for p in POINTS for c in CHANNELS) || error("missing selected point/channel")
    outer_index=Dict{Tuple{String,Symbol},Any}()
    if !isempty(outer)
        length(outer)==length(POINTS)*length(OUTER_CHANNELS) || error("incomplete outer channel study")
        for row in outer
            key=(String(row.point_id),Symbol(row.channel))
            haskey(outer_index,key) && error("duplicate outer channel result $key")
            outer_index[key]=row
        end
        all(haskey(outer_index,(p,c)) for p in POINTS for c in OUTER_CHANNELS) || error("missing outer point/channel")
    end
    observable_rows=NamedTuple[]
    tau_rows=NamedTuple[]
    replay_rows=NamedTuple[]
    configurations=Tuple{String,Union{Nothing,Symbol}}[("original",nothing),("base",:base_rate),("tight1",:tight1_rate),
        ("tight2",:tight2_rate),("reference",:reference_rate)]
    isempty(outer) || append!(configurations,[("outer_joint1",:joint1_rate),("outer_joint2",:joint2_rate)])
    for point in POINTS
        row=only(r for r in source if String(r.point_id)==point)
        qp=(m=(u=Float64(row.m_u),d=Float64(row.m_d),s=Float64(row.m_s)),
            μ=(u=Float64(row.muq_fm),d=Float64(row.muq_fm),s=Float64(row.muq_fm)))
        tp=(T=Float64(row.T_fm),Φ=Float64(row.Phi),Φbar=Float64(row.Phibar),ξ=Float64(row.xi))
        densities=(;(s=>Float64(getproperty(row,Symbol("n_",s))) for s in SPECIES)...)
        original=Dict(Symbol(r.channel)=>Float64(r.rate) for r in source_rates if String(r.point_id)==point)
        Set(keys(original))==Set(RT.REQUIRED_RATE_KEYS_FOR_TAU) || error("formal channel set differs")
        for (name,column) in configurations
            rates=copy(original)
            if column!==nothing
                for channel in CHANNELS
                    source_column=startswith(name,"outer_") ? :reference_rate : column
                    rates[channel]=Float64(getproperty(index[(point,channel)],source_column))
                end
                if startswith(name,"outer_")
                    for channel in OUTER_CHANNELS
                        rates[channel]=Float64(getproperty(outer_index[(point,channel)],column))
                    end
                end
            end
            inverse=RT.relaxation_rates(densities,rates)
            all(x->isfinite(x) && x>0,inverse) || error("invalid reconstructed inverse tau")
            tau=(;(s=>1/getproperty(inverse,s) for s in SPECIES)...)
            eta=TC.shear_viscosity(qp,tp;tau=tau,config=integration)
            sigma=TC.electric_conductivity(qp,tp;tau=tau,config=integration)
            eta_s=eta/Float64(row.s_fm3inv)
            sigma_T=sigma/tp.T
            all(x->isfinite(x) && x>0,(eta,sigma,eta_s,sigma_T)) || error("invalid terminal observable")
            if name=="original"
                tau_error=maximum(abs(getproperty(tau,s)/Float64(getproperty(row,Symbol("tau_",s)))-1) for s in SPECIES)
                eta_error=eta/Float64(row.eta)-1
                sigma_error=sigma/Float64(row.sigma)-1
                maximum(abs.((tau_error,eta_error,sigma_error)))<=1e-10 || error("formal terminal replay failed at $point")
                push!(replay_rows,(point_id=point,tau_max_relative=tau_error,
                    eta_relative=eta_error,sigma_relative=sigma_error,passed=true))
            end
            info=(point_id=point,window=String(row.window),T_MeV=Float64(row.T_MeV),
                muB_MeV=Float64(row.muB_MeV),xi=Float64(row.xi),configuration=name)
            push!(observable_rows,merge(info,(eta=eta,sigma=sigma,eta_over_s=eta_s,
                sigma_over_T=sigma_T,sigma_over_T_over_eta_over_s=sigma_T/eta_s,
                eta_relative_to_formal=eta/Float64(row.eta)-1,
                sigma_relative_to_formal=sigma/Float64(row.sigma)-1)))
            for species in SPECIES
                formal=Float64(getproperty(row,Symbol("tau_",species)))
                push!(tau_rows,merge(info,(species=String(species),tau_fm=getproperty(tau,species),
                    inverse_tau_inv_fm=getproperty(inverse,species),
                    relative_to_formal=getproperty(tau,species)/formal-1)))
            end
        end
    end
    shape_rows=NamedTuple[]
    for window in unique(r.window for r in observable_rows), (name,_) in configurations,
        observable in (:sigma_over_T,:eta_over_s,:sigma_over_T_over_eta_over_s)
        scan=sort(filter(r->r.window==window && r.configuration==name,observable_rows);by=r->r.xi)
        length(scan)==3 || error("local shape needs three retained points")
        left,center,right=scan
        values=[getproperty(r,observable) for r in scan]
        x=Float64[r.xi for r in scan]
        chord=values[1]+(x[2]-x[1])/(x[3]-x[1])*(values[3]-values[1])
        slopes=((values[2]-values[1])/(x[2]-x[1]),(values[3]-values[2])/(x[3]-x[2]))
        curvature=2*(slopes[2]-slopes[1])/(x[3]-x[1])
        push!(shape_rows,(window=window,configuration=name,observable=String(observable),
            left_xi=x[1],center_xi=x[2],right_xi=x[3],left_value=values[1],center_value=values[2],
            right_value=values[3],center_above_neighbor_chord=values[2]-chord,
            relative_center_above_neighbor_chord=values[2]/chord-1,
            second_divided_difference=curvature,
            center_is_local_max=values[2]>max(values[1],values[3]),
            sampled_max_xi=x[argmax(values)]))
    end
    CSV.write(joinpath(OUTPUT,"formal_replay.csv"),replay_rows)
    CSV.write(joinpath(OUTPUT,"observables.csv"),observable_rows)
    CSV.write(joinpath(OUTPUT,"relaxation_times.csv"),tau_rows)
    CSV.write(joinpath(OUTPUT,"local_shapes.csv"),shape_rows)
    reference=filter(r->r.configuration=="reference",observable_rows)
    files=Dict(name=>digest(joinpath(OUTPUT,name)) for name in readdir(OUTPUT))
    write_json(joinpath(OUTPUT,"summary.json"),(status="completed_partial_channel_s_resolution_analysis",
        source_sha="2e5c411c716d5cbc524c3fc6cfb335fe05b39444",script_sha256=digest(@__FILE__),
        julia_version=string(VERSION),threads=Threads.nthreads(),blas_threads=BLAS.get_num_threads(),
        finished_utc=string(now(UTC)),points=POINTS,replaced_channels=CHANNELS,
        outer_replaced_channels=isempty(outer) ? Symbol[] : collect(OUTER_CHANNELS),
        all_outer_converged=isempty(outer) ? nothing : all(r.outer_converged for r in outer),
        unchanged_channels=length(RT.REQUIRED_RATE_KEYS_FOR_TAU)-length(CHANNELS),
        terminal_profile=Models.TransportSamplingProfiles.effective_config(profile),
        max_sigma_relative_to_formal=maximum(abs(r.sigma_relative_to_formal) for r in reference),
        max_eta_relative_to_formal=maximum(abs(r.eta_relative_to_formal) for r in reference),
        input_hashes=Dict(r=>digest(joinpath(CASE,r,"channels.csv")) for r in runs),files=files,
        scope="Only the explicitly listed channel s-interpolation and optional outer-grid changes. Fixed original densities, background, t quadrature and terminal quadrature. No zeta derivative replay."))
    println("completed transport propagation for ",length(POINTS)," points")
end

main()
