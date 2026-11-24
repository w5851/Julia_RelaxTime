#!/usr/bin/env julia
# Replay the paired direction and blocking checks quoted in Sec. II C.
# julia --project=SOURCE replay_ward_moller_rate_controls.jl SOURCE ARCHIVE GRID REPLAY OUT
# ARCHIVE contains the original direction outputs and the angular blocking runs.
include(joinpath(@__DIR__, "replay_ward_moller_rate_assembly.jl"))

function wc_rates(file)
    rows = CSV.File(file)
    rates = Dict(Symbol(r.channel) => Float64(r.rate) for r in rows)
    length(rates) == length(rows) || error("duplicate channel in $file")
    wr_validate_rates(rates)
end

wc_change(a,b) = 100*(b/a-1)

function wc_main(args)
    length(args) == 5 || error("expected SOURCE ARCHIVE GRID REPLAY OUT")
    _, archive, grid, replay, out = abspath.(args)
    ispath(out) && error("use a new output directory")
    BLAS.set_num_threads(1)
    config = MM_SP.build_runtime(mm_profile()).transport_config
    inputs = Dict{String,String}()
    tracked(file) = (inputs[replace(relpath(file, archive), '\\' => '/')] = mm_digest(file); file)
    function direction(folder)
        dir = joinpath(archive,folder)
        grid_file = tracked(joinpath(dir,"transport.csv"))
        records = CSV.File(grid_file)
        rates = Dict(String(row.point_id) => wc_rates(tracked(joinpath(dir,String(row.probe_id),"rates.csv")))
            for row in records)
        first(wr_grid(grid_file,rates,config))
    end
    aligned = direction("direction_aligned_v2")
    sphere = Dict(r.probe_id=>r for r in direction("direction_sphere_v2"))
    observables = (Symbol.("tau_", WR_SPECIES)..., :eta_over_s, :zeta_over_s, :sigma_over_T)
    dir_rows = [(probe_id=r.probe_id,T_MeV=wr_float(r,:T_MeV),xi=wr_float(r,:xi),
        observable=String(k),aligned=Float64(r[k] isa AbstractString ? parse(Float64,r[k]) : r[k]),
        sphere_average=Float64(sphere[r.probe_id][k] isa AbstractString ?
            parse(Float64,sphere[r.probe_id][k]) : sphere[r.probe_id][k]))
        for r in aligned for k in observables]
    dir_rows = [merge(r,(change_percent=wc_change(r.aligned,r.sphere_average),)) for r in dir_rows]

    original = Dict(String(r.point_id)=>r for r in CSV.File(joinpath(grid,"mechanism_transport.csv")))
    corrected = Dict(String(r.point_id)=>r for r in CSV.File(joinpath(replay,"mechanism_transport.csv")))
    maps = (("final_blocking_low_v1","dense126_-30","a0654"),
        ("final_blocking_low_v1","dense126_-40","a0653"),
        ("final_blocking_mid_v1","dense138_170","a0776"),
        ("final_blocking_mid_v1","dense138_180","a0777"))
    block = map(maps) do (folder,probe,id)
        avg = only(CSV.File(tracked(joinpath(archive,folder,probe,"transport.csv"));types=String))
        base, ref = original[id], corrected[id]
        rates_file = tracked(joinpath(archive,folder,probe,"rates.csv"))
        mm_digest(rates_file) == mm_digest(joinpath(grid,"blocking_average_rates_$id.csv")) ||
            error("blocking rate provenance mismatch at $id")
        rates = wc_rates(rates_file)
        dens = NamedTuple{WR_SPECIES}(Tuple(getproperty(base,Symbol("n_",s)) for s in WR_SPECIES))
        tau0, J = wr_species(avg,"tau_"), wr_species(avg,"J_")
        _, tau, _ = wr_assemble(dens,rates,tau0,id;rtol=1e-10)
        @assert sum(tau0[s]*J[s] for s in WR_SPECIES) ≈ wr_float(avg,:sigma_over_T)
        sigma_T = sum(tau[s]*J[s] for s in WR_SPECIES)
        qp = (m=(u=base.m_u,d=base.m_d,s=base.m_s),μ=(u=base.muq_fm,d=base.muq_fm,s=base.muq_fm))
        tp = (T=base.T_fm,Φ=base.Phi,Φbar=base.Phibar,ξ=base.xi)
        eta(t) = WR_TC.shear_viscosity(qp,tp; tau=t,config)/base.s_fm3inv
        @assert wr_rel(eta(tau0),wr_float(avg,:eta_over_s)) < 1e-10
        (point_id=id,T_MeV=base.T_MeV,xi=base.xi,
            eta_change_percent=wc_change(ref.eta_over_s,eta(tau)),
            sigma_change_percent=wc_change(ref.sigma_over_T,sigma_T),
            tau_ubar_change_percent=wc_change(ref.tau_ubar,tau.ubar))
    end
    mkpath(out)
    CSV.write(joinpath(out,"direction_comparison.csv"),dir_rows)
    CSV.write(joinpath(out,"final_blocking_comparison.csv"),collect(block))
    mm_json(joinpath(out,"rate_controls.json"),(
        host=get(ENV,"COMPUTERNAME",Sys.MACHINE),julia=string(VERSION),
        driver_sha256=mm_digest(@__FILE__),inputs=inputs,
        files=Dict(basename(p)=>mm_digest(p) for p in readdir(out;join=true))))
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && wc_main(ARGS)
