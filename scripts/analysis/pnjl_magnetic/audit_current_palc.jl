#!/usr/bin/env julia

# Audit saved current-kernel roots without running a new branch scan.
# Separate stationarity, same-node condensate self consistency, and drift
# relative to the independent Fermi-resolved formula. Compare the actual AD
# Jacobian with centered finite differences at selected states.
using JSON, LinearAlgebra, StaticArrays, ForwardDiff
include(joinpath(@__DIR__, "MagneticReference.jl"))
using .MagneticReference

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main, :Models)
    include(joinpath(ROOT, "src", "models", "Models.jl"))
end
const MT = Models.MagneticThermodynamics
const MI = MT.MagneticIntegrals
const RP = MagneticReference.Parameters()
const H = MagneticReference.HBARC

function same_node_gap(x, mu, T, eB, p_num)
    M = MagneticReference.masses(x[1:3], eB, RP)
    loop = zeros(3)
    for f in 1:3
        q = MagneticReference.CHARGES[f]
        _, vac_loop = MagneticReference.vacuum(M[f], q*eB, RP,
            MagneticReference.Controls())
        therm_loop = ForwardDiff.derivative(M[f]) do mass
            MI.omegat_flavor_landau(mass, mu[f], T, x[4], x[5], q, eB;
                n_max=32, p_num=p_num, pz_max=40.0)
        end
        loop[f] = vac_loop + therm_loop
    end
    return loop - x[1:3]
end

function audit_root(c, r, muB, source, family, index; jacobian=false)
    eB = c["eB_GeV2"] * (1000/H)^2
    T = c["T_MeV"]/H
    mu = SVector{3,Float64}(fill(muB/(3H), 3))
    x = SVector{5,Float64}(r["x_state"])
    conf = MT.MagneticConfig(eB_fm2=eB, n_max=32, p_num=c["p_num"],
        pz_max=40.0)
    model = Models.PNJLMagneticModel(Models.PNJLModel(), conf)
    residual_fn = y -> Models.magnetic_gap_residual(model, y, T, mu;
        p_num=c["p_num"], n_max=32, pz_max=40.0, t_num=4)
    gradient = residual_fn(x)
    Jm = MagneticReference.mass_jacobian(x[1:3], eB, RP)
    gap = same_node_gap(x, mu, T, eB, c["p_num"])
    resolved = MagneticReference.evaluate(x, T, mu, eB, RP,
        MagneticReference.Controls(order=48))
    actual_rho = MT.calculate_magnetic_rho(x, mu, T, conf)
    out = Dict{String,Any}(
        "source"=>source, "family"=>family, "index"=>index,
        "eB_GeV2"=>c["eB_GeV2"], "muB_MeV"=>muB,
        "p_num"=>c["p_num"], "x_state"=>collect(x),
        "stationarity_norm"=>norm(gradient),
        "same_node_condensate_gap_norm"=>norm(gap),
        "chain_rule_error"=>norm(gradient[1:3]-Jm'*gap),
        "mass_jacobian_sigma_min"=>minimum(svdvals(Jm)),
        "resolved_condensate_gap_norm"=>norm(resolved.condensate_gap),
        "resolved_stationarity_norm"=>norm(resolved.stationarity),
        "resolved_rho_delta"=>norm(actual_rho-resolved.rho),
        "rhoB_over_rho0"=>sum(actual_rho)/(3*.16),
        "physical"=>r["physical"],
    )
    if jacobian
        Jad = ForwardDiff.jacobian(y->collect(residual_fn(y)), collect(x))
        errors = Any[]
        for step in (1e-4, 6e-6, 1e-6, 1e-7)
            Jfd = zeros(5,5)
            for j in 1:5
                h = step*max(1.0,abs(x[j]))
                xp,xm = collect(x),collect(x)
                xp[j]+=h
                xm[j]-=h
                Jfd[:,j] = (residual_fn(xp)-residual_fn(xm))/(2h)
            end
            push!(errors,(;step,relative_error=norm(Jfd-Jad)/max(norm(Jad),eps())))
        end
        out["jacobian_checks"] = errors
        out["residual_jacobian_sigma_min"] = minimum(svdvals(Jad))
    end
    return out
end

function main(args=ARGS)
    length(args)>=2 || throw(ArgumentError("output JSON followed by saved current-kernel input JSON files required"))
    records = Any[]
    for file in args[2:end]
        for c in JSON.parsefile(file)["cases"]
            haskey(c,"targets") || continue
            # Every retained target root, regardless of whether PALC and the
            # catalogue agree. These are cheap fixed-state evaluations.
            for t in c["targets"]
                for family in ("palc_roots","catalog_roots")
                    for (i,r) in enumerate(t[family])
                        push!(records,audit_root(c,r,t["muB_MeV"],basename(file),family,i))
                    end
                end
            end
            if c["eB_GeV2"] == 0.4
                pts = c["points"]
                selected = unique(vcat(1, length(pts),
                    argmin([pt["muB_MeV"] for pt in pts]),
                    argmax([pt["muB_MeV"] for pt in pts])))
                for i in selected
                    push!(records,audit_root(c,pts[i],pts[i]["muB_MeV"],basename(file),"branch",i;jacobian=true))
                end
            end
        end
    end
    open(args[1],"w") do io
        JSON.print(io,(;julia_version=string(VERSION),records),2)
    end
    println("audited ",length(records)," fixed states -> ",args[1])
end

if abspath(PROGRAM_FILE)==@__FILE__
    main()
end
