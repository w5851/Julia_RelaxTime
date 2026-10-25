#!/usr/bin/env julia
using CSV
if !isdefined(Main,:ThresholdUpgradeCommon)
    include(joinpath(@__DIR__,"ThresholdUpgradeCommon.jl"))
end
using .ThresholdUpgradeCommon

function density_validation(inputdir,outdir)
    rows=collect(CSV.File(joinpath(inputdir,"review_points.csv")))
    isempty(rows) && error("empty density validation selection")
    configs=((name="finite10_p64_c8",nodes=64,angles=8,scale=10.0,semi=false),
        (name="finite10_p128_c32",nodes=128,angles=32,scale=10.0,semi=false),
        (name="finite20_p96_c32",nodes=96,angles=32,scale=20.0,semi=false),
        (name="semi32_scale10_c20",nodes=32,angles=20,scale=10.0,semi=true),
        (name="semi32_scale5_c16",nodes=32,angles=16,scale=5.0,semi=true),
        (name="semi64_scale5_c32",nodes=64,angles=32,scale=5.0,semi=true),
        (name="semi128_scale10_c64",nodes=128,angles=64,scale=10.0,semi=true))
    start_run(outdir,inputdir,@__FILE__,(purpose="external tau density and thermal-rate normalization",
        configurations=configs,relative_target=0.01,
        reference="RS change of variable: n(xi)=n(0)/sqrt(1+xi), semi-infinite radial p256 scale5",
        source="fixed retained masses, chemical potentials and Polyakov loop"))
    results=NamedTuple[];failures=NamedTuple[]
    rp,rw=AS._cached_semi_infinite_momentum_grid(256,5.0)
    species=(:u,:d,:s,:ubar,:dbar,:sbar)
    for row in rows
        try
            qp=(m=(u=Float64(row.m_u),d=Float64(row.m_d),s=Float64(row.m_s)),
                μ=(u=Float64(row.muq_fm),d=Float64(row.muq_fm),s=Float64(row.muq_fm)))
            T,Phi,Phibar,xi=Float64(row.T_fm),Float64(row.Phi),Float64(row.Phibar),Float64(row.xi)
            for s in species
                m,mu=AS.get_mass(s,qp),AS.get_mu(s,qp)
                reference=AS.number_density(s,m,mu,T,Phi,Phibar,0.0;
                    p_grid=rp,p_w=rw,angle_nodes=1)/sqrt(1+xi)
                for config in configs
                    p,w=config.semi ? AS._cached_semi_infinite_momentum_grid(config.nodes,config.scale) :
                        Main.RelaxTime.GaussLegendre.gauleg(0.0,config.scale,config.nodes)
                    n=AS.number_density(s,m,mu,T,Phi,Phibar,xi;p_grid=p,p_w=w,angle_nodes=config.angles)
                    isfinite(n) && n>0 && isfinite(reference) && reference>0 || error("invalid density")
                    push!(results,(point_id=String(row.point_id),muB_MeV=Float64(row.muB_MeV),
                        alpha_T=Float64(row.alpha_T),xi=xi,species=String(s),configuration=config.name,
                        density=n,reference=reference,relative_difference=n/reference-1,
                        retained_density=Float64(getproperty(row,Symbol(:n_,s)))))
                end
            end
        catch err
            err isa InterruptException && rethrow()
            push!(failures,(point_id=String(row.point_id),error=sprint(showerror,err,catch_backtrace())))
        end
    end
    write_rows(joinpath(outdir,"densities.csv"),results)
    write_rows(joinpath(outdir,"failures.csv"),failures)
    checks=filter(r->r.configuration=="semi128_scale10_c64",results)
    finish_run(outdir;points=length(rows),comparisons=length(results),failures=length(failures),
        rs_identity_checked=length(checks)==6length(rows) && isempty(failures) &&
            all(r->abs(r.relative_difference)<1e-9,checks))
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)==2 || error("Usage: run_transport_density_validation.jl INPUT_DIR NEW_OUTPUT_DIR")
    density_validation(abspath(ARGS[1]),abspath(ARGS[2]))
end
