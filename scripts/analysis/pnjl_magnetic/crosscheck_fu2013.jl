#!/usr/bin/env julia
#= Fu, PRD88 014009 (2013), arXiv:1306.5804 Eqs.5--11, Sec.III.
Independent smooth-regulated three-flavor PNJL with polynomial Polyakov U.
Target: |dMu/dT| peak 224 -> 240 MeV at eB=0 -> 20 m_pi^2, mu_f=0.
The paper's fitted m_pi=135 MeV supplies the explicit field conversion.
=#
include(joinpath(@__DIR__,"crosscheck_literature.jl"))

function polynomial_polyakov(Phi,Phibar,T,p)
    b2=6.75-1.95p.T0/T+2.625(p.T0/T)^2-7.44(p.T0/T)^3
    U=T^4*(-b2*Phi*Phibar/2-.75*(Phi^3+Phibar^3)/6+7.5*(Phi*Phibar)^2/4)
    du=T^4*(-b2*Phibar/2-.75Phi^2/2+7.5Phi*Phibar^2/2)
    db=T^4*(-b2*Phi/2-.75Phibar^2/2+7.5Phibar*Phi^2/2)
    U,du,db
end

function evaluate_fu(x,T,eB,p,ctl)
    phi=x[1:3];Phi=x[4];Phibar=x[5]
    M=MR.masses(phi,eB,p);G=MR.coupling(eB,p)
    poly,du,db=polynomial_polyakov(Phi,Phibar,T,p)
    loop=zeros(3);vac=0.;therm=0.
    for f in 1:3
        if eB>0
            qB=MR.CHARGES[f]*eB
            vf,lf=MR.vacuum(M[f],qB,p,ctl);vac+=vf;loop[f]=lf
            for n in 0:ctl.n_max
                v=MR.thermal_level(M[f],0.,T,Phi,Phibar,qB,n,ctl,p.colors)
                v.*=qB/(2pi^2)*(n==0 ? 1 : 2)
                therm+=v[1];loop[f]+=v[2];du+=v[4];db+=v[5]
            end
        else
            edges=sort!(unique([0.,.8p.cutoff,p.cutoff,1.2p.cutoff,2p.cutoff,4p.cutoff,ctl.pz_max]))
            v=MR.integrate(edges,ctl.order) do q
                E=sqrt(q^2+M[f]^2);fc=1/(1+(q/p.cutoff)^(2ctl.cutoff_N))
                lp,fp,dp,dbp=MR.occupation(-E/T,Phi,Phibar)
                lm,fm,dbm,dm=MR.occupation(-E/T,Phibar,Phi)
                [-p.colors*fc*E, -p.colors*fc*M[f]/E,
                    -T*(lp+lm),p.colors*M[f]/E*(fp+fm),-T*(dp+dm),-T*(dbp+dbm)] .* (q^2/pi^2)
            end
            vac+=v[1];loop[f]+=v[2]+v[4];therm+=v[3];du+=v[5];db+=v[6]
        end
    end
    gap=loop-phi
    chi=2G*sum(abs2,phi)-4p.K*prod(phi)
    (;omega=vac+therm+poly+chi,masses=M,residual=vcat(gap,du,db))
end

function solve_fu(T,eB,p,ctl,seed)
    function f(y)
        x=vcat(y[1:3],y[4],y[4])
        ev=evaluate_fu(x,T,eB,p,ctl)
        ev.residual[1:4]
    end
    r=solve_fd(f,seed;h=1e-5,ftol=2e-12)
    r.residual<1e-8 || error("Fu equilibrium failed $(r.residual)")
    ev=evaluate_fu(vcat(r.x[1:3],r.x[4],r.x[4]),T,eB,p,ctl)
    (;x=r.x,masses=ev.masses,residual=r.residual)
end

function fu_peak(eB,p,ctl,guess;delta_MeV=.25,seed=[-.9,-.9,-1.8,.5])
    delta=delta_MeV/H
    cache=Ref(copy(seed))
    function second(T)
        c=solve_fu(T,eB,p,ctl,cache[]);cache[]=c.x
        a=solve_fu(T-delta,eB,p,ctl,c.x);b=solve_fu(T+delta,eB,p,ctl,c.x)
        (a.masses[1]+b.masses[1]-2c.masses[1])/delta^2
    end
    lo=(guess-5)/H;hi=(guess+5)/H;fa=second(lo);fb=second(hi)
    fa*fb<0 || error("Fu peak not bracketed: $fa $fb")
    for _ in 1:23
        mid=(lo+hi)/2;fm=second(mid)
        if signbit(fm)==signbit(fa)
            lo=mid;fa=fm
        else
            hi=mid;fb=fm
        end
    end
    T=(lo+hi)/2;c=solve_fu(T,eB,p,ctl,cache[])
    a=solve_fu(T-delta,eB,p,ctl,c.x);b=solve_fu(T+delta,eB,p,ctl,c.x)
    (;T_MeV=T*H,masses_MeV=c.masses.*H,Phi=c.x[4],x_state=vcat(c.x[1:3],c.x[4],c.x[4]),
        residual=c.residual,minus_dMu_dT=(a.masses[1]-b.masses[1])/(2delta),
        order=ctl.order,n_max=ctl.n_max,delta_MeV)
end

function main_fu(output)
    p=literature_params()
    out=Dict{String,Any}("generated_at"=>string(now()),"julia_version"=>string(VERSION),
        "source"=>"arXiv:1306.5804v1 Eqs.5-11 and Sec.III",
        "m_pi_MeV"=>135.,"T0_MeV"=>270.,"route"=>"landau_legacy N=10 polynomial Polyakov",
        "cases"=>Any[])
    for (field,target) in ((0.,224.),(20.,240.))
        eB=field*(135/H)^2
        rows=Any[]
        for (order,nmax,d) in ((32,96,.5),(48,160,.25),(64,240,.125))
            ctl=MR.Controls(;order,n_max=nmax,route=:landau_legacy)
            r=fu_peak(eB,p,ctl,target;delta_MeV=d)
            push!(rows,r)
            println("Fu field=",field," peak=",r.T_MeV," slope=",r.minus_dMu_dT," order=",order);flush(stdout)
        end
        push!(out["cases"],(;field_mpi2=field,eB_GeV2=field*.135^2,literature_T_MeV=target,checks=rows))
        open(io->JSON.print(io,out,2),output,"w")
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    main_fu(only(ARGS))
end
