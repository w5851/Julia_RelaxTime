# Ward-matched pi/K poles on frozen zero-density backgrounds.
using CSV, JSON3, SHA, LinearAlgebra, NLsolve
const SPECTRA_ROOT = normpath(joinpath(@__DIR__,"..",".."))
if !isdefined(Main,:Models)
    include(joinpath(SPECTRA_ROOT,"src","models","Models.jl"))
end

function spectral_bisection(f,a,b;tol=1e-10)
    fa,fb=f(a),f(b)
    isfinite(fa) && isfinite(fb) && signbit(fa)!=signbit(fb) || error("invalid spectral bracket")
    while b-a>tol
        x=(a+b)/2;fx=f(x)
        isfinite(fx) || error("nonfinite spectral residual")
        if signbit(fx)==signbit(fa);a=x;fa=fx;else;b=x;end
    end
    return (a+b)/2
end

function spectral_pole(meson,st,previous)
    rt=Main.RelaxTime
    qp,tp,K=st.quark_params,st.thermo_params,st.K_coeffs
    threshold=rt.MottTransition.mott_threshold_mass(meson,qp)
    lower=meson===:pi ? 0.003 : abs(qp.m.s-qp.m.u)+0.002
    f(m,w)=rt.MesonMass.meson_mass_equation(meson,m,w,0.0,qp,tp,K)
    a,b=real(f(lower,0.0)),real(f(threshold,0.0))
    if isfinite(a) && isfinite(b) && signbit(a)!=signbit(b)
        mass=spectral_bisection(m->real(f(m,0.0)),lower,threshold)
        residual=abs(f(mass,0.0))
        residual<2e-7 || error("bound pole residual: $meson $residual")
        return (mass=mass,width=0.0,residual=residual,status="bound")
    end
    default=rt.MesonMass.default_meson_mass_guess(meson,qp)
    guesses=previous===nothing ? Tuple{Float64,Float64}[] : [(previous.mass,previous.width)]
    append!(guesses,[(default,0.0),(max(threshold*1.03,default),0.05),
        (max(threshold*1.1,default),0.2),(max(threshold*1.5,default),0.6)])
    for seed in guesses
        try
            sol=nlsolve((out,x)->begin
                out[1],out[2]=reim(f(x[1],x[2]));nothing
            end,collect(seed);iterations=80,ftol=1e-10,xtol=1e-11)
            mass,width=sol.zero
            residual=abs(f(mass,width))
            if mass>lower && width>=-1e-8 && residual<2e-7
                return (mass=mass,width=max(0.0,width),residual=residual,
                    status=mass<threshold ? "bound" : "resonance")
            end
        catch e
            e isa InterruptException && rethrow()
        end
    end
    error("no valid real-axis pole: $meson xi=$(tp.ξ) T=$(tp.T)")
end

function spectral_main(args)
    length(args)==4 || error("INPUT OUTPUT PROFILE XI_LIST")
    input,outdir,profile_name,xi_arg=args
    profile=Models.TransportSamplingProfiles.load_profile(profile_name)
    xis=parse.(Float64,split(xi_arg,','))
    rows=collect(CSV.File(input))
    mkpath(outdir)
    output=NamedTuple[]
    hc=Main.Constants_PNJL.ħc_MeV_fm
    for xi in xis
        by=sort!([r for r in rows if abs(r.xi-xi)<1e-10];by=r->r.T_MeV)
        length(by)==201 || error("incomplete frozen spectral scan")
        previous=Dict{Symbol,NamedTuple}()
        for row in by
            qp=(m=(u=row.m_u_MeV/hc,d=row.m_u_MeV/hc,s=row.m_s_MeV/hc),μ=(u=0.0,d=0.0,s=0.0))
            tp=(T=row.T_MeV/hc,Φ=Float64(row.Phi),Φbar=Float64(row.Phibar),ξ=Float64(row.xi))
            state=Models.TransportSamplingProfiles.fixed_point_state(profile,qp,tp)
            for meson in (:pi,:K)
                pole=spectral_pole(meson,state,get(previous,meson,nothing))
                previous[meson]=pole
                threshold=Main.RelaxTime.MottTransition.mott_threshold_mass(meson,qp)
                inverse=Main.RelaxTime.MesonMass.meson_mass_equation(meson,threshold,0.0,0.0,
                    state.quark_params,state.thermo_params,state.K_coeffs)
                push!(output,(T_MeV=row.T_MeV,muB_MeV=0.0,xi=xi,meson=String(meson),
                    mass_MeV=pole.mass*hc,width_MeV=pole.width*hc,pole_residual=pole.residual,
                    status=pole.status,threshold_MeV=threshold*hc,threshold_inverse=real(inverse),
                    m_u_MeV=row.m_u_MeV,m_s_MeV=row.m_s_MeV,Phi=row.Phi,Phibar=row.Phibar,
                    sampling_profile=profile.id))
            end
            if row.T_MeV%25==0
                println("SPECTRUM xi=",xi," T=",row.T_MeV);flush(stdout)
            end
        end
        CSV.write(joinpath(outdir,"spectra.csv.writing"),output)
        mv(joinpath(outdir,"spectra.csv.writing"),joinpath(outdir,"spectra.csv");force=true)
    end
    write(joinpath(outdir,"manifest.json"),JSON3.write((input_sha256=bytes2hex(sha256(read(input))),
        profile=Models.TransportSamplingProfiles.effective_config(profile),xi_values=xis,rows=length(output),
        background="frozen p96/cos24, thermal upper limit25 inverse fm",source_commit=get(ENV,"PNJL_SOURCE_TREE_SHA256","unrecorded"))))
end

if abspath(PROGRAM_FILE)==@__FILE__
    spectral_main(ARGS)
end
