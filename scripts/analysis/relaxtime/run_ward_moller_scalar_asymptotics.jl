#!/usr/bin/env julia
# Continuous near-threshold positive-square integral; no transport cache is read.
# julia --project=SOURCE run_ward_moller_scalar_asymptotics.jl SOURCE OUT [fine|reference]
include(joinpath(@__DIR__, "run_ward_moller_mott_response.jl"))

const SC_CS = MM_RT.TotalCrossSection
const SC_PROP = MM_RT.MesonPropagator
const SC_XMAX = 0.01 # invariant excess s-4M_u^2, in fm^-2

# At fixed s, the two-particle measure in total momentum Q and CM angles
# contains E1 E2 beta/(4E). Multiplying by the medium Møller flux gives
# x/(8E), with x=s-4m^2 and E=sqrt(s+Q^2). The common azimuth and
# normalization are already included in the production cache prefactor.
function sc_moller_weight(state, x, prefactor; p=64, c=16, angles=12, boltzmann=false)
    qp, tp = state.quark_params, state.thermo_params
    m = qp.m.u
    root_s = sqrt(4m^2+x)
    k = sqrt(x)/2
    Qs, ws = MM_AS._cached_semi_infinite_momentum_grid(p, 20.0)
    cs, wc = MM_AS._cached_interval_grid(-1.0, 1.0, c)
    bs, wb = MM_AS._cached_interval_grid(-1.0, 1.0, angles)
    phis, wp = MM_AS._cached_interval_grid(0.0, 2pi, angles)
    μ1, μ2 = MM_AS.get_mu(:u, qp), MM_AS.get_mu(:ubar, qp)
    total = 0.0
    for (Q, wQ) in zip(Qs, ws)
        E = hypot(root_s, Q)
        γ = E/root_s
        for (ct, wt) in zip(cs, wc), (cb, wcb) in zip(bs, wb), (ϕ, wϕ) in zip(phis, wp)
            product = if boltzmann
                exp(-E/tp.T)
            else
                boost = Q*k/root_s*cb
                E1, E2 = E/2+boost, E/2-boost
                transverse = k*sqrt((1-cb^2)*(1-ct^2))*cos(ϕ)
                z1, z2 = (Q/2+γ*k*cb)*ct-transverse, (Q/2-γ*k*cb)*ct+transverse
                p1, p2 = sqrt(max(0.0, E1^2-m^2)), sqrt(max(0.0, E2^2-m^2))
                c1, c2 = p1>0 ? clamp(z1/p1,-1,1) : 0.0, p2>0 ? clamp(z2/p2,-1,1) : 0.0
                MM_AS.distribution_with_anisotropy(:u,p1,m,μ1,tp.T,tp.Φ,tp.Φbar,tp.ξ,c1) *
                    MM_AS.distribution_with_anisotropy(:ubar,p2,m,μ2,tp.T,tp.Φ,tp.Φbar,tp.ξ,c2)
            end
            total += wQ*Q^2/E*wt*wcb*wϕ*product
        end
    end
    prefactor*total/8
end

function sc_boltzmann_oracle(state, x, prefactor; p=64)
    s = 4state.quark_params.m.u^2+x
    q, w = MM_AS._cached_semi_infinite_momentum_grid(p, 20.0)
    pi*prefactor*sum(wi*Qi^2/hypot(sqrt(s),Qi)*exp(-hypot(sqrt(s),Qi)/state.thermo_params.T)
        for (Qi,wi) in zip(q,w))
end

function sc_weight_table(state, prefactor; fine=false, degree=8)
    opts = fine ? (p=96,c=24,angles=20) : (p=64,c=16,angles=12)
    xs = [SC_XMAX/2*(1-cospi(j/degree)) for j in 0:degree]
    ys = [sc_moller_weight(state,x,prefactor; opts...) for x in xs]
    weights = [(-1.0)^j*(j in (0,degree) ? 0.5 : 1.0) for j in 0:degree]
    function weight(x)
        i = findfirst(==(x), xs)
        i !== nothing && return ys[i]
        a = weights ./ (x .- xs)
        sum(a .* ys)/sum(a)
    end
    xcheck = 0.003731
    defect = abs(weight(xcheck)/sc_moller_weight(state,xcheck,prefactor; opts...)-1)
    defect < 1e-8 || error("continuous Møller weight interpolation failed: $defect")
    (; weight, defect)
end

function sc_occupation(y, m, μ, tp)
    E = sqrt(m^2+y)
    T,P,Pb,xi = BigFloat.((tp.T,tp.Φ,tp.Φbar,tp.ξ))
    function occ(chemical,p,pb)
        z = exp(-(E-chemical)/T)
        num = p*z+2pb*z^2+z^3
        den = 1+3p*z+3pb*z^2+z^3
        d = (-(p*z+4pb*z^2+3z^3)*den+num*(3p*z+6pb*z^2+3z^3))/(T*den^2)
        num/den + xi*y/(6E)*d
    end
    (1-occ(μ,P,Pb)-occ(-μ,Pb,P))/E
end

function sc_pair_plan(st, flavor, n)
    qp,tp = st.quark_params, st.thermo_params
    m, μ = BigFloat.((getproperty(qp.m,flavor),getproperty(qp.μ,flavor)))
    L = BigFloat(Main.Constants_PNJL.Λ_inv_fm)
    z,w = SC_CS.standard_nodes_weights(n)
    ys = [(L*(BigFloat(zi)+1)/2)^2 for zi in z]
    ws = L/2 .* BigFloat.(w)
    Fs = [sc_occupation(y,m,μ,tp) for y in ys]
    (; m, μ, tp, L, ys, ws, Fs, F0=sc_occupation(big"0",m,μ,tp), I0=sum(ws.*Fs))
end

# Same subtraction and on-shell discontinuity as PairLoopQuadrature,
# represented at high precision near a vanishing mixed determinant.
function sc_pair_loop(plan, q2)
    if q2 == 0
        return complex(2plan.I0, big"0")
    elseif q2 < 0
        κ = sqrt(-q2)
        reg = sum(w*(F-plan.F0)/(y-q2) for (y,w,F) in zip(plan.ys,plan.ws,plan.Fs))
        return complex(2plan.I0+2q2*reg-2κ*plan.F0*atan(plan.L/κ), big"0")
    end
    q = sqrt(q2)
    F = sc_occupation(q2,plan.m,plan.μ,plan.tp)
    reg = sum(w*(v-F)/(y-q2) for (y,w,v) in zip(plan.ys,plan.ws,plan.Fs))
    re = 2plan.I0+2q2*reg+q*F*log1p(-2q/(plan.L+q))
    complex(re,big(pi)*q*F)
end

function sc_components(st, x, plans)
    u,s = plans
    s_hat = 4u.m^2+x
    qp,kc = st.quark_params,st.K_coeffs
    pref = -BigFloat(Main.Constants_PNJL.N_color)/(8big(pi)^2)
    pu = pref*(2BigFloat(qp.A.u)-x*sc_pair_loop(u,x/4))
    ps = pref*(2BigFloat(qp.A.s)+(4s.m^2-s_hat)*sc_pair_loop(s,(s_hat-4s.m^2)/4))
    dK = BigFloat(kc.det_K_minus)
    M00 = BigFloat(kc.K0_minus)-4dK/3*(pu+2ps)
    M08 = BigFloat(kc.K08_minus)+4dK/3*sqrt(BigFloat(2))* (pu-ps)
    M88 = BigFloat(kc.K8_minus)-4dK/3*(2pu+ps)
    det = M00*M88-M08^2
    j0,j8 = sqrt(BigFloat(2)/3),inv(sqrt(BigFloat(3)))
    num = 2dK*(j0*(M00*j0+M08*j8)+j8*(M08*j0+M88*j8))
    det,num
end

function sc_kernel(st; n=128, degree=14)
    setprecision(BigFloat,192) do
        plans = (sc_pair_plan(st,:u,n),sc_pair_plan(st,:s,n))
        function components(y)
            x = BigFloat(SC_XMAX)*y
            det,num = sc_components(st,x,plans)
            # Im(det) and Im(numerator) share the light scalar cut factor x^(3/2).
            z = x == 0 ? big"1e-36" : x
            imagdet,imagnum = x == 0 ? map(imag,sc_components(st,z,plans)) : (imag(det),imag(num))
            [real(det), imagdet/(z*sqrt(z)), real(num), imagnum/(z*sqrt(z))]
        end
        ys = [(1-cos(big(pi)*j/degree))/2 for j in 0:degree]
        v = [y^j for y in ys, j in 0:degree]
        coeffs = v \ reduce(vcat,[transpose(components(y)) for y in ys])
        err = maximum(maximum(abs.(components(y).- [evalpoly(y,col) for col in eachcol(coeffs)]))
            for y in BigFloat.((0.0137,0.0731,0.2193,0.4387,0.6713,0.8921,0.9783)))
        err < big"1e-12" || error("pair kernel representation unresolved: $err")
        c = coeffs[:,1]
        dc = [j*c[j+1] for j in 1:degree]
        yroot = -c[1]/c[2]
        for _ in 1:15
            dy = evalpoly(yroot,c)/evalpoly(yroot,dc)
            yroot -= dy
            abs(dy)<big"1e-43" && break
        end
        abs(evalpoly(yroot,c)) < big"1e-40" || error("scalar resonance zero not resolved")
        # Translate polynomial to the root; suppress only its solved zero.
        real_coeffs = Float64.([j==0 ? big"0" :
            sum(c[k+1]*binomial(k,j)*yroot^(k-j) for k in j:degree) for j in 0:degree])
        xr = Float64(yroot*BigFloat(SC_XMAX))
        slope = real_coeffs[2]/SC_XMAX
        width = xr>0 ? abs(xr^1.5*Float64(evalpoly(yroot,coeffs[:,2]))/slope) : 0.0
        finite_coeffs = [Float64.(col) for col in eachcol(coeffs)]
        (; xr,width,real_coeffs, coeffs=finite_coeffs, interpolation_defect=Float64(err))
    end
end

function sc_block(state, s)
    qp,tp = state.quark_params,state.thermo_params
    z,w = SC_CS.standard_nodes_weights(24)
    m,μ = qp.m.u,qp.μ.u
    E = sqrt(s)/2
    sum(wi*SC_CS.combined_final_state_factor(:u,:ubar,E,E,m,m,μ,μ,
        tp.T,tp.Φ,tp.Φbar,tp.ξ,zi) for (zi,wi) in zip(z,w))/2
end

function sc_cross_square(st, x, dx, kernel)
    y, x32 = x/SC_XMAX, x*sqrt(x)
    c = kernel.coeffs
    D = complex(evalpoly(dx/SC_XMAX,kernel.real_coeffs),x32*evalpoly(y,c[2]))
    N = complex(evalpoly(y,c[3]),x32*evalpoly(y,c[4]))
    s = 4st.quark_params.m.u^2+x
    abs2(N/D)*x^2/(16pi*s)*sc_block(st,s)
end

function sc_gl(f,a,b,n)
    z,w = SC_CS.standard_nodes_weights(n)
    (b-a)/2*sum(wi*f((a+b)/2+(b-a)/2*zi) for (zi,wi) in zip(z,w))
end

function sc_quad(f,a,b; tol=1e-21, rtol=2e-11, depth=0)
    coarse,fine = sc_gl(f,a,b,16),sc_gl(f,a,b,32)
    defect = abs(fine-coarse)
    defect <= tol+rtol*abs(fine) && return (fine,defect)
    depth<18 || error("continuous integral not resolved at [$a,$b]")
    mid = (a+b)/2
    left = sc_quad(f,a,mid;tol=tol/2,rtol,depth=depth+1)
    right = sc_quad(f,mid,b;tol=tol/2,rtol,depth=depth+1)
    left .+ right
end

function sc_integrate(f, r, w; rtol=2e-11)
    edges = [0.0,SC_XMAX]
    peak,error = 0.0,0.0
    if 0<r<SC_XMAX
        w>0 || error("nonpositive scalar width")
        a,b = max(0.0,r-64w),min(SC_XMAX,r+64w)
        tlo,thi = atan((a-r)/w),atan((b-r)/w)
        peak,error = sc_quad(t->begin
            dx = w*tan(t)
            f(r+dx,dx)*w/cos(t)^2
        end,tlo,thi;rtol)
        # Geometric tails in the shifted variable retain the full integral.
        tails = vcat([r-w*2.0^j for j in 6:55 if 0<r-w*2.0^j<a],
                     [r+w*2.0^j for j in 6:55 if b<r+w*2.0^j<SC_XMAX])
        append!(edges, [a,b,tails...])
        sort!(unique!(edges))
        intervals = [(l,u) for (l,u) in zip(edges[1:end-1],edges[2:end])
            if u<=a || l>=b]
    else
        s = max(abs(r),1e-12)
        append!(edges,[s*2.0^j for j in -3:55 if s*2.0^j<SC_XMAX])
        sort!(unique!(edges))
        intervals = collect(zip(edges[1:end-1],edges[2:end]))
    end
    value = peak
    for (a,b) in intervals
        v,e = sc_quad(x->f(x,x-r),a,b;rtol)
        value += v
        error += e
    end
    (; I=value, peak, quadrature_defect=error)
end

function sc_prefactor(runtime, st)
    problem = MM_AS._prepare_thermal_cache_problem(:uubar_to_uubar,st.quark_params,
        st.thermo_params,st.K_coeffs;runtime.tau_kwargs...)
    problem.sampling.relative_velocity === :medium_moller || error("wrong scalar flux")
    problem.plan.prefactor
end

# Independent local coefficient from the cutoff-subtracted pair kernel at
# the threshold, without fitting the integrated I_S or the peak areas.
function sc_limit(bg, center, runtime)
    st = bg(center).state
    qp,tp,kc = st.quark_params,st.thermo_params,st.K_coeffs
    function components(state,x)
        q,t,k = state.quark_params,state.thermo_params,state.K_coeffs
        s = 4q.m.u^2+x
        p = map((:u,:s)) do flavor
            m,μ,A = getproperty(q.m,flavor),getproperty(q.μ,flavor),getproperty(q.A,flavor)
            plan = MM_PAIR.PairLoopPlan(m,μ,t;n_nodes=256)
            -Main.Constants_PNJL.N_color/(8pi^2) * (2A+(flavor===:u ? -x : 4m^2-s)*MM_PAIR.pair_b0_pv(plan,s))
        end
        M = MM_EC.mixing_matrix_elements(p...,k,:S)
        D = M.M00*M.M88-M.M08^2
        j0,j8 = SC_PROP._current_components(:u,:ubar)
        N = 2M.det_K*(j0*(M.M00*j0+M.M08*j8)+j8*(M.M08*j0+M.M88*j8))
        D,N
    end
    d0,num = components(st,0.0)
    function limit(du, dx)
        a = real(components(bg(center+du).state,0.0)[1]-components(bg(center-du).state,0.0)[1])/(2du)
        d = components(st,dx)[1]
        b,c = real(d-d0)/dx, imag(d)/dx^1.5
        a,b,c
    end
    a1,b1,c1 = limit(2e-5,2e-5)
    a2,b2,c2 = limit(1e-5,1e-5)
    a,b,c = (4a2-a1)/3,2b2-b1,2c2-c1
    ratio = a/b
    ratio != 0 && c != 0 || error("degenerate scalar threshold")
    pref = sc_prefactor(runtime,st)
    g0 = sc_moller_weight(st,0.0,pref;p=96,c=24,angles=20)
    s0 = 4qp.m.u^2
    F0 = abs2(num)*g0*sc_block(st,s0)/(16pi*s0)
    expected = 3pi/4*F0*abs(ratio)^1.5/abs(b*c)
    checks = [let q = sc_moller_weight(st,x,pref;boltzmann=true)
       q/sc_boltzmann_oracle(st,x,pref)-1 end for x in (0.0,1e-4,0.01)]
    maximum(abs,checks) < 3e-11 || error("medium flux Jacobian failed the angular oracle")
    (; a,b,c,F0,expected,resonance_sign=-sign(ratio),
        sampling_check=checks, coefficient_defects=(a1/a-1,b1/b-1,c1/c-1))
end

function sc_curvature(getrow, u, step)
    nodes = [getrow(u+j*step) for j in -2:2]
    terms = (-1.0,16.0,-30.0,16.0,-1.0)
    # Subtract the center before accumulating the cancelling stencil.
    estimate = sum(w*(r.I-nodes[3].I) for (w,r) in zip(terms,nodes))/(12step^2)
    defect = sum(abs(w)*r.quadrature_defect for (w,r) in zip(terms,nodes))/(12step^2)
    estimate,defect
end

function sc_main(args)
    length(args) in (2,3) || error("SOURCE OUT [fine|reference]")
    outdir,tier = abspath(args[2]),get(args,3,"fine")
    tier in ("fine","reference") || error("invalid scalar precision tier")
    folder = joinpath(outdir,"scalar_"*tier)
    ispath(folder) && error("refusing to overwrite $folder")
    mkpath(folder)
    LinearAlgebra.BLAS.set_num_threads(1)
    runtime = mm_runtime(mm_profile(),"production")
    target = mm_target(:sigma)
    center_path = joinpath(outdir,"center_sigma.json")
    location = if isfile(center_path)
        record = JSON3.read(read(center_path,String))
        record.metadata.core_commit == MM_SHA || error("unexpected center provenance")
        record.location
    else
        mm_location(runtime,target,folder)
    end
    center = Float64(location.xi)
    cache = Dict{Float64,NamedTuple}()
    bg(x) = get!(cache,x) do
        mm_equilibrium(runtime,target,x;seed=location.center_state)
    end
    lim = sc_limit(bg,center,runtime)
    nodes = Dict{Float64,NamedTuple}()
    reftier = tier=="reference"
    println("SCALAR xi_M=$center sign=$(lim.resonance_sign) target=$(lim.expected)")
    flush(stdout)
    function getrow(u)
        get!(nodes,u) do
            st = bg(center+u).state
            kernel = sc_kernel(st;n=reftier ? 256 : 128,degree=reftier ? 18 : 14)
            kernel.xr > 0 || error("failed to retain the resonance branch at $u")
            w = sc_weight_table(st,sc_prefactor(runtime,st);fine=reftier,degree=reftier ? 10 : 8)
            integ = sc_integrate((x,dx)->x*w.weight(x)*sc_cross_square(st,x,dx,kernel),kernel.xr,kernel.width;
                rtol=reftier ? 2e-12 : 2e-11)
            (; u, xr=kernel.xr, width=kernel.width, integ..., kernel.interpolation_defect, weight_defect=w.defect)
        end
    end
    rows = NamedTuple[]
    steps = 2e-4 ./ (2.0 .^ (reftier ? (1:2:9) : (0:10)))
    for h in steps
        u = lim.resonance_sign*h
        coarse,_ = sc_curvature(getrow,u,h/8)
        fine,quadrature = sc_curvature(getrow,u,h/16)
        r = getrow(u)
        row = (; h, u, r.xr,r.width,r.I,r.peak, curvature=fine, scaled_curvature=sqrt(h)*fine,
            step_defect=abs(fine-coarse),quadrature_defect=quadrature,
            relative_to_limit=sqrt(h)*fine/lim.expected-1)
        push!(rows,row)
        CSV.write(joinpath(folder,"curvatures.csv"),rows)
        CSV.write(joinpath(folder,"integrals.csv"),sort!(collect(values(nodes));by=x->x.u))
        println("SCALAR $tier h=$h scaled=$(row.scaled_curvature) step=$(row.step_defect) peak=$(r.peak)")
        flush(stdout)
    end
    mm_json(joinpath(folder,"summary.json"),(; status="complete",metadata=mm_metadata(runtime),
        scalar_driver_sha256=mm_digest(@__FILE__),center,lim,tier,xmax_inv_fm2=SC_XMAX,
        loops=reftier ? 256 : 128,scope="continuous uubar elastic s-channel mixed scalar positive square",
        determinant_floor=false,files=Dict(name=>mm_digest(joinpath(folder,name))
            for name in ("curvatures.csv","integrals.csv"))))
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    sc_main(ARGS)
end
