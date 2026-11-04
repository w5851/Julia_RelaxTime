#!/usr/bin/env julia
# Diagnostic numerical identity checks, executed only on quote-remote.
using Test
include(joinpath(@__DIR__,"run_transport_mott_probe.jl"))
include(joinpath(@__DIR__,"transport_mott_thresholds.jl"))
include(joinpath(@__DIR__,"transport_mott_kaon_endpoint.jl"))

function direct_symmetric_kaon_inverse(state,mass;n_nodes=512)
    qp,tp,kc=state.qp,state.tp,state.kc
    p,w=Main.RelaxTime.GaussLegendre.gauleg(0.0,Main.Constants_PNJL.Λ_inv_fm,n_nodes)
    b0=0.0
    for (f,o) in ((:u,:s),(:s,:u))
        m,mother,mu=getproperty(qp.m,f),getproperty(qp.m,o),getproperty(qp.μ,f)
        a=mass^2+m^2-mother^2
        for (pp,ww) in zip(p,w)
            E=sqrt(m^2+pp^2)
            occupation=PAIR_OTHER.pair_occupation(pp^2,m,mu,tp.T,tp.Φ,tp.Φbar,tp.ξ)
            b0+=-2ww*pp^2*occupation*(1/(a+2mass*E)+1/(a-2mass*E))
        end
    end
    pol=-Main.Constants_PNJL.N_color/(8pi^2)*(qp.A.u+qp.A.s+((qp.m.u-qp.m.s)^2-mass^2)*b0)
    1-4MP._meson_coupling(:K,kc)*pol
end

profile=Models.TransportSamplingProfiles.load_profile(joinpath(INPUT_DIR,"base.toml"))
source=collect(CSV.File(joinpath(INPUT_DIR,"source_rows.csv")))
targets=JSON3.read(read(joinpath(INPUT_DIR,"targets.json"),String)).targets
ids=Set(String(t.point_left) for t in targets)
rows=NamedTuple[]
@testset "Other-Mott endpoint identities" begin
    for row in filter(r->String(r.point_id) in ids,source)
        st=frozen_state(row,profile,32,16)
        threshold=st.qp.m.u+st.qp.m.s
        endpoint128=kaon_endpoint_inverse(st;n_nodes=128)
        endpoint256=kaon_endpoint_inverse(st;n_nodes=256)
        @test abs(endpoint128-endpoint256)<1e-10
        for offset in (-.01,-.001)
            mass=threshold+offset
            exact=kaon_pair_inverse(st,mass;n_nodes=256)
            direct=direct_symmetric_kaon_inverse(st,mass;n_nodes=512)
            @test abs(exact-direct)<1e-8
            native=MP.meson_mass_equation(:K,mass,0.0,0.0,st.qp,st.tp,st.kc)
            push!(rows,merge(point_fields(row),(offset_inv_fm=offset,endpoint128=real(endpoint128),
                endpoint256=real(endpoint256),stable_re=real(exact),stable_im=imag(exact),
                direct_re=direct,native_re=real(native),native_im=imag(native),
                stable_direct_abs_difference=abs(exact-direct),native_stable_abs_difference=abs(exact-native))))
        end
        qp=merge(st.qp,(m=(u=st.qp.m.u,d=st.qp.m.u,s=st.qp.m.u),A=(u=st.qp.A.u,d=st.qp.A.u,s=st.qp.A.u)))
        equal=merge(st,(qp=qp,))
        mass=2qp.m.u-.001
        expected=1-4MP._meson_coupling(:K,st.kc)*equal_pair_polarization(equal,:u,:P,mass;n_nodes=256)
        @test abs(kaon_pair_inverse(equal,mass;n_nodes=256)-expected)<1e-10
    end
end
CSV.write(joinpath(OUTPUT_DIR,"endpoint_controls.csv"),rows)
write_json(joinpath(OUTPUT_DIR,"summary.json"),(status="passed_diagnostic_identity_checks",points=length(ids),
    direct_quadrature_nodes=512,pair_nodes=256,source_rows_sha256=digest(joinpath(INPUT_DIR,"source_rows.csv")),
    max_stable_direct_difference=maximum(r.stable_direct_abs_difference for r in rows),
    max_native_stable_difference=maximum(r.native_stable_abs_difference for r in rows),
    production_kernel_changed=false,script_sha256=digest(@__FILE__)))
