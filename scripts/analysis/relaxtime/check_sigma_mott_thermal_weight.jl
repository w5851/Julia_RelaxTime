#!/usr/bin/env julia
# Independent original-particle-variable check of the continuum Jacobian.
include(joinpath(@__DIR__, "run_sigma_mott_thermal_audit.jl"))

function boltzmann_original_variables(state,x,prefactor;n=128)
    m,T=state.qp.m.u,state.tp.T
    s=4m^2+x
    p,w=AS._cached_semi_infinite_momentum_grid(n,10.0)
    A=(s-2m^2)/2
    B=sqrt(s*x)/2
    total=0.0
    for (momentum,weight) in zip(p,w)
        E1=sqrt(m^2+momentum^2)
        E2lo=(A*E1-momentum*B)/m^2
        E2hi=(A*E1+momentum*B)/m^2
        delta=E2hi-E2lo
        # Analytic integral of E2 exp(-E2/T); expm1 protects small x.
        inner=T*exp(-E2lo/T)*((E2lo+T)*(-expm1(-delta/T))-delta*exp(-delta/T))
        total+=weight*momentum*exp(-E1/T)*inner
    end
    prefactor*4pi*sqrt(x/s)/x*total
end

function check_sigma_weight()
    profile=SP_OTHER.load_profile(joinpath(INPUT_DIR,"s6.toml"))
    rows=collect(CSV.File(joinpath(dirname(OUTPUT_DIR),"spectrum_v2","background.csv")))
    row=only(r for r in rows if String(r.point_id)=="sig_c000")
    state=frozen_state(row,profile,32,16)
    # Density normalization is common and arbitrary for the Jacobian comparison.
    prefactor=1.0
    checks=NamedTuple[]
    for x in (1e-8,1e-6,1e-4,0.001,0.01,0.1,1.0)
        direct=boltzmann_original_variables(state,x,prefactor;n=128)
        pair=sigma_boltzmann_g_oracle(state,x,prefactor;p_nodes=128)
        relative=direct/pair-1
        abs(relative)<2e-8 || error("Original p1,p2 angular-delta Jacobian test failed at x=$x: $relative")
        push!(checks,(x_inv_fm2=x,original_particle_variables_g=direct,pair_variables_g=pair,
            relative_error=relative))
    end
    CSV.write(joinpath(OUTPUT_DIR,"jacobian_comparison.csv"),checks)
    write_json(joinpath(OUTPUT_DIR,"summary.json"),(status="completed_diagnostic",checks=length(checks),
        maximum_relative_error=maximum(abs(r.relative_error) for r in checks),
        formula="Original isotropic Boltzmann angular delta integral: W/x=prefactor*4pi*sqrt(x/s)/x*integral p1 dp1 exp(-E1/T) integral[E2lo,E2hi] E2 dE2 exp(-E2/T).",
        calculation_sha="2e5c411c716d5cbc524c3fc6cfb335fe05b39444",script_sha256=digest(@__FILE__),
        production_kernel_changed=false,
        files=Dict(f=>digest(joinpath(OUTPUT_DIR,f)) for f in readdir(OUTPUT_DIR))))
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    check_sigma_weight()
end
