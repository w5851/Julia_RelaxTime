module PrecompileRegistry

using StaticArrays
using ForwardDiff
using NLsolve
using ..Models: create_model, solve_constraint, FixedMu, HADRON_SEED_5
using ..Models: thermo_derivatives, mass_derivatives, bulk_viscosity_coefficients
using ..Models: run_scan_pipeline
using ..Models: chi1_B, chi2_B, chi3_B, chi4_B, chi2_Q, chi2_S, chi11_BQ, chi11_BS, chi11_QS
using ..Models: chi_BQS, conserved_charge_susceptibility, cumulant_B, cumulant_BQS
using ..Models: baryon_Ssigma, baryon_kappa_sigma2, flavor_pressure_derivatives
using Main.Constants_PNJL: ħc_MeV_fm

const _VALID_CAPABILITIES = Set([
    :gap_solver_ad,
    :thermo_derivatives_ad,
    :conserved_charge_highorder,
    :ad_shape_stabilization,
    :scan_pipeline_cli,
])

const _PROFILE_CAPABILITIES = Dict(
    :smoke => [:gap_solver_ad],
    :test => [:gap_solver_ad, :thermo_derivatives_ad],
    :scan => [:gap_solver_ad, :thermo_derivatives_ad, :scan_pipeline_cli],
    :core => [:gap_solver_ad, :thermo_derivatives_ad, :conserved_charge_highorder, :ad_shape_stabilization, :scan_pipeline_cli],
    :all => [:gap_solver_ad, :thermo_derivatives_ad, :conserved_charge_highorder, :ad_shape_stabilization, :scan_pipeline_cli],
    :full => [:gap_solver_ad, :thermo_derivatives_ad, :conserved_charge_highorder, :ad_shape_stabilization, :scan_pipeline_cli],
)

@inline function list_precompile_capabilities()
    return sort!(collect(_VALID_CAPABILITIES))
end

@inline function list_precompile_profile(profile::Symbol)
    haskey(_PROFILE_CAPABILITIES, profile) || throw(ArgumentError("unknown precompile profile: $(profile)"))
    return copy(_PROFILE_CAPABILITIES[profile])
end

@inline function _safe_run!(name::Symbol, f::Function; strict::Bool)
    try
        f()
    catch err
        if strict
            rethrow(err)
        end
        @warn "precompile capability failed" capability=name exception=(err, catch_backtrace())
    end
    return nothing
end

function _cap_gap_solver_ad()
    model = create_model(:PNJL)
    T_fm = 150.0 / ħc_MeV_fm
    mu_fm = 250.0 / ħc_MeV_fm

    solve_constraint(
        model,
        FixedMu(),
        T_fm;
        μ_fm=mu_fm,
        seed_guess=HADRON_SEED_5,
        p_num=6,
        t_num=3,
        residual_norm_max=1e-4,
        iterations=24,
    )

    return nothing
end

function _cap_thermo_derivatives_ad()
    T_fm = 150.0 / ħc_MeV_fm
    mu_fm = 250.0 / ħc_MeV_fm

    thermo_derivatives(T_fm, mu_fm; xi=0.0, p_num=6, t_num=3)
    mass_derivatives(T_fm, mu_fm; xi=0.0, p_num=6, t_num=3)
    bulk_viscosity_coefficients(T_fm, mu_fm; xi=0.0, p_num=6, t_num=3)

    return nothing
end

function _cap_conserved_charge_highorder()
    T_fm = 0.58
    muB_fm = 0.21
    muQ_fm = 0.04
    muS_fm = 0.01
    V = 64.0
    kwargs = (; xi=0.0, p_num=8, t_num=4)

    # Match high-order unit test hot paths (large-grid signatures).
    T_hi = 140.0 / ħc_MeV_fm
    muB_hi = 360.0 / ħc_MeV_fm
    chi1_B(T_hi, muB_hi; xi=0.0, p_num=48, t_num=12)
    chi2_B(T_hi, muB_hi; xi=0.0, p_num=48, t_num=12)
    chi3_B(T_hi, muB_hi; xi=0.0, p_num=48, t_num=12)
    chi4_B(T_hi, muB_hi; xi=0.0, p_num=48, t_num=12)

    T_map = 130.0 / ħc_MeV_fm
    muB_map = 300.0 / ħc_MeV_fm
    muq_map = muB_map / 3
    thermo_derivatives(T_map, muq_map; xi=0.0, p_num=48, t_num=12)
    chi1_B(T_map, muB_map; xi=0.0, p_num=48, t_num=12)

    chi2_B(T_fm, muB_fm; kwargs...)
    chi3_B(T_fm, muB_fm; kwargs...)
    chi4_B(T_fm, muB_fm; kwargs...)
    chi2_Q(T_fm, muB_fm, muQ_fm, muS_fm; kwargs...)
    chi2_S(T_fm, muB_fm, muQ_fm, muS_fm; kwargs...)
    chi11_BQ(T_fm, muB_fm, muQ_fm, muS_fm; kwargs...)
    chi11_BS(T_fm, muB_fm, muQ_fm, muS_fm; kwargs...)
    chi11_QS(T_fm, muB_fm, muQ_fm, muS_fm; kwargs...)

    chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(0, 3, 0), kwargs...)
    chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(0, 0, 4), kwargs...)
    conserved_charge_susceptibility(T_fm, muB_fm, 0.0, 0.0; orders=(2, 0, 0), kwargs...)
    cumulant_B(T_fm, muB_fm, V; order=4, kwargs...)
    cumulant_BQS(T_fm, muB_fm, muQ_fm, muS_fm, V; orders=(1, 1, 0), kwargs...)
    baryon_Ssigma(T_fm, muB_fm; kwargs...)
    baryon_kappa_sigma2(T_fm, muB_fm; kwargs...)

    # Shape stabilization for StaticArrays + AD wrappers
    mu_vec = SVector{3, Float64}(muB_fm, muQ_fm, muS_fm)
    flavor_pressure_derivatives(T_fm, mu_vec; order=2, kwargs...)

    return nothing
end

function _cap_ad_shape_stabilization()
    # Stabilize ForwardDiff + StaticArrays specializations for 3/5-dimensional
    # AD paths used by high-order susceptibilities and solver wrappers.
    f5 = x -> SVector(
        x[1]^2 + x[2],
        x[2]^2 + x[3],
        x[3]^2 + x[4],
        x[4]^2 + x[5],
        x[5]^2 + x[1],
    )
    ForwardDiff.jacobian(f5, rand(5))

    f3 = x -> SVector(
        x[1]^2 + x[2],
        x[2]^2 + x[3],
        x[3]^2 + x[1],
    )
    ForwardDiff.jacobian(f3, rand(3))

    # Trigger NLSolversBase.OnceDifferentiable(..., ForwardDiff.Chunk{5})
    function residual5!(F, x)
        @inbounds for i in 1:5
            F[i] = x[i]^2 - 1.0
        end
        return nothing
    end
    nlsolve(residual5!, ones(5); autodiff=:forward, method=:newton, iterations=1, ftol=1e-6, xtol=1e-6)

    return nothing
end

function _cap_scan_pipeline_cli()
    out_root = mktempdir()

    run_scan_pipeline(
        :tmu;
        model_kind=:PNJL,
        T_values=[150.0],
        mu_values=[0.0, 100.0],
        xi_values=[0.0],
        output_path=joinpath(out_root, "tmu.csv"),
        overwrite=true,
    )

    run_scan_pipeline(
        :trho;
        model_kind=:PNJL,
        T_values=[150.0],
        rho_values=[0.1, 0.2],
        xi_values=[0.0],
        output_path=joinpath(out_root, "trho.csv"),
        overwrite=true,
    )

    return nothing
end

function run_precompile_capability(capability::Symbol; strict::Bool=false)
    capability in _VALID_CAPABILITIES || throw(ArgumentError("unknown precompile capability: $(capability)"))

    if capability === :gap_solver_ad
        _safe_run!(capability, _cap_gap_solver_ad; strict=strict)
    elseif capability === :thermo_derivatives_ad
        _safe_run!(capability, _cap_thermo_derivatives_ad; strict=strict)
    elseif capability === :conserved_charge_highorder
        _safe_run!(capability, _cap_conserved_charge_highorder; strict=strict)
    elseif capability === :ad_shape_stabilization
        _safe_run!(capability, _cap_ad_shape_stabilization; strict=strict)
    elseif capability === :scan_pipeline_cli
        _safe_run!(capability, _cap_scan_pipeline_cli; strict=strict)
    else
        throw(ArgumentError("unhandled precompile capability: $(capability)"))
    end

    return nothing
end

function run_precompile_profile(profile::Symbol=:smoke; strict::Bool=false)
    capabilities = get(_PROFILE_CAPABILITIES, profile, nothing)
    capabilities === nothing && throw(ArgumentError("unknown precompile profile: $(profile)"))

    for capability in capabilities
        run_precompile_capability(capability; strict=strict)
    end

    return nothing
end

end # module PrecompileRegistry
