module ScanQuality

export assess_point_quality

const _TAU_FIELDS = (:u, :d, :s, :ubar, :dbar, :sbar)

function assess_point_quality(tau::NamedTuple, eta_over_s::Float64, sigma_over_T::Float64; ratio_threshold::Float64=6.0)
    for field in _TAU_FIELDS
        val = getfield(tau, field)
        if !isfinite(val) || val <= 0.0
            return true, string("tau_invalid_", field), val
        end
    end

    ratio = max(tau.u, tau.ubar) / min(tau.u, tau.ubar)
    if !isfinite(ratio) || ratio > ratio_threshold
        return true, "tau_u_ubar_ratio_high", ratio
    end

    if !isfinite(eta_over_s)
        return true, "eta_over_s_nonfinite", eta_over_s
    end

    if !isfinite(sigma_over_T)
        return true, "sigma_over_T_nonfinite", sigma_over_T
    end

    return false, "ok", ratio
end

end
