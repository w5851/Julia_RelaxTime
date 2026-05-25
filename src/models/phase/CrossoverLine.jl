struct CrossoverResult
	found::Bool
	T_crossover::Union{Nothing, Float64}
	rho::Union{Nothing, Float64}
	method::Symbol
	derivative_value::Union{Nothing, Float64}
	iterations::Int
	details::Dict{Symbol, Any}
end

CrossoverResult(; method::Symbol=:peak) = CrossoverResult(false, nothing, nothing, method, nothing, 0, Dict{Symbol, Any}())

const HBARC_MEV_FM = 197.327

function _build_crossover_evaluator(μ_fm::Float64, xi::Real, p_num::Int, t_num::Int,
		model_kind::Symbol, solver_backend::Symbol)
	model = create_model(model_kind)
	mu_vec = normalize_mu_vec(μ_fm)

	function evaluate_point(T::Float64; order::Int=1, var_idx::Int=1, need_rho::Bool=false)
		derivative_order = order >= 2 ? 2 : 1
		result = solve_pnjl_with_derivatives(
			T,
			μ_fm;
			order=derivative_order,
			xi=xi,
			p_num=p_num,
			t_num=t_num,
			thermo_backend=:models,
			solver_backend=solver_backend,
			derivative_backend=:taylordiff,
		)
		x = result.x
		dx_dT = result.dx_dT
		d1 = Float64(dx_dT[var_idx])
		d2 = NaN
		if order >= 2
			d2x_dT2 = result.d2x_dT2
			d2 = Float64(d2x_dT2[var_idx])
		end

		rho = nothing
		if need_rho
			_, rho_norm, _, _ = model_thermo(model, x, mu_vec, T;
				xi=xi,
				p_num=p_num,
				t_num=t_num,
			)
			rho = Float64(rho_norm)
		end

		return (d1=d1, d2=d2, rho=rho)
	end

	return evaluate_point
end

function _detect_crossover_peak(T_min::Float64, T_max::Float64,
		evaluate_point::Function, n_scan::Int, tol::Real, max_iter::Int, var_idx::Int)
	T_vals = collect(range(T_max, T_min; length=n_scan))
	abs_derivs = Float64[]
	derivs = Float64[]

	for T in T_vals
		try
			probe = evaluate_point(T; order=1, var_idx=var_idx, need_rho=false)
			push!(derivs, probe.d1)
			push!(abs_derivs, abs(probe.d1))
		catch
			push!(derivs, NaN)
			push!(abs_derivs, NaN)
		end
	end

	valid_mask = .!isnan.(abs_derivs)
	!any(valid_mask) && return CrossoverResult(method=:peak)

	local_maxima = Int[]
	for i in 2:(length(abs_derivs) - 1)
		if valid_mask[i] && valid_mask[i - 1] && valid_mask[i + 1]
			if abs_derivs[i] > abs_derivs[i - 1] && abs_derivs[i] > abs_derivs[i + 1]
				push!(local_maxima, i)
			end
		end
	end

	peak_idx = if isempty(local_maxima)
		valid_indices = findall(valid_mask)
		valid_indices[argmax(abs_derivs[valid_mask])]
	else
		local_maxima[argmax(abs_derivs[local_maxima])]
	end

	left_idx = max(1, peak_idx - 1)
	right_idx = min(length(T_vals), peak_idx + 1)
	a, b = T_vals[left_idx], T_vals[right_idx]
	if a < b
		a, b = b, a
	end

	ϕ = (sqrt(5) - 1) / 2
	iterations = 0

	for iter in 1:max_iter
		iterations = iter
		a - b < tol && break

		c = b + ϕ * (a - b)
		d = a - ϕ * (a - b)

		fc = try abs(evaluate_point(c; order=1, var_idx=var_idx, need_rho=false).d1) catch; NaN end
		fd = try abs(evaluate_point(d; order=1, var_idx=var_idx, need_rho=false).d1) catch; NaN end

		(isnan(fc) || isnan(fd)) && break
		if fc > fd
			b = d
		else
			a = c
		end
	end

	T_crossover = (a + b) / 2
	final_probe = try evaluate_point(T_crossover; order=1, var_idx=var_idx, need_rho=true) catch; nothing end
	final_deriv = final_probe === nothing ? NaN : final_probe.d1
	final_rho = final_probe === nothing ? nothing : final_probe.rho

	all_peaks = isempty(local_maxima) ? Tuple{Float64, Float64}[] : [(T_vals[i], abs_derivs[i]) for i in local_maxima]
	details = Dict{Symbol, Any}(
		:T_range => (T_min, T_max),
		:n_scan => n_scan,
		:variable_index => var_idx,
		:n_peaks_found => length(local_maxima),
		:all_peaks => all_peaks,
		:scan_data => collect(zip(T_vals, derivs)),
	)

	return CrossoverResult(true, T_crossover, final_rho, :peak, final_deriv, iterations, details)
end

function _refine_inflection_brent(eval_d2::Function, T_a::Float64, T_b::Float64, f_a::Float64, f_b::Float64,
		tol::Float64, max_iter::Int)
	if !(isfinite(f_a) && isfinite(f_b)) || f_a * f_b > 0
		return (root=(T_a + T_b) / 2, iterations=0, converged=false, method=:invalid_bracket, final_d2=NaN)
	end

	a, b = T_a, T_b
	fa, fb = f_a, f_b
	if abs(fa) < abs(fb)
		a, b = b, a
		fa, fb = fb, fa
	end

	c = a
	fc = fa
	d = b - a
	e = d
	iterations = 0
	converged = false

	for iter in 1:max_iter
		iterations = iter
		if fb * fc > 0
			c = a
			fc = fa
			d = b - a
			e = d
		end

		if abs(fc) < abs(fb)
			a, b, c = b, c, b
			fa, fb, fc = fb, fc, fb
		end

		tol1 = 2 * eps(Float64) * abs(b) + 0.5 * tol
		xm = 0.5 * (c - b)
		if abs(xm) <= tol1 || fb == 0.0
			converged = true
			break
		end

		if abs(e) >= tol1 && abs(fa) > abs(fb)
			s = fb / fa
			local p::Float64
			local q::Float64
			if a == c
				p = 2 * xm * s
				q = 1 - s
			else
				q1 = fa / fc
				r = fb / fc
				p = s * (2 * xm * q1 * (q1 - r) - (b - a) * (r - 1))
				q = (q1 - 1) * (r - 1) * (s - 1)
			end
			if p > 0
				q = -q
			end
			p = abs(p)
			if 2 * p < min(3 * xm * q - abs(tol1 * q), abs(e * q))
				e = d
				d = p / q
			else
				d = xm
				e = d
			end
		else
			d = xm
			e = d
		end

		a = b
		fa = fb
		if abs(d) > tol1
			b += d
		else
			b += sign(xm) * tol1
		end
		fb = eval_d2(b)
		if !isfinite(fb)
			break
		end
	end

	if converged
		return (root=b, iterations=iterations, converged=true, method=:brent, final_d2=fb)
	end

	lo, hi = min(T_a, T_b), max(T_a, T_b)
	f_lo = eval_d2(lo)
	f_hi = eval_d2(hi)
	if !(isfinite(f_lo) && isfinite(f_hi)) || f_lo * f_hi > 0
		return (root=(lo + hi) / 2, iterations=iterations, converged=false, method=:fallback_failed, final_d2=NaN)
	end

	for k in 1:max_iter
		mid = 0.5 * (lo + hi)
		f_mid = eval_d2(mid)
		iterations += 1
		if !isfinite(f_mid)
			return (root=mid, iterations=iterations, converged=false, method=:fallback_nan, final_d2=NaN)
		end
		if abs(hi - lo) < tol || f_mid == 0.0
			return (root=mid, iterations=iterations, converged=true, method=:bisection_fallback, final_d2=f_mid)
		end
		if f_lo * f_mid < 0
			hi, f_hi = mid, f_mid
		else
			lo, f_lo = mid, f_mid
		end
	end

	root = 0.5 * (lo + hi)
	return (root=root, iterations=iterations, converged=false, method=:bisection_fallback, final_d2=eval_d2(root))
end

function _detect_crossover_inflection(T_min::Float64, T_max::Float64,
		evaluate_point::Function, n_scan::Int, tol::Real, max_iter::Int, var_idx::Int)
	T_vals = collect(range(T_max, T_min; length=n_scan))
	d1_vals, d2_vals = Float64[], Float64[]

	for T in T_vals
		try
			probe = evaluate_point(T; order=2, var_idx=var_idx, need_rho=false)
			push!(d1_vals, probe.d1)
			push!(d2_vals, probe.d2)
		catch
			push!(d1_vals, NaN)
			push!(d2_vals, NaN)
		end
	end

	valid_d1 = filter(!isnan, d1_vals)
	max_abs_d1 = isempty(valid_d1) ? 1.0 : maximum(abs, valid_d1)
	d1_threshold = max_abs_d1 * 0.15

	valid_inflections = Tuple{Int, Float64}[]
	for i in 1:(length(d2_vals) - 1)
		if !isnan(d2_vals[i]) && !isnan(d2_vals[i + 1])
			if d2_vals[i] < 0 && d2_vals[i + 1] > 0
				d1_at_inflection = (abs(d1_vals[i]) + abs(d1_vals[i + 1])) / 2
				d1_at_inflection >= d1_threshold && push!(valid_inflections, (i, d1_at_inflection))
			end
		end
	end

	sign_change_idx = nothing
	if !isempty(valid_inflections)
		_, best_pos = findmax(x -> x[2], valid_inflections)
		sign_change_idx = valid_inflections[best_pos][1]
	end

	if sign_change_idx === nothing
		details = Dict{Symbol, Any}(
			:T_range => (T_min, T_max),
			:n_scan => n_scan,
			:variable_index => var_idx,
			:reason => "no_valid_sign_change",
			:d1_threshold => d1_threshold,
			:refine_method => :none,
			:scan_data => collect(zip(T_vals, d1_vals, d2_vals)),
		)
		return CrossoverResult(false, nothing, nothing, :inflection, nothing, 0, details)
	end

	T1, T2 = T_vals[sign_change_idx], T_vals[sign_change_idx + 1]
	f1, f2 = d2_vals[sign_change_idx], d2_vals[sign_change_idx + 1]
	if T1 > T2
		T1, T2 = T2, T1
		f1, f2 = f2, f1
	end

	eval_d2 = T -> begin
		probe = evaluate_point(T; order=2, var_idx=var_idx, need_rho=false)
		return probe.d2
	end

	refine = _refine_inflection_brent(eval_d2, T1, T2, f1, f2, Float64(tol), max_iter)
	T_crossover = refine.root
	iterations = refine.iterations
	final_probe = try evaluate_point(T_crossover; order=2, var_idx=var_idx, need_rho=true) catch; nothing end
	final_d1 = final_probe === nothing ? NaN : final_probe.d1
	final_d2 = final_probe === nothing ? NaN : final_probe.d2
	final_rho = final_probe === nothing ? nothing : final_probe.rho

	if T_crossover < T1 || T_crossover > T2
		T_crossover = 0.5 * (T1 + T2)
	end

	all_inflections = Tuple{Float64, Float64, Bool}[]
	for i in 1:(length(d2_vals) - 1)
		if !isnan(d2_vals[i]) && !isnan(d2_vals[i + 1])
			if d2_vals[i] < 0 && d2_vals[i + 1] > 0
				d1_at_inflection = (abs(d1_vals[i]) + abs(d1_vals[i + 1])) / 2
				is_valid = d1_at_inflection >= d1_threshold
				push!(all_inflections, (T_vals[i], T_vals[i + 1], is_valid))
			end
		end
	end

	details = Dict{Symbol, Any}(
		:T_range => (T_min, T_max),
		:n_scan => n_scan,
		:variable_index => var_idx,
		:final_d2 => final_d2,
		:d1_threshold => d1_threshold,
		:n_inflections_found => length(all_inflections),
		:all_inflections => all_inflections,
		:refine_method => refine.method,
		:refine_converged => refine.converged,
		:scan_data => collect(zip(T_vals, d1_vals, d2_vals)),
	)

	return CrossoverResult(true, T_crossover, final_rho, :inflection, final_d1, iterations, details)
end

function detect_crossover(μ_fm::Real, T_range::Tuple{Real, Real};
		method::Symbol=:peak,
		variable::Symbol=:phi_u,
		xi::Real=0.0,
		n_scan::Int=20,
		tol::Real=1e-4,
		max_iter::Int=20,
		p_num::Int=24,
		t_num::Int=12,
		model_kind::Symbol=:PNJL,
		solver_backend::Symbol=:models)
	T_min, T_max = Float64.(T_range)
	T_min < T_max || return CrossoverResult(method=method)

	var_idx = variable == :phi_u ? 1 : (variable == :Phi ? 4 : 1)
	evaluate_point = _build_crossover_evaluator(Float64(μ_fm), xi, p_num, t_num, model_kind, solver_backend)

	result = if method == :peak
		_detect_crossover_peak(T_min, T_max, evaluate_point, n_scan, tol, max_iter, var_idx)
	elseif method == :inflection
		_detect_crossover_inflection(T_min, T_max, evaluate_point, n_scan, tol, max_iter, var_idx)
	else
		error("Unknown method: $method. Use :peak or :inflection")
	end

	return result
end

function scan_crossover_line(mu_range::Tuple{Real, Real, Int}, T_range::Tuple{Real, Real};
		method::Symbol=:peak,
		variable::Symbol=:phi_u,
		xi::Real=0.0,
		model_kind::Symbol=:PNJL,
		solver_backend::Symbol=:models,
		kwargs...)
	μ_min, μ_max, n_mu = mu_range
	μ_vals = range(Float64(μ_min), Float64(μ_max); length=n_mu)

	results = NamedTuple{(:mu_fm, :T_crossover_fm, :rho, :converged, :derivative),
		Tuple{Float64, Union{Nothing, Float64}, Union{Nothing, Float64}, Bool, Union{Nothing, Float64}}}[]

	for μ in μ_vals
		result = detect_crossover(μ, T_range;
			method=method,
			variable=variable,
			xi=xi,
			model_kind=model_kind,
			solver_backend=solver_backend,
			kwargs...)
		push!(results, (
			mu_fm=μ,
			T_crossover_fm=result.T_crossover,
			rho=result.rho,
			converged=result.found,
			derivative=result.derivative_value,
		))
	end
	return results
end

function build_crossover_line(; mu_max_MeV::Real,
		T_min_MeV::Real,
		T_max_MeV::Real,
		xi::Real=0.0,
		n_mu::Int=12,
		method::Symbol=:peak,
		variable::Symbol=:phi_u,
		model_kind::Symbol=:PNJL,
		solver_backend::Symbol=:models)
	mu_max_MeV <= 0 && return NamedTuple[]
	T_min_MeV < T_max_MeV || return NamedTuple[]

	μ_max_fm = Float64(mu_max_MeV) / HBARC_MEV_FM
	T_min_fm = Float64(T_min_MeV) / HBARC_MEV_FM
	T_max_fm = Float64(T_max_MeV) / HBARC_MEV_FM

	raw = scan_crossover_line((0.0, μ_max_fm, max(3, n_mu)), (T_min_fm, T_max_fm);
		method=method,
		variable=variable,
		xi=xi,
		model_kind=model_kind,
		solver_backend=solver_backend)

	rows = NamedTuple[]
	for row in raw
		T_mev = row.T_crossover_fm === nothing ? NaN : Float64(row.T_crossover_fm) * HBARC_MEV_FM
		mu_mev = Float64(row.mu_fm) * HBARC_MEV_FM
		push!(rows, (
			mu_MeV=mu_mev,
			T_crossover_MeV=T_mev,
			rho=row.rho === nothing ? NaN : Float64(row.rho),
			method=String(method),
			converged=Bool(row.converged),
			derivative=row.derivative === nothing ? NaN : Float64(row.derivative),
			variable=String(variable),
		))
	end
	return rows
end
