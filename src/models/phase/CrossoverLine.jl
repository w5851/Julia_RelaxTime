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

function _compute_rho_at_crossover(T_fm::Float64, μ_fm::Float64, xi::Real, p_num::Int, t_num::Int, model_kind::Symbol)
	try
		deriv = solve_pnjl_with_derivatives(T_fm, μ_fm;
			order=1,
			xi=xi,
			p_num=p_num,
			t_num=t_num,
			thermo_backend=:models,
			solver_backend=:models)
		model = create_model(model_kind)
		_, rho_norm, _, _ = model_thermo(model, deriv.x, normalize_mu_vec(μ_fm), T_fm;
			xi=xi, p_num=p_num, t_num=t_num)
		return Float64(rho_norm)
	catch
		return nothing
	end
end

function _detect_crossover_peak(T_min::Float64, T_max::Float64,
		compute_derivatives::Function, n_scan::Int, tol::Real, max_iter::Int, var_idx::Int)
	T_vals = collect(range(T_max, T_min; length=n_scan))
	abs_derivs = Float64[]
	derivs = Float64[]

	for T in T_vals
		try
			dφ_dT, _ = compute_derivatives(T)
			push!(derivs, dφ_dT)
			push!(abs_derivs, abs(dφ_dT))
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

		fc = try dφ_dT, _ = compute_derivatives(c); abs(dφ_dT) catch; NaN end
		fd = try dφ_dT, _ = compute_derivatives(d); abs(dφ_dT) catch; NaN end

		(isnan(fc) || isnan(fd)) && break
		if fc > fd
			b = d
		else
			a = c
		end
	end

	T_crossover = (a + b) / 2
	final_deriv, _ = try compute_derivatives(T_crossover) catch; (NaN, NaN) end

	all_peaks = isempty(local_maxima) ? Tuple{Float64, Float64}[] : [(T_vals[i], abs_derivs[i]) for i in local_maxima]
	details = Dict{Symbol, Any}(
		:T_range => (T_min, T_max),
		:n_scan => n_scan,
		:variable_index => var_idx,
		:n_peaks_found => length(local_maxima),
		:all_peaks => all_peaks,
		:scan_data => collect(zip(T_vals, derivs)),
	)

	return CrossoverResult(true, T_crossover, nothing, :peak, final_deriv, iterations, details)
end

function _detect_crossover_inflection(T_min::Float64, T_max::Float64,
		compute_derivatives::Function, n_scan::Int, tol::Real, max_iter::Int, var_idx::Int)
	T_vals = collect(range(T_max, T_min; length=n_scan))
	d1_vals, d2_vals = Float64[], Float64[]

	for T in T_vals
		try
			dφ_dT, d2φ_dT2 = compute_derivatives(T)
			push!(d1_vals, dφ_dT)
			push!(d2_vals, d2φ_dT2)
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
			:scan_data => collect(zip(T_vals, d1_vals, d2_vals)),
		)
		return CrossoverResult(false, nothing, nothing, :inflection, nothing, 0, details)
	end

	a, b = T_vals[sign_change_idx], T_vals[sign_change_idx + 1]
	if a < b
		a, b = b, a
	end
	fa, fb = d2_vals[sign_change_idx], d2_vals[sign_change_idx + 1]

	iterations = 0
	for iter in 1:max_iter
		iterations = iter
		a - b < tol && break
		mid = (a + b) / 2
		_, f_mid = try compute_derivatives(mid) catch; (NaN, NaN) end
		isnan(f_mid) && break
		if fa * f_mid < 0
			b, fb = mid, f_mid
		else
			a, fa = mid, f_mid
		end
	end

	T_crossover = (a + b) / 2
	final_d1, final_d2 = try compute_derivatives(T_crossover) catch; (NaN, NaN) end

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
		:scan_data => collect(zip(T_vals, d1_vals, d2_vals)),
	)

	return CrossoverResult(true, T_crossover, nothing, :inflection, final_d1, iterations, details)
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

	function compute_derivatives(T::Float64)
		result = solve_pnjl_with_derivatives(T, Float64(μ_fm);
			order=2,
			xi=xi,
			p_num=p_num,
			t_num=t_num,
			thermo_backend=:models,
			solver_backend=solver_backend)
		return result.dx_dT[var_idx], result.d2x_dT2[var_idx]
	end

	result = if method == :peak
		_detect_crossover_peak(T_min, T_max, compute_derivatives, n_scan, tol, max_iter, var_idx)
	elseif method == :inflection
		_detect_crossover_inflection(T_min, T_max, compute_derivatives, n_scan, tol, max_iter, var_idx)
	else
		error("Unknown method: $method. Use :peak or :inflection")
	end

	if result.found && result.T_crossover !== nothing
		rho = _compute_rho_at_crossover(result.T_crossover, Float64(μ_fm), xi, p_num, t_num, model_kind)
		return CrossoverResult(result.found, result.T_crossover, rho, result.method,
			result.derivative_value, result.iterations, result.details)
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
		method::Symbol=:inflection,
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
