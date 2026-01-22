---
title: RelaxTime Fortran 对照/诊断脚本归档
archived: true
original: scripts/relaxtime/* (临时对照/实验脚本)
archived_date: 2026-01-22
---

本文件集中归档与 Fortran 参考实现对照相关的临时脚本。它们通常依赖外部 Fortran 工作区或特定 CSV 中间产物，因此不作为仓库默认 workflow 的一部分。

说明：

- 这些脚本已从 `scripts/relaxtime` 清理。
- 如需复现，可从本归档文档中把对应脚本代码复制回 `.jl` 文件再运行。

提示：

- 本文内嵌的脚本源码与注释中可能仍包含历史路径形如 `scripts/relaxtime/*.jl` 的用法示例；由于脚本文件已从仓库移除，复现时请将该路径替换为你保存的本地临时脚本路径（例如 `D:/tmp/<script>.jl`）。

## 对照入口（meson.dat）

脚本：`compare_meson_masses_with_fortran.jl`

用途：抽样对比 Julia 介子质量与 Fortran `quark_phase/meson.dat`。

关键点：

- 内置 `meson.dat` 列映射（T, π, K, η′, η, σπ, σK, σ′, σ）。
- 混合通道支持 swap 匹配与不变量（`sum(m^2)`, `prod(m^2)`）诊断。

```julia
#+#+#+#+
# Compare Julia meson masses vs Fortran reference output.
#
# NOTE: This script depends on an external Fortran workspace (not tracked by this repo).
# It is archived for reference/reproducibility and is not part of the default workflow.
#
# Usage (PowerShell):
#   julia --project=. .\scripts\relaxtime\archived\compare_meson_masses_with_fortran.jl \
#     --fortran-meson-file "D:/Desktop/fortran代码/输运系数/RelaxTime/PNJL-mott-mu_T/PNJL-mu-T/quark_phase/meson.dat" \
#     --t-list 150,180,200,210,250,300

include("../../../src/Constants_PNJL.jl")
include("../../../src/pnjl/workflows/MesonMassWorkflow.jl")

using .Constants_PNJL: ħc_MeV_fm
using .MesonMassWorkflow: solve_gap_and_meson_point
using Printf

function _parse_args(argv)
	args = Dict{String,String}()
	i = 1
	while i <= length(argv)
		a = argv[i]
		if startswith(a, "--")
			key = a[3:end]
			if i == length(argv)
				args[key] = "true"
			else
				nxt = argv[i+1]
				if startswith(nxt, "--")
					args[key] = "true"
				else
					args[key] = nxt
					i += 1
				end
			end
		end
		i += 1
	end
	return args
end

function _parse_t_list(s::AbstractString)
	isempty(strip(s)) && return Float64[]
	return [parse(Float64, strip(x)) for x in split(s, ',') if !isempty(strip(x))]
end

function _parse_symbol_list(s::AbstractString)
	isempty(strip(s)) && return Symbol[]
	out = Symbol[]
	for x in split(s, ',')
		t = strip(x)
		isempty(t) && continue
		push!(out, Symbol(t))
	end
	return out
end

"""Parse Fortran quark_phase/meson.dat.

Format observed (9 floats per line, all in MeV):

	T, M_pi, M_K, M_eta_prime, M_eta, M_sigma_pi, M_sigma_K, M_sigma_prime, M_sigma

Note: The scalar channels are easy to mislabel across codes; this mapping is
chosen to match the Fortran reference outputs used in this repo.
"""
function read_fortran_meson_dat(path::AbstractString)
	rows = Vector{NamedTuple}()
	open(path, "r") do io
		for ln in eachline(io)
			s = strip(ln)
			isempty(s) && continue
			parts = split(s)
			length(parts) < 9 && continue
			vals = parse.(Float64, parts[1:9])
			T = vals[1]
			push!(rows, (
				T_MeV=T,
				M_pi=vals[2],
				M_K=vals[3],
				M_eta_prime=vals[4],
				M_eta=vals[5],
				M_sigma_pi=vals[6],
				M_sigma_K=vals[7],
				M_sigma_prime=vals[8],
				M_sigma=vals[9],
			))
		end
	end
	return rows
end

function _relerr(j::Real, f::Real)
	if isfinite(j) && isfinite(f) && f != 0
		return (j - f) / f
	end
	return NaN
end

function _print_one(name::AbstractString, jf::Real, ff::Real)
	if isfinite(ff) && ff != 0 && isfinite(jf)
		rel = (jf - ff) / ff
		println(rpad(name, 9), " Julia=", lpad(_fmt(jf), 10), "  Fortran=", lpad(_fmt(ff), 10), "  rel=", @sprintf("%+.2f%%", 100 * rel))
	else
		println(rpad(name, 9), " Julia=", lpad(_fmt(jf), 10), "  Fortran=", lpad(_fmt(ff), 10))
	end
	return nothing
end

"""Read a scan_csv_v1 file and return the first row matching (T_MeV, muB_MeV, xi).

Assumes scan output columns like:
  T_MeV, muB_MeV, xi, M_pi, M_K, M_eta, ... (masses in fm^-1)
"""
function read_scan_csv_row(path::AbstractString; T_MeV::Real, muB_MeV::Real, xi::Real)
	header = nothing
	open(path, "r") do io
		for ln in eachline(io)
			s = strip(ln)
			isempty(s) && continue
			startswith(s, "#") && continue
			header = split(s, ',')
			break
		end
	end
	header === nothing && error("no CSV header found in scan file: $path")
	col = Dict{String,Int}()
	for (i, name) in enumerate(header)
		col[strip(name)] = i
	end

	function getf(fields, key)
		i = get(col, key, 0)
		i == 0 && return NaN
		v = strip(fields[i])
		isempty(v) && return NaN
		return try
			parse(Float64, v)
		catch
			NaN
		end
	end

	open(path, "r") do io
		for ln in eachline(io)
			s = strip(ln)
			isempty(s) && continue
			startswith(s, "#") && continue
			startswith(s, "T_MeV,") && continue
			fields = split(s, ',')

			T0 = getf(fields, "T_MeV")
			muB0 = getf(fields, "muB_MeV")
			xi0 = getf(fields, "xi")
			if isfinite(T0) && isfinite(muB0) && isfinite(xi0) && abs(T0 - T_MeV) < 1e-9 && abs(muB0 - muB_MeV) < 1e-9 && abs(xi0 - xi) < 1e-12
				return (fields=fields, col=col)
			end
		end
	end
	return nothing
end

function _nearest_by_T(rows, T_MeV::Real)
	best = nothing
	best_d = Inf
	for r in rows
		d = abs(r.T_MeV - T_MeV)
		if d < best_d
			best = r
			best_d = d
		end
	end
	return best
end

function _nearest_index_by_T(rows, T_MeV::Real)
	best_i = 0
	best_d = Inf
	for (i, r) in enumerate(rows)
		d = abs(r.T_MeV - T_MeV)
		if d < best_d
			best_i = i
			best_d = d
		end
	end
	return best_i
end

function _print_fortran_neighbors(rows, T_MeV::Real; window::Int=2)
	window <= 0 && return
	i0 = _nearest_index_by_T(rows, T_MeV)
	i0 == 0 && return
	lo = max(1, i0 - window)
	hi = min(length(rows), i0 + window)
	println("[fortran] neighbor rows (T, eta', eta)")
	for i in lo:hi
		r = rows[i]
		println(@sprintf("  T=%7.1f  eta'=%9.3f  eta=%9.3f", r.T_MeV, r.M_eta_prime, r.M_eta))
	end
	return nothing
end

function _mix_invariants(m1::Real, m2::Real)
	if !(isfinite(m1) && isfinite(m2))
		return (sum2=NaN, prod2=NaN)
	end
	a = Float64(m1)^2
	b = Float64(m2)^2
	return (sum2=a + b, prod2=a * b)
end

function _print_mix_invariants(label::AbstractString, j1::Real, j2::Real, f1::Real, f2::Real)
	ji = _mix_invariants(j1, j2)
	fi = _mix_invariants(f1, f2)
	if isfinite(ji.sum2) && isfinite(fi.sum2) && fi.sum2 != 0
		println(@sprintf("[%s] sum(m^2): Julia=%10.3e  Fortran=%10.3e  rel=%+.2f%%", label, ji.sum2, fi.sum2, 100 * (ji.sum2 - fi.sum2) / fi.sum2))
	else
		println(@sprintf("[%s] sum(m^2): Julia=%10.3e  Fortran=%10.3e", label, ji.sum2, fi.sum2))
	end
	if isfinite(ji.prod2) && isfinite(fi.prod2) && fi.prod2 != 0
		println(@sprintf("[%s] prod(m^2): Julia=%10.3e  Fortran=%10.3e  rel=%+.2f%%", label, ji.prod2, fi.prod2, 100 * (ji.prod2 - fi.prod2) / fi.prod2))
	else
		println(@sprintf("[%s] prod(m^2): Julia=%10.3e  Fortran=%10.3e", label, ji.prod2, fi.prod2))
	end
	return nothing
end

function _fmt(x)
	return isfinite(x) ? @sprintf("%.3f", x) : "NaN"
end

function main(argv)
	args = _parse_args(argv)

	fortran_file = get(args, "fortran-meson-file", "")
	isempty(fortran_file) && error("missing --fortran-meson-file")

	scan_file = get(args, "julia-scan-file", "")

	t_list = _parse_t_list(get(args, "t-list", "150,180,200,210,250,300"))
	isempty(t_list) && error("empty --t-list")

	xi = parse(Float64, get(args, "xi", "0.0"))
	muB_MeV = parse(Float64, get(args, "muB", "0.0"))
	mu_MeV = muB_MeV / 3.0

	meson_list = _parse_symbol_list(get(args, "mesons", "pi,K,eta,eta_prime,sigma_pi,sigma_K,sigma,sigma_prime"))
	mesons = Tuple(meson_list)

	rows = read_fortran_meson_dat(fortran_file)
	isempty(rows) && error("no data parsed from fortran file: $fortran_file")

	neighbor_window = parse(Int, get(args, "neighbor-window", "2"))

	println("Compare (muB=$(muB_MeV) MeV, xi=$(xi))")
	println("Fortran: $fortran_file")
	println("Columns: T, pi, K, eta', eta, sigma_pi, sigma_K, sigma', sigma (MeV)")

	for T_MeV in t_list
		ref = _nearest_by_T(rows, T_MeV)
		ref === nothing && continue
		abs(ref.T_MeV - T_MeV) > 1e-6 && println("[warn] T=$(T_MeV) MeV: using nearest Fortran row at T=$(ref.T_MeV) MeV")

		_print_fortran_neighbors(rows, ref.T_MeV; window=neighbor_window)

		j = if !isempty(scan_file)
			row = read_scan_csv_row(scan_file; T_MeV=ref.T_MeV, muB_MeV=muB_MeV, xi=xi)
			row === nothing && error("no matching row in scan file for T=$(ref.T_MeV), muB=$(muB_MeV), xi=$(xi): $scan_file")
			fields = row.fields
			col = row.col
			function getm(key)
				i = get(col, key, 0)
				i == 0 && return NaN
				v = strip(fields[i])
				isempty(v) && return NaN
				mf = try parse(Float64, v) catch; NaN end
				return isfinite(mf) ? mf * ħc_MeV_fm : NaN
			end
			(
				M_pi=getm("M_pi"),
				M_K=getm("M_K"),
				M_eta=getm("M_eta"),
				M_eta_prime=getm("M_eta_prime"),
				M_sigma_pi=getm("M_sigma_pi"),
				M_sigma_K=getm("M_sigma_K"),
				M_sigma=getm("M_sigma"),
				M_sigma_prime=getm("M_sigma_prime"),
			)
		else
			p_num = parse(Int, get(args, "p-num", "16"))
			t_num = parse(Int, get(args, "t-num", "8"))

			T_fm = T_MeV / ħc_MeV_fm
			mu_fm = mu_MeV / ħc_MeV_fm

			res = solve_gap_and_meson_point(
				T_fm,
				mu_fm;
				xi=xi,
				mesons=mesons,
				p_num=p_num,
				t_num=t_num,
				mass_kwargs=(;),
			)

			function jmass(m)
				r = get(res.meson_results, m, nothing)
				r === nothing && return NaN
				return Float64(r.mass) * ħc_MeV_fm
			end

			(
				M_pi=jmass(:pi),
				M_K=jmass(:K),
				M_eta=jmass(:eta),
				M_eta_prime=jmass(:eta_prime),
				M_sigma_pi=jmass(:sigma_pi),
				M_sigma_K=jmass(:sigma_K),
				M_sigma=jmass(:sigma),
				M_sigma_prime=jmass(:sigma_prime),
			)
		end

		f = (
			M_pi=ref.M_pi,
			M_K=ref.M_K,
			M_eta=ref.M_eta,
			M_eta_prime=ref.M_eta_prime,
			M_sigma_pi=ref.M_sigma_pi,
			M_sigma_K=ref.M_sigma_K,
			M_sigma=ref.M_sigma,
			M_sigma_prime=ref.M_sigma_prime,
		)

		println("\nT=$(ref.T_MeV) MeV")

		:pi in mesons && _print_one("pi", j.M_pi, f.M_pi)
		:K in mesons && _print_one("K", j.M_K, f.M_K)

		# Mixing eigenstates may swap ordering across implementations.
		if (:eta in mesons) || (:eta_prime in mesons)
			# Try both assignments and pick the one with smaller total mismatch.
			e1 = abs(_relerr(j.M_eta, f.M_eta)) + abs(_relerr(j.M_eta_prime, f.M_eta_prime))
			e2 = abs(_relerr(j.M_eta, f.M_eta_prime)) + abs(_relerr(j.M_eta_prime, f.M_eta))
			if e2 < e1
				:eta in mesons && _print_one("eta", j.M_eta, f.M_eta_prime)
				:eta_prime in mesons && _print_one("eta'", j.M_eta_prime, f.M_eta)
				println("[note] eta/eta' matched with swap")
			else
				:eta in mesons && _print_one("eta", j.M_eta, f.M_eta)
				:eta_prime in mesons && _print_one("eta'", j.M_eta_prime, f.M_eta_prime)
			end

			_print_mix_invariants("eta", j.M_eta, j.M_eta_prime, f.M_eta, f.M_eta_prime)
		end

		:sigma_pi in mesons && _print_one("sigma_pi", j.M_sigma_pi, f.M_sigma_pi)
		:sigma_K in mesons && _print_one("sigma_K", j.M_sigma_K, f.M_sigma_K)

		if (:sigma in mesons) || (:sigma_prime in mesons)
			s1 = abs(_relerr(j.M_sigma, f.M_sigma)) + abs(_relerr(j.M_sigma_prime, f.M_sigma_prime))
			s2 = abs(_relerr(j.M_sigma, f.M_sigma_prime)) + abs(_relerr(j.M_sigma_prime, f.M_sigma))
			if s2 < s1
				:sigma in mesons && _print_one("sigma", j.M_sigma, f.M_sigma_prime)
				:sigma_prime in mesons && _print_one("sigma'", j.M_sigma_prime, f.M_sigma)
				println("[note] sigma/sigma' matched with swap")
			else
				:sigma in mesons && _print_one("sigma", j.M_sigma, f.M_sigma)
				:sigma_prime in mesons && _print_one("sigma'", j.M_sigma_prime, f.M_sigma_prime)
			end

			_print_mix_invariants("sigma", j.M_sigma, j.M_sigma_prime, f.M_sigma, f.M_sigma_prime)
		end
	end
end

main(ARGS)
```

## A/K 系数对比

脚本：`compare_ak_fortran.jl`

用途：对比 Fortran 输出的 A、G、K 与 Julia 计算结果，生成差值 CSV。

```julia
# 对比 Fortran 与 Julia 的 A / K 系数。
#
# 输入 CSV：data/outputs/results/ak_fortran_output.csv
# 输出 CSV：data/outputs/results/ak_compare.csv

using CSV

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "OneLoopIntegrals.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))
include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))

using .GaussLegendre: gauleg
using .OneLoopIntegrals: A
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .Constants_PNJL: G_fm2, K_fm5

function main()
	input_path = joinpath(PROJECT_ROOT, "data", "outputs", "results", "ak_fortran_output.csv")
	output_path = joinpath(PROJECT_ROOT, "data", "outputs", "results", "ak_compare.csv")

	rows = CSV.File(input_path)
	results = NamedTuple[]

	nodes_p, weights_p = gauleg(0.0, 20.0, 16)

	for row in rows
		T = row.T
		mu_u = row.mu_u
		mu_s = row.mu_s
		m_u = row.m_u
		m_s = row.m_s
		Phi = row.Phi
		Phibar = row.Phibar

		A_u_j = A(m_u, mu_u, T, Phi, Phibar, nodes_p, weights_p)
		A_s_j = A(m_s, mu_s, T, Phi, Phibar, nodes_p, weights_p)
		G_u_j = calculate_G_from_A(A_u_j, m_u)
		G_s_j = calculate_G_from_A(A_s_j, m_s)
		K_j = calculate_effective_couplings(G_fm2, K_fm5, G_u_j, G_s_j)

		push!(results, (
			T=row.T,
			mu_u=row.mu_u,
			mu_d=row.mu_d,
			mu_s=row.mu_s,
			m_u=row.m_u,
			m_d=row.m_d,
			m_s=row.m_s,
			Phi=row.Phi,
			Phibar=row.Phibar,
			A_u_fortran=row.A_u,
			A_s_fortran=row.A_s,
			A_u_julia=A_u_j,
			A_s_julia=A_s_j,
			dA_u=A_u_j - row.A_u,
			dA_s=A_s_j - row.A_s,
			G_u_fortran=row.G_u,
			G_s_fortran=row.G_s,
			G_u_julia=G_u_j,
			G_s_julia=G_s_j,
			dG_u=G_u_j - row.G_u,
			dG_s=G_s_j - row.G_s,
			K123_plus_fortran=row.K123_plus,
			K4567_plus_fortran=row.K4567_plus,
			K123_plus_julia=K_j.K123_plus,
			K4567_plus_julia=K_j.K4567_plus,
			dK123=K_j.K123_plus - row.K123_plus,
			dK4567=K_j.K4567_plus - row.K4567_plus,
		))
	end

	CSV.write(output_path, results)
end

main()
```

## B0 对比（Julia vs Fortran 输出）

脚本：`compare_b0_fortran.jl`

用途：对比 Julia 的 `OneLoopIntegrals.B0` 与 Fortran 输出的 B0（实部/虚部）。

```julia
"""
对比 Julia 的 B0 与 Fortran 输出。

输入 CSV 需包含列：
λ, k, m1, m2, μ1, μ2, T, Φ, Φbar, B0_re_fortran, B0_im_fortran

输出 CSV 追加计算结果与差值：
B0_re_julia, B0_im_julia, dB0_re, dB0_im

用法：
julia scripts/relaxtime/compare_b0_fortran.jl input.csv [output.csv]
"""

using CSV

include("../../src/relaxtime/OneLoopIntegrals.jl")

using .OneLoopIntegrals: B0

function _require(row, key::Symbol, alt::Union{Nothing,Symbol}=nothing)
	if hasproperty(row, key)
		return getproperty(row, key)
	elseif alt !== nothing && hasproperty(row, alt)
		return getproperty(row, alt)
	end
	error("Missing column: $(key)")
end

function _row_result(row)
	λ = _require(row, :λ, :lam)
	k = _require(row, :k)
	m1 = _require(row, :m1)
	m2 = _require(row, :m2)
	μ1 = _require(row, :μ1, :mu1)
	μ2 = _require(row, :μ2, :mu2)
	T = _require(row, :T)
	Φ = _require(row, :Phi)
	Φbar = _require(row, :Phibar)

	B0_re, B0_im = B0(λ, k, m1, μ1, m2, μ2, T; Φ=Φ, Φbar=Φbar)

	return (
		B0_re_julia = B0_re,
		B0_im_julia = B0_im,
		dB0_re = B0_re - _require(row, :B0_re_fortran),
		dB0_im = B0_im - _require(row, :B0_im_fortran),
	)
end

function main()
	if length(ARGS) < 1
		error("Usage: julia compare_b0_fortran.jl input.csv [output.csv]")
	end
	input_path = ARGS[1]
	output_path = length(ARGS) >= 2 ? ARGS[2] : nothing

	rows = CSV.File(input_path)
	results = NamedTuple[]
	for row in rows
		result = _row_result(row)
		push!(results, result)
	end

	if output_path === nothing
		for r in results
			println(r)
		end
	else
		CSV.write(output_path, results)
	end
end

main()
```

## B0 交叉验证（QuadGK 主值积分）

脚本：`compare_b0_quadgk.jl`

用途：使用 QuadGK 的分段主值积分独立计算 B0，用于判断 Julia/Fortran B0 哪个更接近数值收敛极限。

```julia
# 使用独立数值方法（QuadGK 分段主值积分）交叉验证 B0。
#
# 输入 CSV 需包含列：
# lam,k,m1,m2,mu1,mu2,T,Phi,Phibar,B0_re_fortran,B0_im_fortran
#
# 输出 CSV 追加：
# B0_re_quadgk,B0_im_quadgk,dB0_re_julia_quadgk,dB0_im_julia_quadgk,dB0_re_fortran_quadgk,dB0_im_fortran_quadgk
#
# 用法：
# julia scripts/relaxtime/compare_b0_quadgk.jl input.csv [output.csv]

using CSV
using QuadGK

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "relaxtime", "OneLoopIntegrals.jl"))

using .OneLoopIntegrals: B0, EPS_K, EPS_SEGMENT, PV_GAP_REL,
	energy_cutoff, internal_momentum, distribution_value_b0, distribution_integral_b0,
	singularity_k_zero, singularity_k_positive,
	real_integrand_k_zero, real_integrand_k_positive

const QTOL = parse(Float64, get(ENV, "B0_QUADGK_RTOL", "1e-9"))

function _require(row, key::Symbol, alt::Union{Nothing,Symbol}=nothing)
	if hasproperty(row, key)
		return getproperty(row, key)
	elseif alt !== nothing && hasproperty(row, alt)
		return getproperty(row, alt)
	end
	error("Missing column: $(key)")
end

function _quadgk_split(integrand, a::Float64, b::Float64, split_points::Vector{Float64})
	pts = sort([a; split_points; b])
	total = 0.0
	for i in 1:(length(pts)-1)
		x1 = pts[i]
		x2 = pts[i+1]
		if x2 <= x1
			continue
		end
		val, _ = quadgk(integrand, x1, x2; rtol=QTOL)
		total += val
	end
	return total
end

function tilde_B0_k_zero_quadgk(sign_flag::Symbol, λ::Float64, m::Float64, m_prime::Float64,
	μ::Float64, T::Float64, Φ::Float64, Φbar::Float64)
	m_pos = max(m, 0.0)
	m_prime_pos = max(m_prime, 0.0)
	Emin = m_pos
	Emax = energy_cutoff(m_pos)
	denominator_term = (λ ^ 2 + m_pos ^ 2 - m_prime_pos ^ 2) / 2.0
	singularity = singularity_k_zero(λ, Emin, Emax, denominator_term)

	integrand(E) = real_integrand_k_zero(sign_flag, λ, m_pos, denominator_term, μ, T, Φ, Φbar, E)

	real_part = 0.0
	if isempty(singularity)
		real_part, _ = quadgk(integrand, Emin, Emax; rtol=QTOL)
	else
		E0 = singularity[1]
		gap = max(PV_GAP_REL * (Emax - Emin), 1e-8)
		left_end = max(Emin, E0 - gap)
		right_start = min(Emax, E0 + gap)
		if left_end > Emin
			val, _ = quadgk(integrand, Emin, left_end; rtol=QTOL)
			real_part += val
		end
		if right_start < Emax
			val, _ = quadgk(integrand, right_start, Emax; rtol=QTOL)
			real_part += val
		end
	end

	# 解析虚部（残数项）
	imag_part = 0.0
	if !isempty(singularity)
		E0 = singularity[1]
		p0 = internal_momentum(E0, m_pos)
		imag_part = 2.0 * π * p0 * distribution_value_b0(sign_flag, E0, μ, T, Φ, Φbar)
	end

	return real_part * 2.0, imag_part / λ
end

function tilde_B0_k_positive_quadgk(sign_flag::Symbol, λ::Float64, k::Float64, m::Float64, m_prime::Float64,
	μ::Float64, T::Float64, Φ::Float64, Φbar::Float64)
	m_pos = max(m, 0.0)
	m_prime_pos = max(m_prime, 0.0)
	Emin = m_pos
	Emax = energy_cutoff(m_pos)

	integrand(E) = real_integrand_k_positive(sign_flag, λ, k, m_pos, m_prime_pos, μ, T, Φ, Φbar, E)
	intervals, _ = singularity_k_positive(λ, k, m_pos, m_prime_pos, Emin, Emax)

	split_points = Float64[]
	for (a, b) in intervals
		push!(split_points, a)
		push!(split_points, b)
	end

	real_part = _quadgk_split(integrand, Emin, Emax, split_points)

	# 独立数值计算虚部（用分布函数积分）
	imag_part = 0.0
	if !isempty(intervals)
		for (E1, E2) in intervals
			val, _ = quadgk(E -> distribution_value_b0(sign_flag, E, μ, T, Φ, Φbar), E1, E2; rtol=QTOL)
			imag_part += val
		end
		imag_part *= π * sign(λ)
	end

	return real_part / k, imag_part / k
end

function tilde_B0_quadgk(sign_flag::Symbol, λ::Float64, k::Float64, m::Float64, m_prime::Float64,
	μ::Float64, T::Float64, Φ::Float64, Φbar::Float64)
	if abs(k) < EPS_K
		return tilde_B0_k_zero_quadgk(sign_flag, λ, m, m_prime, μ, T, Φ, Φbar)
	else
		return tilde_B0_k_positive_quadgk(sign_flag, λ, k, m, m_prime, μ, T, Φ, Φbar)
	end
end

function B0_quadgk(λ::Float64, k::Float64, m1::Float64, μ1::Float64, m2::Float64, μ2::Float64, T::Float64;
	Φ::Float64=0.0, Φbar::Float64=0.0)
	term1 = tilde_B0_quadgk(:plus, -λ, k, m1, m2, μ1, T, Φ, Φbar)
	term2 = tilde_B0_quadgk(:minus, λ, k, m1, m2, μ1, T, Φ, Φbar)
	term3 = tilde_B0_quadgk(:plus, λ, k, m2, m1, μ2, T, Φ, Φbar)
	term4 = tilde_B0_quadgk(:minus, -λ, k, m2, m1, μ2, T, Φ, Φbar)

	real_part = term1[1] - term2[1] + term3[1] - term4[1]
	imag_part = term1[2] - term2[2] + term3[2] - term4[2]
	return real_part, imag_part
end

function _row_result(row)
	lam = _require(row, :lam, :λ)
	k = row.k
	m1 = row.m1
	m2 = row.m2
	mu1 = _require(row, :mu1, :μ1)
	mu2 = _require(row, :mu2, :μ2)
	T = row.T
	Phi = row.Phi
	Phibar = row.Phibar

	b0_j_re, b0_j_im = B0(lam, k, m1, mu1, m2, mu2, T; Φ=Phi, Φbar=Phibar)
	b0_q_re, b0_q_im = B0_quadgk(lam, k, m1, mu1, m2, mu2, T; Φ=Phi, Φbar=Phibar)

	return (
		B0_re_quadgk = b0_q_re,
		B0_im_quadgk = b0_q_im,
		dB0_re_julia_quadgk = b0_j_re - b0_q_re,
		dB0_im_julia_quadgk = b0_j_im - b0_q_im,
		dB0_re_fortran_quadgk = b0_q_re - row.B0_re_fortran,
		dB0_im_fortran_quadgk = b0_q_im - row.B0_im_fortran,
	)
end

function main()
	if length(ARGS) < 1
		error("Usage: julia compare_b0_quadgk.jl input.csv [output.csv]")
	end
	input_path = ARGS[1]
	output_path = length(ARGS) >= 2 ? ARGS[2] : nothing

	rows = CSV.File(input_path)
	results = NamedTuple[]
	for row in rows
		result = _row_result(row)
		push!(results, result)
	end

	if output_path === nothing
		for r in results
			println(r)
		end
	else
		CSV.write(output_path, results)
	end
end

main()
```

## MesonMass/MottTransition 对比（π/K，基于 CSV 输入）

脚本：`compare_meson_mass_fortran.jl`

用途：给定 Fortran 导出的输入点（T、μ、m、Φ 等）与 π/K 的 (M, Γ) 初值，对比 Julia 求解结果与残差。

```julia
"""
对比 MesonMass/MottTransition 结果与 Fortran 输出。

输入 CSV 需包含列：
T_fm, mu_u, mu_d, mu_s, m_u, m_d, m_s, Phi, Phibar, xi,
M_pi_fortran, Gamma_pi_fortran, M_K_fortran, Gamma_K_fortran

输出 CSV 追加计算结果与差值：
M_pi_julia, Gamma_pi_julia, dM_pi, dG_pi, M_K_julia, Gamma_K_julia, dM_K, dG_K

用法：
julia scripts/relaxtime/compare_meson_mass_fortran.jl input.csv [output.csv]
"""

using CSV

include("../../src/Constants_PNJL.jl")
include("../../src/relaxtime/EffectiveCouplings.jl")
include("../../src/relaxtime/MesonMass.jl")
include("../../src/relaxtime/MottTransition.jl")

using .MesonMass: solve_meson_mass, ensure_quark_params_has_A

function _require(row, key::Symbol)
	hasproperty(row, key) || error("Missing column: $(key)")
	return getproperty(row, key)
end

function _row_result(row)
	quark_params = (
		m = (u=_require(row, :m_u), d=_require(row, :m_d), s=_require(row, :m_s)),
		μ = (u=_require(row, :mu_u), d=_require(row, :mu_d), s=_require(row, :mu_s))
	)
	thermo_params = (
		T=_require(row, :T_fm),
		Φ=_require(row, :Phi),
		Φbar=_require(row, :Phibar),
		ξ=_require(row, :xi)
	)

	quark_params = ensure_quark_params_has_A(quark_params, thermo_params)

	pi_res = solve_meson_mass(:pi, quark_params, thermo_params;
							  initial_mass=_require(row, :M_pi_fortran),
							  initial_gamma=_require(row, :Gamma_pi_fortran))
	k_res = solve_meson_mass(:K, quark_params, thermo_params;
							 initial_mass=_require(row, :M_K_fortran),
							 initial_gamma=_require(row, :Gamma_K_fortran))

	return (
		M_pi_julia = pi_res.mass,
		Gamma_pi_julia = pi_res.gamma,
		dM_pi = pi_res.mass - _require(row, :M_pi_fortran),
		dG_pi = pi_res.gamma - _require(row, :Gamma_pi_fortran),
		M_K_julia = k_res.mass,
		Gamma_K_julia = k_res.gamma,
		dM_K = k_res.mass - _require(row, :M_K_fortran),
		dG_K = k_res.gamma - _require(row, :Gamma_K_fortran),
		converged_pi = pi_res.converged,
		converged_K = k_res.converged,
		residual_pi = pi_res.residual_norm,
		residual_K = k_res.residual_norm,
	)
end

function main()
	if length(ARGS) < 1
		error("Usage: julia compare_meson_mass_fortran.jl input.csv [output.csv]")
	end
	input_path = ARGS[1]
	output_path = length(ARGS) >= 2 ? ARGS[2] : nothing

	rows = CSV.File(input_path)
	results = NamedTuple[]
	for row in rows
		result = _row_result(row)
		push!(results, result)
	end

	if output_path === nothing
		for r in results
			println(r)
		end
	else
		CSV.write(output_path, results)
	end
end

main()
```

## 实验：A 热项尾部混合积分

脚本：`experiment_A_thermal_hybrid.jl`

用途：评估 A 的热项在高温参数区对动量上限截断的敏感性，并给出主区间/尾项的耗时与误差估计。

```julia
# 实验脚本：A 的热项采用“主区间 GL + 尾部变换积分(quadgk)”的混合策略。
#
# 目的：评估在不修改模块实现的前提下，
# - A_th_main = 4 * ∫_0^{p_c} dp  p^2/E * (f + fbar)
# - A_th_tail = 4 * ∫_{p_c}^{∞} dp p^2/E * (f + fbar)
# 并输出尾项误差估计/采样次数、主区间/尾项耗时。
#
# 默认输入：data/outputs/results/ak_compare.csv（取前 N 行）
# 默认输出：data/outputs/results/a_thermal_hybrid_experiment.csv
#
# 运行示例：
#   julia --project=. scripts/relaxtime/experiment_A_thermal_hybrid.jl
#
# 可选环境变量：
# - A_HYBRID_INPUT   输入 CSV 路径（默认 data/outputs/results/ak_compare.csv）
# - A_HYBRID_OUT     输出 CSV 路径（默认 data/outputs/results/a_thermal_hybrid_experiment.csv）
# - A_HYBRID_NROWS   读取前几行（默认 10；<=0 表示全量）
# - A_HYBRID_PC      主区间截断 p_c（默认 20.0）
# - A_HYBRID_NMAIN   主区间 GL 节点数（默认 128）
# - A_HYBRID_RTOL    尾项 quadgk 相对容差（默认 1e-10）
# - A_HYBRID_ATOL    尾项 quadgk 绝对容差（默认 0.0）
# - A_HYBRID_MAXEVAL 尾项 quadgk 最大评估次数（默认 20000）

using CSV
using QuadGK

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "OneLoopIntegrals.jl"))

using .Constants_PNJL: Λ_inv_fm
using .GaussLegendre: gauleg
using .OneLoopIntegrals: A

# -------------------------------
# PNJL 分布（用于尾部积分）：不对 exp_term 做下界 clamp
#
# 说明：现有 PNJL 分布实现对 exp_term 做了 min clamp(1e-200)。
# 在做 [p_c,∞) 积分时，严格意义上会让分布在极大 p 处不再继续衰减到 0，
# 造成“尾部积分到无穷”在数学上不再绝对收敛。
#
# 本实验脚本只在尾部/混合积分里使用“无下界 clamp”的稳定实现，允许自然下溢到 0。
@inline function _exp_pos_clamped(x::Float64)
	y = exp(x)
	if !isfinite(y)
		return 1e200
	end
	return y > 1e200 ? 1e200 : (y < 0.0 ? 0.0 : y)
end

@fastmath function pnjl_quark_distribution_no_min(E::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64)
	β = 1 / T
	exp_term = _exp_pos_clamped(-(E - μ) * β)
	exp_term2 = exp_term * exp_term
	exp_term3 = exp_term2 * exp_term

	numerator = Φ * exp_term + 2 * Φbar * exp_term2 + exp_term3
	denominator = 1 + 3 * Φ * exp_term + 3 * Φbar * exp_term2 + exp_term3
	return numerator / denominator
end

@fastmath function pnjl_antiquark_distribution_no_min(E::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64)
	β = 1 / T
	exp_term = _exp_pos_clamped(-(E + μ) * β)
	exp_term2 = exp_term * exp_term
	exp_term3 = exp_term2 * exp_term

	numerator = Φbar * exp_term + 2 * Φ * exp_term2 + exp_term3
	denominator = 1 + 3 * Φbar * exp_term + 3 * Φ * exp_term2 + exp_term3
	return numerator / denominator
end

@inline @fastmath function thermal_integrand_p(p::Float64, m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64)
	# p^2/E * (f + fbar)
	E = sqrt(p * p + m * m)
	f = pnjl_quark_distribution_no_min(E, μ, T, Φ, Φbar)
	fbar = pnjl_antiquark_distribution_no_min(E, μ, T, Φ, Φbar)
	val = p * p / E * (f + fbar)
	return isfinite(val) ? val : 0.0
end

function A_thermal_GL(m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64; pmax::Float64, n::Int)
	nodes, weights = gauleg(0.0, pmax, n)
	integral = 0.0
	@inbounds @simd for i in eachindex(nodes)
		integral += weights[i] * thermal_integrand_p(nodes[i], m, μ, T, Φ, Φbar)
	end
	return 4.0 * integral
end

function A_thermal_tail_quadgk(m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64;
	p_c::Float64, rtol::Float64, atol::Float64, maxevals::Int)

	eval_count = Ref{Int}(0)

	function integrand_t(t::Float64)
		eval_count[] += 1

		if t <= 0.0
			p = p_c
			return thermal_integrand_p(p, m, μ, T, Φ, Φbar)
		elseif t >= 1.0
			return 0.0
		end

		inv = 1.0 - t
		# p = p_c + t/(1-t)
		p = p_c + t / inv
		base = thermal_integrand_p(p, m, μ, T, Φ, Φbar)
		if base == 0.0
			return 0.0
		end

		jac = 1.0 / (inv * inv)
		val = base * jac
		return isfinite(val) ? val : 0.0
	end

	val, err = quadgk(integrand_t, 0.0, 1.0; rtol=rtol, atol=atol, maxevals=maxevals)
	return 4.0 * val, 4.0 * err, eval_count[]
end

function A_vacuum_analytic(m::Float64)
	# 与 OneLoopIntegrals.A 保持一致：A_vac = 4 * (- const_integral_term_A(m))
	return 4.0 * (-OneLoopIntegrals.const_integral_term_A(m))
end

@inline function rel_diff(a::Float64, b::Float64)
	denom = max(abs(b), 1e-12)
	return abs(a - b) / denom
end

# -------------------------------
# I/O
function main()
	input_path = get(ENV, "A_HYBRID_INPUT", joinpath(PROJECT_ROOT, "data", "outputs", "results", "ak_compare.csv"))
	output_path = get(ENV, "A_HYBRID_OUT", joinpath(PROJECT_ROOT, "data", "outputs", "results", "a_thermal_hybrid_experiment.csv"))

	nrows = parse(Int, get(ENV, "A_HYBRID_NROWS", "10"))
	p_c = parse(Float64, get(ENV, "A_HYBRID_PC", "20.0"))
	n_main = parse(Int, get(ENV, "A_HYBRID_NMAIN", "128"))
	rtol_tail = parse(Float64, get(ENV, "A_HYBRID_RTOL", "1e-10"))
	atol_tail = parse(Float64, get(ENV, "A_HYBRID_ATOL", "0.0"))
	maxevals_tail = parse(Int, get(ENV, "A_HYBRID_MAXEVAL", "20000"))

	rows_in = CSV.File(input_path)
	results = NamedTuple[]

	idx = 0
	for row in rows_in
		idx += 1
		if nrows > 0 && idx > nrows
			break
		end

		T = Float64(row.T)
		Φ = Float64(row.Phi)
		Φbar = Float64(row.Phibar)

		m_u = Float64(row.m_u)
		m_s = Float64(row.m_s)
		μ_u = Float64(row.mu_u)
		μ_s = Float64(row.mu_s)

		A_u_fortran = hasproperty(row, :A_u_fortran) ? Float64(row.A_u_fortran) : NaN
		A_s_fortran = hasproperty(row, :A_s_fortran) ? Float64(row.A_s_fortran) : NaN

		# baseline：项目当前推荐的默认配置 (p<=10, n=32)
		nodes10, weights10 = gauleg(0.0, 10.0, 32)

		A_u_10_32 = A(m_u, μ_u, T, Φ, Φbar, nodes10, weights10)
		A_s_10_32 = A(m_s, μ_s, T, Φ, Φbar, nodes10, weights10)

		# hybrid：热项主区间 + 尾项
		t0 = time_ns()
		A_u_th_main = A_thermal_GL(m_u, μ_u, T, Φ, Φbar; pmax=p_c, n=n_main)
		A_s_th_main = A_thermal_GL(m_s, μ_s, T, Φ, Φbar; pmax=p_c, n=n_main)
		t1 = time_ns()

		t2 = time_ns()
		A_u_th_tail, A_u_tail_err, A_u_tail_neval = A_thermal_tail_quadgk(m_u, μ_u, T, Φ, Φbar;
			p_c=p_c, rtol=rtol_tail, atol=atol_tail, maxevals=maxevals_tail)
		A_s_th_tail, A_s_tail_err, A_s_tail_neval = A_thermal_tail_quadgk(m_s, μ_s, T, Φ, Φbar;
			p_c=p_c, rtol=rtol_tail, atol=atol_tail, maxevals=maxevals_tail)
		t3 = time_ns()

		A_u_vac = A_vacuum_analytic(m_u)
		A_s_vac = A_vacuum_analytic(m_s)

		A_u_hybrid = A_u_vac + A_u_th_main + A_u_th_tail
		A_s_hybrid = A_s_vac + A_s_th_main + A_s_th_tail

		push!(results, (
			idx=idx,
			T=T,
			mu_u=μ_u,
			mu_s=μ_s,
			m_u=m_u,
			m_s=m_s,
			Phi=Φ,
			Phibar=Φbar,
			p_c=p_c,
			n_main=n_main,
			rtol_tail=rtol_tail,
			atol_tail=atol_tail,
			maxevals_tail=maxevals_tail,
			A_u_fortran=A_u_fortran,
			A_s_fortran=A_s_fortran,
			A_u_10_32=A_u_10_32,
			A_s_10_32=A_s_10_32,
			A_u_vac=A_u_vac,
			A_s_vac=A_s_vac,
			A_u_th_main=A_u_th_main,
			A_s_th_main=A_s_th_main,
			A_u_th_tail=A_u_th_tail,
			A_s_th_tail=A_s_th_tail,
			A_u_tail_err=A_u_tail_err,
			A_s_tail_err=A_s_tail_err,
			A_u_tail_neval=A_u_tail_neval,
			A_s_tail_neval=A_s_tail_neval,
			A_u_hybrid=A_u_hybrid,
			A_s_hybrid=A_s_hybrid,
			dA_u_fortran_vs_hybrid=(isfinite(A_u_fortran) ? (A_u_hybrid - A_u_fortran) : NaN),
			dA_s_fortran_vs_hybrid=(isfinite(A_s_fortran) ? (A_s_hybrid - A_s_fortran) : NaN),
			rel_u_10_vs_hybrid=rel_diff(A_u_10_32, A_u_hybrid),
			rel_s_10_vs_hybrid=rel_diff(A_s_10_32, A_s_hybrid),
			t_main_s=1e-9 * (t1 - t0),
			t_tail_s=1e-9 * (t3 - t2),
		))
	end

	CSV.write(output_path, results)
	println("Wrote: ", output_path)
	println("Rows: ", length(results))
end

main()
```

## 实验：B0 注入残差对比

脚本：`experiment_b0_injection.jl`

用途：在给定参数点用 Fortran B0 替换 Julia B0，比较极点方程残差变化，用于归因差异是否主要来自 B0。

```julia
"""
实验脚本：用 Fortran B0 替换 Julia B0（仅在给定参数点）并比较介子极点方程残差。

输入：
- data/outputs/results/meson_fortran_input.csv
- data/outputs/results/b0_fortran_output.csv

输出：
- data/outputs/results/meson_b0_injection_experiment.csv

说明：
- 若 Fortran B0 缺少某些参数点（例如负 λ），会回退到 Julia B0，并在输出中标记。
"""

using CSV

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "OneLoopIntegrals.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))

using .Constants_PNJL: N_color, G_fm2, K_fm5
using .GaussLegendre: gauleg
using .OneLoopIntegrals: A, B0
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings

const FACTOR = -N_color / (8π^2)

function _key(vals...; digits::Int=10)
	return Tuple(round(Float64(v); digits=digits) for v in vals)
end

function _load_b0_map(path::String; digits::Int=10)
	rows = CSV.File(path)
	map = Dict{Tuple,Tuple{Float64,Float64}}()
	for row in rows
		lam = hasproperty(row, :lam) ? row.lam : row.λ
		key = _key(lam, row.k, row.m1, row.m2, row.mu1, row.mu2, row.T, row.Phi, row.Phibar; digits=digits)
		map[key] = (row.B0_re_fortran, row.B0_im_fortran)
	end
	return map
end

function _lookup_b0(map, lam, k, m1, m2, mu1, mu2, T, Phi, Phibar; digits::Int=10)
	key = _key(lam, k, m1, m2, mu1, mu2, T, Phi, Phibar; digits=digits)
	return get(map, key, nothing)
end

function _polarization_from_b0(channel::Symbol, k0::Float64, gamma::Float64, k_norm::Float64,
							   m1::Float64, m2::Float64, mu1::Float64, mu2::Float64,
							   T::Float64, Phi::Float64, Phibar::Float64,
							   A1::Float64, A2::Float64, num_s_quark::Int,
							   b0_re::Float64, b0_im::Float64)
	λ = k0 + mu1 - mu2

	prefactor = k_norm^2 - λ^2
	if channel == :P
		prefactor += (m1 - m2)^2
	elseif channel == :S
		prefactor += (m1 + m2)^2
	else
		error("Unsupported channel: $channel")
	end

	prefactor += (gamma^2) / 4.0

	real_part = A1 + A2 + prefactor * b0_re - gamma * λ * b0_im
	imag_part = prefactor * b0_im + gamma * λ * b0_re

	return FACTOR * real_part, FACTOR * imag_part
end

function _meson_residuals(row, b0_map)
	T = row.T_fm
	Phi = row.Phi
	Phibar = row.Phibar

	m_u = row.m_u
	m_d = row.m_d
	m_s = row.m_s
	mu_u = row.mu_u
	mu_d = row.mu_d
	mu_s = row.mu_s

	# A 与 K
	nodes_p, weights_p = gauleg(0.0, 10.0, 32)
	A_u = A(m_u, mu_u, T, Phi, Phibar, nodes_p, weights_p)
	A_s = A(m_s, mu_s, T, Phi, Phibar, nodes_p, weights_p)
	G_u = calculate_G_from_A(A_u, m_u)
	G_s = calculate_G_from_A(A_s, m_s)
	K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

	# --- π ---
	M_pi = row.M_pi_fortran
	G_pi = row.Gamma_pi_fortran
	b0_pi = _lookup_b0(b0_map, M_pi, 0.0, m_u, m_u, mu_u, mu_d, T, Phi, Phibar)
	if b0_pi === nothing
		b0_pi_re, b0_pi_im = B0(M_pi, 0.0, m_u, mu_u, m_u, mu_d, T; Φ=Phi, Φbar=Phibar)
		pi_b0_source = "julia"
	else
		b0_pi_re, b0_pi_im = b0_pi
		pi_b0_source = "fortran"
	end
	Π_pi_re, Π_pi_im = _polarization_from_b0(:P, M_pi, G_pi, 0.0, m_u, m_u, mu_u, mu_d, T, Phi, Phibar,
											 A_u, A_u, 0, b0_pi_re, b0_pi_im)
	f_pi = 1.0 - 4.0 * K_coeffs.K123_plus * ComplexF64(Π_pi_re, Π_pi_im)

	# --- K ---
	M_K = row.M_K_fortran
	G_K = row.Gamma_K_fortran
	b0_K = _lookup_b0(b0_map, M_K, 0.0, m_u, m_s, mu_u, mu_s, T, Phi, Phibar)
	# 对称化需要 λ 与 -λ；若缺失则回退 Julia
	b0_K_neg = _lookup_b0(b0_map, -M_K, 0.0, m_u, m_s, mu_u, mu_s, T, Phi, Phibar)

	if b0_K === nothing
		b0_K_re, b0_K_im = B0(M_K, 0.0, m_u, mu_u, m_s, mu_s, T; Φ=Phi, Φbar=Phibar)
		K_b0_source = "julia"
	else
		b0_K_re, b0_K_im = b0_K
		K_b0_source = "fortran"
	end

	if b0_K_neg === nothing
		b0_Kn_re, b0_Kn_im = B0(-M_K, 0.0, m_u, mu_u, m_s, mu_s, T; Φ=Phi, Φbar=Phibar)
		K_b0_source = K_b0_source == "fortran" ? "mixed" : "julia"
	else
		b0_Kn_re, b0_Kn_im = b0_K_neg
		if K_b0_source == "fortran"
			K_b0_source = "fortran"
		else
			K_b0_source = "mixed"
		end
	end

	# num_s_quark=1 的对称化
	b0_K_re = 0.5 * (b0_K_re + b0_Kn_re)
	b0_K_im = 0.5 * (b0_K_im + b0_Kn_im)

	Π_K_re, Π_K_im = _polarization_from_b0(:P, M_K, G_K, 0.0, m_u, m_s, mu_u, mu_s, T, Phi, Phibar,
										   A_u, A_s, 1, b0_K_re, b0_K_im)
	f_K = 1.0 - 4.0 * K_coeffs.K4567_plus * ComplexF64(Π_K_re, Π_K_im)

	# Julia B0 计算残差对照
	b0_pi_j_re, b0_pi_j_im = B0(M_pi, 0.0, m_u, mu_u, m_u, mu_d, T; Φ=Phi, Φbar=Phibar)
	Π_pi_j_re, Π_pi_j_im = _polarization_from_b0(:P, M_pi, G_pi, 0.0, m_u, m_u, mu_u, mu_d, T, Phi, Phibar,
												 A_u, A_u, 0, b0_pi_j_re, b0_pi_j_im)
	f_pi_j = 1.0 - 4.0 * K_coeffs.K123_plus * ComplexF64(Π_pi_j_re, Π_pi_j_im)

	b0_K_j_re, b0_K_j_im = B0(M_K, 0.0, m_u, mu_u, m_s, mu_s, T; Φ=Phi, Φbar=Phibar)
	b0_Kn_j_re, b0_Kn_j_im = B0(-M_K, 0.0, m_u, mu_u, m_s, mu_s, T; Φ=Phi, Φbar=Phibar)
	b0_K_j_re = 0.5 * (b0_K_j_re + b0_Kn_j_re)
	b0_K_j_im = 0.5 * (b0_K_j_im + b0_Kn_j_im)
	Π_K_j_re, Π_K_j_im = _polarization_from_b0(:P, M_K, G_K, 0.0, m_u, m_s, mu_u, mu_s, T, Phi, Phibar,
											   A_u, A_s, 1, b0_K_j_re, b0_K_j_im)
	f_K_j = 1.0 - 4.0 * K_coeffs.K4567_plus * ComplexF64(Π_K_j_re, Π_K_j_im)

	return (
		f_pi_fortran_re = real(f_pi),
		f_pi_fortran_im = imag(f_pi),
		f_pi_julia_re = real(f_pi_j),
		f_pi_julia_im = imag(f_pi_j),
		f_K_fortran_re = real(f_K),
		f_K_fortran_im = imag(f_K),
		f_K_julia_re = real(f_K_j),
		f_K_julia_im = imag(f_K_j),
		pi_b0_source = pi_b0_source,
		K_b0_source = K_b0_source,
	)
end

function main()
	b0_path = joinpath(PROJECT_ROOT, "data", "outputs", "results", "b0_fortran_output.csv")
	input_path = joinpath(PROJECT_ROOT, "data", "outputs", "results", "meson_fortran_input.csv")
	output_path = joinpath(PROJECT_ROOT, "data", "outputs", "results", "meson_b0_injection_experiment.csv")

	b0_map = _load_b0_map(b0_path)
	rows = CSV.File(input_path)
	results = NamedTuple[]

	for row in rows
		res = _meson_residuals(row, b0_map)
		push!(results, merge((T_fm=row.T_fm,), res))
	end

	CSV.write(output_path, results)
end

main()
```

## 实验：全替换（Fortran A/K + Fortran B0）

脚本：`experiment_full_fortran_params.jl`

用途：使用 Fortran A/K 与 Fortran B0 全替换，计算 π/K 的极点方程残差，并给出交叉组合残差用于定位差异来源。

```julia
# 实验脚本：使用 Fortran A/K 与 Fortran B0 全替换，计算介子极点方程残差。
#
# 输入：
# - data/outputs/results/meson_fortran_input.csv
# - data/outputs/results/b0_fortran_output.csv
# - data/outputs/results/ak_fortran_output.csv
#
# 输出：
# - data/outputs/results/meson_full_fortran_params_experiment.csv

using CSV

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "OneLoopIntegrals.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))

using .Constants_PNJL: N_color, G_fm2, K_fm5
using .Constants_PNJL: Λ_inv_fm
using .GaussLegendre: gauleg
using .OneLoopIntegrals: A, B0
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings

const FACTOR = -N_color / (8π^2)

function _key(vals...; digits::Int=10)
	return Tuple(round(Float64(v); digits=digits) for v in vals)
end

function _load_b0_map(path::String; digits::Int=10)
	rows = CSV.File(path)
	map = Dict{Tuple,Tuple{Float64,Float64}}()
	for row in rows
		lam = hasproperty(row, :lam) ? row.lam : row.λ
		key = _key(lam, row.k, row.m1, row.m2, row.mu1, row.mu2, row.T, row.Phi, row.Phibar; digits=digits)
		map[key] = (row.B0_re_fortran, row.B0_im_fortran)
	end
	return map
end

function _load_ak_map(path::String; digits::Int=10)
	rows = CSV.File(path)
	map = Dict{Tuple,NamedTuple}()
	for row in rows
		key = _key(row.T, row.mu_u, row.mu_s, row.m_u, row.m_s, row.Phi, row.Phibar; digits=digits)
		map[key] = (
			A_u=row.A_u,
			A_s=row.A_s,
			K123_plus=row.K123_plus,
			K4567_plus=row.K4567_plus,
		)
	end
	return map
end

function _lookup_b0(map, lam, k, m1, m2, mu1, mu2, T, Phi, Phibar; digits::Int=10, atol::Float64=1e-7)
	key = _key(lam, k, m1, m2, mu1, mu2, T, Phi, Phibar; digits=digits)
	hit = get(map, key, nothing)
	hit === nothing || return hit

	target = (lam, k, m1, m2, mu1, mu2, T, Phi, Phibar)
	for (kkey, v) in map
		if all(isapprox(kkey[i], target[i]; atol=atol) for i in 1:length(target))
			return v
		end
	end
	return nothing
end

function _lookup_ak(map, T, mu_u, mu_s, m_u, m_s, Phi, Phibar; digits::Int=10, atol::Float64=1e-7)
	key = _key(T, mu_u, mu_s, m_u, m_s, Phi, Phibar; digits=digits)
	hit = get(map, key, nothing)
	hit === nothing || return hit

	target = (T, mu_u, mu_s, m_u, m_s, Phi, Phibar)
	for (kkey, v) in map
		if all(isapprox(kkey[i], target[i]; atol=atol) for i in 1:length(target))
			return v
		end
	end
	return nothing
end

function _lookup_b0_with_meta(map, lam, k, m1, m2, mu1, mu2, T, Phi, Phibar; digits::Int=10, atol::Float64=1e-7)
	key = _key(lam, k, m1, m2, mu1, mu2, T, Phi, Phibar; digits=digits)
	hit = get(map, key, nothing)
	hit === nothing || return hit, :exact

	target = (lam, k, m1, m2, mu1, mu2, T, Phi, Phibar)
	for (kkey, v) in map
		if all(isapprox(kkey[i], target[i]; atol=atol) for i in 1:length(target))
			return v, :approx
		end
	end
	return nothing, :missing
end

function _lookup_ak_with_meta(map, T, mu_u, mu_s, m_u, m_s, Phi, Phibar; digits::Int=10, atol::Float64=1e-7)
	key = _key(T, mu_u, mu_s, m_u, m_s, Phi, Phibar; digits=digits)
	hit = get(map, key, nothing)
	hit === nothing || return hit, :exact

	target = (T, mu_u, mu_s, m_u, m_s, Phi, Phibar)
	for (kkey, v) in map
		if all(isapprox(kkey[i], target[i]; atol=atol) for i in 1:length(target))
			return v, :approx
		end
	end
	return nothing, :missing
end

function _prefactor(channel::Symbol, k0::Float64, gamma::Float64, k_norm::Float64,
					m1::Float64, m2::Float64, mu1::Float64, mu2::Float64)
	λ = k0 + mu1 - mu2
	prefactor = k_norm^2 - λ^2
	if channel == :P
		prefactor += (m1 - m2)^2
	elseif channel == :S
		prefactor += (m1 + m2)^2
	else
		error("Unsupported channel: $channel")
	end

	prefactor += (gamma^2) / 4.0
	return λ, prefactor
end

function _polarization_from_b0(channel::Symbol, k0::Float64, gamma::Float64, k_norm::Float64,
							   m1::Float64, m2::Float64, mu1::Float64, mu2::Float64,
							   T::Float64, Phi::Float64, Phibar::Float64,
							   A1::Float64, A2::Float64, num_s_quark::Int,
							   b0_re::Float64, b0_im::Float64)
	λ, prefactor = _prefactor(channel, k0, gamma, k_norm, m1, m2, mu1, mu2)

	real_part = A1 + A2 + prefactor * b0_re - gamma * λ * b0_im
	imag_part = prefactor * b0_im + gamma * λ * b0_re

	return FACTOR * real_part, FACTOR * imag_part
end

function _meson_residuals(row, b0_map, ak_map)
	T = row.T_fm
	Phi = row.Phi
	Phibar = row.Phibar

	m_u = row.m_u
	m_d = row.m_d
	m_s = row.m_s
	mu_u = row.mu_u
	mu_d = row.mu_d
	mu_s = row.mu_s

	ak, ak_method = _lookup_ak_with_meta(ak_map, T, mu_u, mu_s, m_u, m_s, Phi, Phibar)
	if ak === nothing
		error("Missing Fortran A/K entry for T=$(T)")
	end

	# --- π (Fortran A/K + Fortran B0) ---
	M_pi = row.M_pi_fortran
	G_pi = row.Gamma_pi_fortran
	b0_pi, b0_pi_method = _lookup_b0_with_meta(b0_map, M_pi, 0.0, m_u, m_u, mu_u, mu_d, T, Phi, Phibar)
	if b0_pi === nothing
		error("Missing Fortran B0 entry for pi")
	end
	b0_pi_re, b0_pi_im = b0_pi

	Π_pi_re, Π_pi_im = _polarization_from_b0(:P, M_pi, G_pi, 0.0, m_u, m_u, mu_u, mu_d, T, Phi, Phibar,
											 ak.A_u, ak.A_u, 0, b0_pi_re, b0_pi_im)
	f_pi_fortran = 1.0 - 4.0 * ak.K123_plus * ComplexF64(Π_pi_re, Π_pi_im)

	# --- K (Fortran A/K + Fortran B0) ---
	M_K = row.M_K_fortran
	G_K = row.Gamma_K_fortran
	b0_K, b0_K_method = _lookup_b0_with_meta(b0_map, M_K, 0.0, m_u, m_s, mu_u, mu_s, T, Phi, Phibar)
	b0_Kn, b0_Kn_method = _lookup_b0_with_meta(b0_map, -M_K, 0.0, m_u, m_s, mu_u, mu_s, T, Phi, Phibar)
	if b0_K === nothing || b0_Kn === nothing
		error("Missing Fortran B0 entry for K")
	end
	b0_K_re, b0_K_im = b0_K
	b0_Kn_re, b0_Kn_im = b0_Kn
	b0_K_re = 0.5 * (b0_K_re + b0_Kn_re)
	b0_K_im = 0.5 * (b0_K_im + b0_Kn_im)

	Π_K_re, Π_K_im = _polarization_from_b0(:P, M_K, G_K, 0.0, m_u, m_s, mu_u, mu_s, T, Phi, Phibar,
										   ak.A_u, ak.A_s, 1, b0_K_re, b0_K_im)
	f_K_fortran = 1.0 - 4.0 * ak.K4567_plus * ComplexF64(Π_K_re, Π_K_im)

	# Julia 基线（Julia A/K + Julia B0）
	# 说明：高温/轻质量参数区对热积分截断更敏感，这里使用更稳健的默认配置。
	nodes_p, weights_p = gauleg(0.0, 20.0, 16)
	A_u_j = A(m_u, mu_u, T, Phi, Phibar, nodes_p, weights_p)
	A_s_j = A(m_s, mu_s, T, Phi, Phibar, nodes_p, weights_p)

	# A 的额外诊断：节点数/上限敏感性（不改模块实现，仅改实验脚本采样）
	nodes_p64, weights_p64 = gauleg(0.0, 10.0, 64)
	A_u_j_64 = A(m_u, mu_u, T, Phi, Phibar, nodes_p64, weights_p64)
	A_s_j_64 = A(m_s, mu_s, T, Phi, Phibar, nodes_p64, weights_p64)

	nodes_p128, weights_p128 = gauleg(0.0, 10.0, 128)
	A_u_j_128 = A(m_u, mu_u, T, Phi, Phibar, nodes_p128, weights_p128)
	A_s_j_128 = A(m_s, mu_s, T, Phi, Phibar, nodes_p128, weights_p128)

	# 参考 Fortran 侧常用的热积分上限 15 fm⁻¹（用于验证 10 是否已收敛）
	nodes_p15, weights_p15 = gauleg(0.0, 15.0, 128)
	A_u_j_15 = A(m_u, mu_u, T, Phi, Phibar, nodes_p15, weights_p15)
	A_s_j_15 = A(m_s, mu_s, T, Phi, Phibar, nodes_p15, weights_p15)

	nodes_pΛ, weights_pΛ = gauleg(0.0, Λ_inv_fm, 128)
	A_u_j_Λ = A(m_u, mu_u, T, Phi, Phibar, nodes_pΛ, weights_pΛ)
	A_s_j_Λ = A(m_s, mu_s, T, Phi, Phibar, nodes_pΛ, weights_pΛ)
	G_u_j = calculate_G_from_A(A_u_j, m_u)
	G_s_j = calculate_G_from_A(A_s_j, m_s)
	K_j = calculate_effective_couplings(G_fm2, K_fm5, G_u_j, G_s_j)

	b0_pi_j_re, b0_pi_j_im = B0(M_pi, 0.0, m_u, mu_u, m_u, mu_d, T; Φ=Phi, Φbar=Phibar)
	Π_pi_j_re, Π_pi_j_im = _polarization_from_b0(:P, M_pi, G_pi, 0.0, m_u, m_u, mu_u, mu_d, T, Phi, Phibar,
												 A_u_j, A_u_j, 0, b0_pi_j_re, b0_pi_j_im)
	f_pi_julia = 1.0 - 4.0 * K_j.K123_plus * ComplexF64(Π_pi_j_re, Π_pi_j_im)

	b0_K_j_re, b0_K_j_im = B0(M_K, 0.0, m_u, mu_u, m_s, mu_s, T; Φ=Phi, Φbar=Phibar)
	b0_Kn_j_re, b0_Kn_j_im = B0(-M_K, 0.0, m_u, mu_u, m_s, mu_s, T; Φ=Phi, Φbar=Phibar)
	b0_K_j_re = 0.5 * (b0_K_j_re + b0_Kn_j_re)
	b0_K_j_im = 0.5 * (b0_K_j_im + b0_Kn_j_im)
	Π_K_j_re, Π_K_j_im = _polarization_from_b0(:P, M_K, G_K, 0.0, m_u, m_s, mu_u, mu_s, T, Phi, Phibar,
											   A_u_j, A_s_j, 1, b0_K_j_re, b0_K_j_im)
	f_K_julia = 1.0 - 4.0 * K_j.K4567_plus * ComplexF64(Π_K_j_re, Π_K_j_im)

	# 交叉组合：Fortran B0 + Julia A
	Π_pi_fb_ja_re, Π_pi_fb_ja_im = _polarization_from_b0(:P, M_pi, G_pi, 0.0, m_u, m_u, mu_u, mu_d, T, Phi, Phibar,
														 A_u_j, A_u_j, 0, b0_pi_re, b0_pi_im)
	f_pi_fb_ja = 1.0 - 4.0 * ak.K123_plus * ComplexF64(Π_pi_fb_ja_re, Π_pi_fb_ja_im)

	Π_K_fb_ja_re, Π_K_fb_ja_im = _polarization_from_b0(:P, M_K, G_K, 0.0, m_u, m_s, mu_u, mu_s, T, Phi, Phibar,
													   A_u_j, A_s_j, 1, b0_K_re, b0_K_im)
	f_K_fb_ja = 1.0 - 4.0 * ak.K4567_plus * ComplexF64(Π_K_fb_ja_re, Π_K_fb_ja_im)

	# 交叉组合：Julia B0 + Fortran A
	Π_pi_jb_fa_re, Π_pi_jb_fa_im = _polarization_from_b0(:P, M_pi, G_pi, 0.0, m_u, m_u, mu_u, mu_d, T, Phi, Phibar,
														 ak.A_u, ak.A_u, 0, b0_pi_j_re, b0_pi_j_im)
	f_pi_jb_fa = 1.0 - 4.0 * ak.K123_plus * ComplexF64(Π_pi_jb_fa_re, Π_pi_jb_fa_im)

	Π_K_jb_fa_re, Π_K_jb_fa_im = _polarization_from_b0(:P, M_K, G_K, 0.0, m_u, m_s, mu_u, mu_s, T, Phi, Phibar,
													   ak.A_u, ak.A_s, 1, b0_K_j_re, b0_K_j_im)
	f_K_jb_fa = 1.0 - 4.0 * ak.K4567_plus * ComplexF64(Π_K_jb_fa_re, Π_K_jb_fa_im)

	λ_pi, pref_pi = _prefactor(:P, M_pi, G_pi, 0.0, m_u, m_u, mu_u, mu_d)
	λ_K, pref_K = _prefactor(:P, M_K, G_K, 0.0, m_u, m_s, mu_u, mu_s)

	return (
		ak_method = String(ak_method),
		b0_pi_method = String(b0_pi_method),
		b0_K_method = String(b0_K_method),
		b0_Kn_method = String(b0_Kn_method),
		A_u_fortran = ak.A_u,
		A_s_fortran = ak.A_s,
		A_u_julia = A_u_j,
		A_s_julia = A_s_j,
		A_u_julia_64 = A_u_j_64,
		A_s_julia_64 = A_s_j_64,
		A_u_julia_128 = A_u_j_128,
		A_s_julia_128 = A_s_j_128,
		A_u_julia_15 = A_u_j_15,
		A_s_julia_15 = A_s_j_15,
		A_u_julia_Lambda = A_u_j_Λ,
		A_s_julia_Lambda = A_s_j_Λ,
		dA_u = A_u_j - ak.A_u,
		dA_s = A_s_j - ak.A_s,
		dA_u_64 = A_u_j_64 - ak.A_u,
		dA_s_64 = A_s_j_64 - ak.A_s,
		dA_u_128 = A_u_j_128 - ak.A_u,
		dA_s_128 = A_s_j_128 - ak.A_s,
		dA_u_15 = A_u_j_15 - ak.A_u,
		dA_s_15 = A_s_j_15 - ak.A_s,
		dA_u_Lambda = A_u_j_Λ - ak.A_u,
		dA_s_Lambda = A_s_j_Λ - ak.A_s,
		K123_fortran = ak.K123_plus,
		K4567_fortran = ak.K4567_plus,
		K123_julia = K_j.K123_plus,
		K4567_julia = K_j.K4567_plus,
		b0_pi_fortran_re = b0_pi_re,
		b0_pi_fortran_im = b0_pi_im,
		b0_pi_julia_re = b0_pi_j_re,
		b0_pi_julia_im = b0_pi_j_im,
		b0_K_fortran_re = b0_K_re,
		b0_K_fortran_im = b0_K_im,
		b0_K_julia_re = b0_K_j_re,
		b0_K_julia_im = b0_K_j_im,
		lambda_pi = λ_pi,
		prefactor_pi = pref_pi,
		lambda_K = λ_K,
		prefactor_K = pref_K,
		f_pi_fortran_re = real(f_pi_fortran),
		f_pi_fortran_im = imag(f_pi_fortran),
		f_pi_fb_ja_re = real(f_pi_fb_ja),
		f_pi_fb_ja_im = imag(f_pi_fb_ja),
		f_pi_jb_fa_re = real(f_pi_jb_fa),
		f_pi_jb_fa_im = imag(f_pi_jb_fa),
		f_pi_julia_re = real(f_pi_julia),
		f_pi_julia_im = imag(f_pi_julia),
		f_K_fortran_re = real(f_K_fortran),
		f_K_fortran_im = imag(f_K_fortran),
		f_K_fb_ja_re = real(f_K_fb_ja),
		f_K_fb_ja_im = imag(f_K_fb_ja),
		f_K_jb_fa_re = real(f_K_jb_fa),
		f_K_jb_fa_im = imag(f_K_jb_fa),
		f_K_julia_re = real(f_K_julia),
		f_K_julia_im = imag(f_K_julia),
	)
end

function main()
	b0_path = joinpath(PROJECT_ROOT, "data", "outputs", "results", "b0_fortran_output.csv")
	ak_path = joinpath(PROJECT_ROOT, "data", "outputs", "results", "ak_fortran_output.csv")
	input_path = joinpath(PROJECT_ROOT, "data", "outputs", "results", "meson_fortran_input.csv")
	output_path = joinpath(PROJECT_ROOT, "data", "outputs", "results", "meson_full_fortran_params_experiment.csv")

	b0_map = _load_b0_map(b0_path)
	ak_map = _load_ak_map(ak_path)

	rows = CSV.File(input_path)
	results = NamedTuple[]

	for row in rows
		res = _meson_residuals(row, b0_map, ak_map)
		push!(results, merge((T_fm=row.T_fm,), res))
	end

	CSV.write(output_path, results)
end

main()
```
