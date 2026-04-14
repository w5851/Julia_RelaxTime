using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const BASELINE_PATH = joinpath(PROJECT_ROOT, "tests", "baselines", "relaxtime", "baseline_t190_mixed_p_chain_v1.csv")

include(joinpath(PROJECT_ROOT, "scripts", "analysis", "relaxtime", "t190_sigma_chain_decomposition_lib.jl"))

function _trapz(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    n == length(y) || error("trapz: length mismatch")
    n <= 1 && return 0.0
    acc = 0.0
    @inbounds for i in 1:(n - 1)
        acc += 0.5 * (y[i] + y[i + 1]) * (x[i + 1] - x[i])
    end
    return acc
end

function _load_baseline(path::String)
    isfile(path) || error("baseline CSV not found: $path")
    lines = readlines(path)
    length(lines) >= 2 || error("baseline CSV is empty: $path")

    metrics = Dict{String, NamedTuple{(:expected, :rtol, :atol), Tuple{Float64, Float64, Float64}}}()
    for line in lines[2:end]
        row = strip(line)
        isempty(row) && continue
        startswith(row, "#") && continue
        cols = split(row, ',')
        length(cols) == 4 || error("invalid baseline row: $line")
        metric = strip(cols[1])
        expected = parse(Float64, strip(cols[2]))
        rtol = parse(Float64, strip(cols[3]))
        atol = parse(Float64, strip(cols[4]))
        metrics[metric] = (expected=expected, rtol=rtol, atol=atol)
    end
    return metrics
end

function _compute_chain_metrics()
    stA = build_state_point(190.0, 0.0, -0.10)
    stB = build_state_point(190.0, 0.0, -0.08)
    process = :uubar_to_ddbar

    thA = process_threshold_info(process, stA.quark_params)
    thB = process_threshold_info(process, stB.quark_params)

    ds_vals = collect(range(1e-6, 2.0; length=36))
    detM_A = Float64[]
    detM_B = Float64[]
    Dmix_A = Float64[]
    Dmix_B = Float64[]

    detM_point_ratio = NaN
    Dmix_point_ratio = NaN

    for (idx, ds) in enumerate(ds_vals)
        sA = thA.s_th + ds
        sB = thB.s_th + ds

        tbA = Main.TotalCrossSection.calculate_t_bounds(sA, thA.mi, thA.mj, thA.mc, thA.md)
        tbB = Main.TotalCrossSection.calculate_t_bounds(sB, thB.mi, thB.mj, thB.mc, thB.md)
        tA = 0.5 * (tbA.t_min + tbA.t_max)
        tB = 0.5 * (tbB.t_min + tbB.t_max)

        pA = decompose_mixed_p_propagator_chain(process, sA, tA, stA.quark_params, stA.thermo_params, stA.K_coeffs)
        pB = decompose_mixed_p_propagator_chain(process, sB, tB, stB.quark_params, stB.thermo_params, stB.K_coeffs)

        push!(detM_A, pA.abs_detM_sq)
        push!(detM_B, pB.abs_detM_sq)
        push!(Dmix_A, pA.abs_D_mixed_P_sq)
        push!(Dmix_B, pB.abs_D_mixed_P_sq)

        if idx == 1
            detM_point_ratio = pB.abs_detM_sq / pA.abs_detM_sq
            Dmix_point_ratio = pB.abs_D_mixed_P_sq / pA.abs_D_mixed_P_sq
        end
    end

    detM_area_ratio = _trapz(ds_vals, detM_B) / _trapz(ds_vals, detM_A)
    Dmix_area_ratio = _trapz(ds_vals, Dmix_B) / _trapz(ds_vals, Dmix_A)

    return Dict(
        "ratio_detM_area_B_over_A" => detM_area_ratio,
        "ratio_Dmixed_area_B_over_A" => Dmix_area_ratio,
        "ratio_detM_point_B_over_A" => detM_point_ratio,
        "ratio_Dmixed_point_B_over_A" => Dmix_point_ratio,
    )
end

@testset "T190 mixed_P chain regression" begin
    baseline = _load_baseline(BASELINE_PATH)
    actual = _compute_chain_metrics()

    for metric in sort(collect(keys(baseline)))
        @test haskey(actual, metric)
        expected = baseline[metric]
        @test isapprox(actual[metric], expected.expected; rtol=expected.rtol, atol=expected.atol)
    end
end
