# CEPDetector.jl 单元测试
#
# 测试内容：
# 1. find_cep 接口存在
# 2. CEPResult 结构
# 3. 内部辅助函数可访问

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const PHASE_REFERENCE_TRHO = joinpath(PROJECT_ROOT, "data", "reference", "pnjl", "trho_scan_xi0.0.csv")

function _load_reference_curves()
    return Models.load_curves_from_trho_csv(PHASE_REFERENCE_TRHO; xi=0.0, min_points=3)
end

function _interpolate_curve(curves::Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}, T_target::Float64)
    temps = sort(collect(keys(curves)))
    T_below = nothing
    T_above = nothing

    for T in temps
        if T < T_target
            T_below = T
        elseif T > T_target && T_above === nothing
            T_above = T
        end
    end

    (T_below === nothing || T_above === nothing) && error("cannot bracket T=$T_target")
    mu_below, rho_below = curves[T_below]
    mu_above, rho_above = curves[T_above]
    length(rho_below) == length(rho_above) || error("nonuniform rho grid around T=$T_target")

    alpha = (T_target - T_below) / (T_above - T_below)
    return mu_below .+ alpha .* (mu_above .- mu_below), rho_below
end

# ============================================================================

@testset "CEPDetector" begin

    @testset "find_cep 接口存在" begin
        @test isdefined(Models, :find_cep)
        @test Models.find_cep isa Function
    end

    @testset "CEPResult 结构" begin
        @test isdefined(Models, :CEPResult)
        res = Models.CEPResult(
            found=false,
            T_cep_MeV=NaN,
            mu_cep_MeV=NaN,
            uncertainty_T_MeV=NaN,
            eval_count=0,
            unknown_count=0,
            reason="test",
            method=:bisection,
        )
        @test res.found == false
        @test res.reason == "test"
        legacy_not_found = Models.CEPResult(found=false, T_cep_MeV=130.0, mu_cep_MeV=295.0)
        @test isnan(legacy_not_found.T_cep_MeV)
        @test isnan(legacy_not_found.mu_cep_MeV)
        ambiguous = Models.CEPResult(
            result_status=:ambiguous,
            T_last_first_order_MeV=130.0,
            T_first_monotone_MeV=131.0,
        )
        @test ambiguous.result_status == :ambiguous
        @test !ambiguous.found
        @test isnan(ambiguous.T_cep_MeV)
        @test ambiguous.ambiguity_width_T_MeV == 1.0
        @test_throws ArgumentError Models.CEPResult(result_status=:resolved)
        resolved = Models.CEPResult(result_status=:resolved, T_cep_MeV=130, mu_cep_MeV=295)
        @test resolved.found
        @test resolved.T_cep_MeV == 130.0
    end

    @testset "find_cep 空 curves dict → 未找到" begin
        curves = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()
        result = Models.find_cep(curves)
        @test result isa Models.CEPResult
        @test result.found == false
    end

    @testset "_classify_s_curve 将临界弱 S 形从 invalid 中分离" begin
        curves = _load_reference_curves()
        mu_weak, rho_weak = _interpolate_curve(curves, 130.9375)
        cres = Models._classify_s_curve(mu_weak, rho_weak)

        @test cres.status == :weak_s_shape
        @test cres.reason == "weak_s_shape_no_sign_change"
        @test cres.mu_transition !== nothing
    end

    @testset "_classify_s_curve 不将明显非弱 no_sign_change 误判为 weak_s_shape" begin
        rho_vals = collect(0.0:0.05:1.2)
        mu_vals = [
            -0.0003281682876595267,
            0.04841824949440648,
            0.0936102611465675,
            0.1283662431809479,
            0.13861566900656133,
            0.10548202553599884,
            0.07231301412417683,
            0.1626490466370995,
            -0.125093482719776,
            -0.3379666334524407,
            -0.2773641542899543,
            -0.06974533828680642,
            0.22780715960196624,
            0.55209268408524,
            0.8480427451509811,
            1.0674944271492053,
            1.1759809660249494,
            1.1777855883649444,
            1.1218862556900857,
            1.0682844150768769,
            1.0498688154790619,
            1.066650699418709,
            1.1044049618930905,
            1.1509234887062252,
            1.2001534364461435,
        ]

        cres = Models._classify_s_curve(mu_vals, rho_vals)

        @test cres.status == :invalid
        @test cres.reason == "no_sign_change"
        @test cres.mu_transition === nothing
    end

    @testset "interpolate CEP 将 weak_s_shape disappearance 保留为 ambiguous 区间" begin
        curves = _load_reference_curves()
        mu_strong, rho_strong = _interpolate_curve(curves, 130.625)
        mu_weak, rho_weak = _interpolate_curve(curves, 130.9375)
        mu_none, rho_none = _interpolate_curve(curves, 131.0)

        seeded_curves = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}(
            130.0 => curves[130.0],
            130.625 => (mu_strong, rho_strong),
            130.9375 => (mu_weak, rho_weak),
            131.0 => (mu_none, rho_none),
        )

        cep = Models.find_cep(seeded_curves; tol=0.01, max_bisect_iter=12, strategy=:interpolate)

        @test !cep.found
        @test cep.result_status == :ambiguous
        @test isnan(cep.T_cep_MeV)
        @test cep.method == :three_state_interpolate
        @test cep.T_last_first_order_MeV == 130.625
        @test isnan(cep.T_first_monotone_MeV)
    end

    @testset "interpolate CEP 的 direct re-evaluate 不伪造单点" begin
        curves = _load_reference_curves()
        mu_strong, rho_strong = _interpolate_curve(curves, 130.625)
        mu_none, rho_none = _interpolate_curve(curves, 131.0)

        seeded_curves = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}(
            130.0 => curves[130.0],
            135.0 => curves[135.0],
        )

        evaluator = function (T_mid::Float64, level::Int)
            @test level == 0
            if isapprox(T_mid, 130.9375; atol=1e-8)
                return mu_strong, rho_strong
            elseif T_mid < 131.0
                return mu_strong, rho_strong
            end
            return mu_none, rho_none
        end

        cep_interp = Models.find_cep(
            copy(seeded_curves);
            tol=0.01,
            max_bisect_iter=12,
            strategy=:interpolate,
        )
        cep_direct_mid = Models.find_cep(
            copy(seeded_curves);
            tol=0.01,
            max_bisect_iter=12,
            strategy=:interpolate,
            evaluate_at_T=evaluator,
        )

        @test !cep_interp.found
        @test cep_interp.result_status == :ambiguous
        @test !cep_direct_mid.found
        @test cep_direct_mid.result_status == :ambiguous
        @test cep_direct_mid.T_last_first_order_MeV >= cep_interp.T_last_first_order_MeV
        @test cep_direct_mid.method == :three_state_interpolate
    end

    @testset "显式 monotone certificate 才能形成高侧证据" begin
        curves = _load_reference_curves()
        mu_s, rho_s = _interpolate_curve(curves, 130.625)
        rho_m = collect(0.0:0.05:1.2)
        mu_m = collect(range(250.0, 310.0; length=length(rho_m)))
        evidence = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}(
            130.625 => (mu_s, rho_s),
            140.0 => (mu_m, rho_m),
        )
        certificate = (T, cres) -> T >= 140.0 && cres.reason == "no_s_shape"
        cep = Models.find_cep(
            evidence;
            tol=0.5,
            max_bisect_iter=8,
            monotone_certificate=certificate,
        )
        @test cep.result_status == :ambiguous
        @test cep.T_last_first_order_MeV >= 130.625
        @test isfinite(cep.T_first_monotone_MeV)
        @test cep.ambiguity_width_T_MeV >= 0.0
    end

    @testset "unknown_budget 只停止收缩而不改写 ambiguous" begin
        curves = _load_reference_curves()
        mu_s, rho_s = _interpolate_curve(curves, 130.625)
        rho_m = collect(0.0:0.05:1.2)
        mu_m = collect(range(250.0, 310.0; length=length(rho_m)))
        evidence = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}(
            130.625 => (mu_s, rho_s),
            140.0 => (mu_m, rho_m),
        )
        cep = Models.find_cep(
            evidence;
            tol=0.5,
            max_bisect_iter=8,
            evaluate_at_T=(T, level) -> nothing,
            monotone_certificate=(T, cres) -> T >= 140.0 && cres.reason == "no_s_shape",
            strategy=:direct,
            unknown_budget=0,
        )
        @test cep.result_status == :ambiguous
        @test cep.unknown_count == 1
        @test occursin("unknown_budget_exceeded", something(cep.reason, ""))
        @test isnan(cep.T_cep_MeV)
    end
end
