# 分析不同区间所需的最小节点数
using Printf
include(joinpath(@__DIR__, "../../../src/relaxtime/OneLoopIntegralsAniso.jl"))
include(joinpath(@__DIR__, "../../../src/relaxtime/OneLoopIntegrals.jl"))
using QuadGK: quadgk

println("=" ^ 70)
println("分析不同区间所需的最小节点数")
println("=" ^ 70)

# 测试参数
test_cases = [
    (name="有根1", λ=-1.0, k=0.01, m=0.3, m_prime=0.3),
    (name="有根2", λ=-1.0, k=0.1, m=0.3, m_prime=0.3),
    (name="有根3", λ=-0.8, k=0.02, m=0.25, m_prime=0.25),
    (name="无根1", λ=-0.5, k=0.05, m=0.3, m_prime=0.3),
    (name="无根2", λ=0.5, k=0.1, m=0.3, m_prime=0.35),
]

ξ = -0.2; T = 0.15; μ = 0.0; Φ = 0.0; Φbar = 0.0
target_err = 1e-4  # 目标精度

println("\n1. 有根情况：各区间所需最小节点数 (目标误差 < 1e-4)")
println("-" ^ 70)

for tc in filter(t -> startswith(t.name, "有根"), test_cases)
    Emin = tc.m
    Emax = OneLoopIntegrals.energy_cutoff(tc.m)
    integrand(E) = OneLoopIntegralsCorrection.real_integrand_k_positive(
        :quark, tc.λ, tc.k, tc.m, tc.m_prime, E, ξ, T, μ, Φ, Φbar)
    
    roots = OneLoopIntegralsCorrection.find_roots_AB(tc.λ, tc.k, tc.m, tc.m_prime, Emin, Emax)
    intervals = OneLoopIntegralsCorrection.build_intervals_from_roots(roots, Emin, Emax)
    
    println("\n$(tc.name): $(length(roots)) 个根, $(length(intervals)) 个区间")
    
    for (idx, (a, b)) in enumerate(intervals)
        # 高精度参考值
        ref, _ = quadgk(integrand, a, b; rtol=1e-12)
        
        # 确定奇点位置
        n_intervals = length(intervals)
        sing_pos = if n_intervals == 1
            isempty(roots) ? OneLoopIntegralsCorrection.SING_NONE : OneLoopIntegralsCorrection.SING_BOTH
        elseif idx == 1
            OneLoopIntegralsCorrection.SING_RIGHT
        elseif idx == n_intervals
            OneLoopIntegralsCorrection.SING_LEFT
        else
            OneLoopIntegralsCorrection.SING_BOTH
        end
        
        # 测试不同节点数
        min_n = 0
        for n in [4, 8, 12, 16, 20, 24, 28, 32]
            xs, wx = OneLoopIntegralsCorrection.hybrid_nodes(a, b, n, sing_pos)
            result = sum(wx .* integrand.(xs))
            err = abs(ref) > 1e-15 ? abs((result - ref) / ref) : abs(result - ref)
            if err < target_err
                min_n = n
                break
            end
        end
        if min_n == 0
            min_n = 32  # 需要 32 或更多
        end
        
        interval_len = b - a
        println(@sprintf("  区间 %d [%.4f, %.4f] (长度=%.4f, %s): 最小 n=%d", 
            idx, a, b, interval_len, sing_pos, min_n))
    end
end

println("\n" * "=" ^ 70)
println("2. 无根情况：所需最小节点数")
println("-" ^ 70)

for tc in filter(t -> startswith(t.name, "无根"), test_cases)
    Emin = tc.m
    Emax = OneLoopIntegrals.energy_cutoff(tc.m)
    integrand(E) = OneLoopIntegralsCorrection.real_integrand_k_positive(
        :quark, tc.λ, tc.k, tc.m, tc.m_prime, E, ξ, T, μ, Φ, Φbar)
    
    ref, _ = quadgk(integrand, Emin, Emax; rtol=1e-12)
    
    min_n = 0
    for n in [4, 8, 12, 16, 20, 24, 28, 32]
        xs, wx = OneLoopIntegralsCorrection.power_left_nodes(Emin, Emax, n)
        result = sum(wx .* integrand.(xs))
        err = abs((result - ref) / ref)
        if err < target_err
            min_n = n
            break
        end
    end
    if min_n == 0
        min_n = 32
    end
    
    println(@sprintf("  %s: 最小 n=%d", tc.name, min_n))
end

println("\n" * "=" ^ 70)
println("3. 区间长度与所需节点数的关系")
println("-" ^ 70)

# 收集所有区间数据
interval_data = []

for tc in test_cases
    Emin = tc.m
    Emax = OneLoopIntegrals.energy_cutoff(tc.m)
    integrand(E) = OneLoopIntegralsCorrection.real_integrand_k_positive(
        :quark, tc.λ, tc.k, tc.m, tc.m_prime, E, ξ, T, μ, Φ, Φbar)
    
    roots = OneLoopIntegralsCorrection.find_roots_AB(tc.λ, tc.k, tc.m, tc.m_prime, Emin, Emax)
    intervals = OneLoopIntegralsCorrection.build_intervals_from_roots(roots, Emin, Emax)
    
    for (idx, (a, b)) in enumerate(intervals)
        ref, _ = quadgk(integrand, a, b; rtol=1e-12)
        
        n_intervals = length(intervals)
        sing_pos = if n_intervals == 1
            isempty(roots) ? OneLoopIntegralsCorrection.SING_NONE : OneLoopIntegralsCorrection.SING_BOTH
        elseif idx == 1
            OneLoopIntegralsCorrection.SING_RIGHT
        elseif idx == n_intervals
            OneLoopIntegralsCorrection.SING_LEFT
        else
            OneLoopIntegralsCorrection.SING_BOTH
        end
        
        min_n = 32
        for n in [4, 8, 12, 16, 20, 24, 28, 32]
            xs, wx = OneLoopIntegralsCorrection.hybrid_nodes(a, b, n, sing_pos)
            result = sum(wx .* integrand.(xs))
            err = abs(ref) > 1e-15 ? abs((result - ref) / ref) : abs(result - ref)
            if err < target_err
                min_n = n
                break
            end
        end
        
        push!(interval_data, (
            case=tc.name, 
            idx=idx, 
            len=b-a, 
            sing=sing_pos, 
            min_n=min_n,
            is_first=(idx==1),
            is_last=(idx==n_intervals),
            is_middle=(idx>1 && idx<n_intervals)
        ))
    end
end

# 按区间类型分组统计
println("\n按区间类型统计:")
for (desc, filter_fn) in [
    ("第一个区间 (含左端点 E=m)", d -> d.is_first && length(OneLoopIntegralsCorrection.find_roots_AB(-1.0, 0.01, 0.3, 0.3, 0.3, 3.0)) > 0),
    ("中间区间 (两端都有根)", d -> d.is_middle),
    ("最后一个区间", d -> d.is_last && !d.is_first),
    ("无根情况 (单区间)", d -> d.is_first && d.is_last),
]
    matching = filter(filter_fn, interval_data)
    if !isempty(matching)
        avg_n = sum(d.min_n for d in matching) / length(matching)
        max_n = maximum(d.min_n for d in matching)
        println(@sprintf("  %s: 平均 n=%.1f, 最大 n=%d", desc, avg_n, max_n))
    end
end

println("\n" * "=" ^ 70)
println("4. 建议的节点分配策略")
println("=" ^ 70)
println("""
基于分析结果，建议：
- 第一个区间（含 E=m 奇点）: 需要较多节点处理 p→0 奇点
- 中间区间（两端都有根）: 通常很短，需要较少节点
- 最后一个区间: 通常最长，但被积函数平滑，需要中等节点
- 无根情况: power_left 变换效果极好，可用较少节点
""")
