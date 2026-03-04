"""
寻找有根的测试用例
"""

include("../../../src/relaxtime/OneLoopIntegralsAniso.jl")
using .OneLoopIntegralsCorrection: tilde_B0_correction_k_positive, STRATEGY_HYBRID

# 扫描参数空间寻找有根的情况
function find_root_cases()
    println("扫描参数空间寻找有根的情况...")
    
    root_cases = []
    
    for λ in [0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0]
        for k in [0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0]
            for m in [0.1, 0.3, 0.5]
                for m_prime in [0.1, 0.3, 0.5]
                    μ, T, Φ, Φbar, ξ = 0.3, 0.15, 0.3, 0.3, 0.1
                    
                    try
                        _, _, diag = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
                            strategy=STRATEGY_HYBRID, diagnostics=true)
                        
                        if diag.n_roots > 0
                            push!(root_cases, (λ=λ, k=k, m=m, m_prime=m_prime, n_roots=diag.n_roots, 
                                n_intervals=diag.n_intervals, roots=diag.roots))
                        end
                    catch e
                        # 忽略错误
                    end
                end
            end
        end
    end
    
    println("\n找到 $(length(root_cases)) 个有根的情况:")
    println("-" ^ 70)
    
    for (i, case) in enumerate(root_cases[1:min(20, length(root_cases))])
        println("Case $i: λ=$(case.λ), k=$(case.k), m=$(case.m), m'=$(case.m_prime)")
        println("        根数=$(case.n_roots), 区间数=$(case.n_intervals)")
        println("        根位置: $(case.roots)")
        println()
    end
    
    return root_cases
end

root_cases = find_root_cases()
