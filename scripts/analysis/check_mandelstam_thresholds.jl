# 计算各种过程的最小s值（阈值）

m_u = 300.0 / 197.327  # ≈ 1.52 fm⁻¹
m_s = 500.0 / 197.327  # ≈ 2.53 fm⁻¹

println("="^70)
println("Mandelstam变量s的阈值条件")
println("="^70)

println("\n裸质量:")
println("  m_u = $m_u fm⁻¹")
println("  m_s = $m_s fm⁻¹")

println("\n各种qq散射过程的阈值s_min = Σm²:")
processes = [
    ("uu→uu", 4*m_u^2),
    ("ud→ud", 2*m_u^2 + 2*m_u^2),  # 等价于4m_u²
    ("us→us", 2*m_u^2 + 2*m_s^2),
    ("dd→dd", 4*m_u^2),             # d和u质量相同
    ("ds→ds", 2*m_u^2 + 2*m_s^2),
    ("ss→ss", 4*m_s^2),
]

global max_threshold = 0.0
for (name, threshold) in processes
    println("  $name: s_min = $threshold fm⁻²")
    global max_threshold = max(max_threshold, threshold)
end

println("\n所有qq过程的最大阈值: $(max_threshold) fm⁻²")
println("建议的测试参数: s = $(ceil(max_threshold) + 5) fm⁻² (略大于阈值)")

# 检查当前测试使用的s=20.0是否足够
s_test = 20.0
println("\n当前测试使用 s = $s_test fm⁻²")
if s_test >= max_threshold
    println("✓ 满足所有阈值条件")
else
    println("✗ 不满足阈值！最大阈值为 $(max_threshold) fm⁻²")
    println("  ss→ss 需要 s ≥ $(4*m_s^2) fm⁻²")
end

# 使用s=20.0时各过程的u值
println("\n使用 s=20.0, t=-0.2 时各过程的 u 值:")
t = -0.2
for (name, sum_m2) in processes
    u = sum_m2 - s_test - t
    status = u < 0 ? "✓" : "✗"
    println("  $name: u = $u fm⁻² $status")
end
