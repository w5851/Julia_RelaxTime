# 寻找 k0²-u = -5.88 的来源

println("="^70)
println("寻找 k0²-u = -5.881889525074166 的来源")
println("="^70)

# 测试参数
s = 20.0
t = -0.2

m_u = 0.3
m_d = 0.3
m_s = 0.5

# 各种过程的u值
println("\n## 各种过程的u值")
println("-"^70)

processes = [
    ("uu→uu", 4*m_u^2),
    ("ss→ss", 4*m_s^2),
    ("ud→ud", 2*m_u^2 + 2*m_d^2),
    ("us→us", 2*m_u^2 + 2*m_s^2),
]

for (name, sum_m2) in processes
    u = sum_m2 - s - t
    println("$name: Σm²=$(sum_m2), u=$(u)")
end

# 检查是否有 k0²-u = -5.88
target = -5.881889525074166

println("\n## 寻找满足 k0²-u = $target 的情况")
println("-"^70)

# 假设k0=0（同种粒子散射）
for (name, sum_m2) in processes
    u = sum_m2 - s - t
    k0_squared_minus_u = 0.0 - u
    if abs(k0_squared_minus_u - target) < 0.01
        println("✓ 找到！$name: k0²-u = 0-($u) = $k0_squared_minus_u")
    end
end

# 检查其他可能：不同的Mandelstam变量组合
println("\n## 检查其他散射道")
println("-"^70)

# ud→ud 过程
println("\nud→ud 过程:")
m1, m2, m3, m4 = m_u, m_d, m_u, m_d
sum_m2 = m1^2 + m2^2 + m3^2 + m4^2
println("  Σm² = $sum_m2")

# 计算三个Mandelstam变量
s_val = s
t_val = t
u_val = sum_m2 - s_val - t_val

println("  s = $s_val")
println("  t = $t_val")  
println("  u = $u_val")

# 对于ud→ud，检查t道和u道
# t道: k0 = |E1 - E3| = |E_u - E_u| = 0
k0_t = 0.0
k0_squared_minus_t = k0_t^2 - t_val
println("\n  t道: k0²-t = $(k0_squared_minus_t)")

# u道: k0 = |E1 - E4| = |E_u - E_d| ≈ 0 (质量相同)
k0_u = 0.0
k0_squared_minus_u = k0_u^2 - u_val
println("  u道: k0²-u = $(k0_squared_minus_u)")

if abs(k0_squared_minus_u - target) < 0.01
    println("  ✓✓✓ 找到匹配！ud→ud的u道 k0²-u = $k0_squared_minus_u")
end

# us→us 过程
println("\nus→us 过程:")
m1, m2, m3, m4 = m_u, m_s, m_u, m_s
sum_m2 = m1^2 + m2^2 + m3^2 + m4^2
println("  Σm² = $sum_m2")

s_val = s
t_val = t
u_val = sum_m2 - s_val - t_val

println("  s = $s_val")
println("  t = $t_val")  
println("  u = $u_val")

# u道: k0 = |E_u - E_s| ≠ 0 (质量不同!)
# 但在质心系中仍然可能接近0
k0_u_approx = 0.0  # 近似
k0_squared_minus_u = k0_u_approx^2 - u_val
println("  u道: k0²-u ≈ $(k0_squared_minus_u)")

if abs(k0_squared_minus_u - target) < 0.01
    println("  ✓✓✓ 找到匹配！us→us的u道 k0²-u ≈ $k0_squared_minus_u")
end

# 测试中的实际调用顺序
println("\n## 测试中的实际调用顺序")
println("-"^70)
println("1. uu→uu: 成功")
println("2. ss→ss: 触发警告 (k0²-u = -5.88)")
println("3. ud→ud: ...")
println("4. us→us: ...")

println("\n猜测：警告可能来自ss→ss过程中计算传播子时的内部调用")
println("      或者来自前一个过程的延迟警告输出")

println("\n" * "="^70)
