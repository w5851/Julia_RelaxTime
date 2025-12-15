## Plan: 在 RTA 下实现输运系数模块

基于仓库现状，实现 η/ζ/σ 的关键依赖（弛豫时间 τ、PNJL 平衡分布函数、热力学导数、数值积分工具）基本齐全；但存在两个高优先级缺口会影响结果可信度：各向异性分布模块中疑似把动量 p 当作能量 E 传入分布函数、以及电导率需要的味电荷常量缺失。计划优先把输入约定（单位、变量）统一清楚，再新增 src/relaxtime 下的输运系数模块，并同步补齐 API 文档与单元测试。

### Steps
1. 盘点并固定依赖接口：对接 [src/relaxtime/RelaxationTime.jl](src/relaxtime/RelaxationTime.jl) 的 `relaxation_times` 输出 `tau`，以及 [src/pnjl/analysis/ThermoDerivatives.jl](src/pnjl/analysis/ThermoDerivatives.jl) 的 `bulk_derivative_coeffs`/`thermo_derivatives` 输出体粘滞所需导数。
2. 核对分布函数自变量与单位：统一用 [src/pnjl/PNJLQuarkDistributions.jl](src/pnjl/PNJLQuarkDistributions.jl) 的 “以能量 E 为自变量”的分布入口，并检查/修正 [src/pnjl/PNJLQuarkDistributions_Aniso.jl](src/pnjl/PNJLQuarkDistributions_Aniso.jl) 中 `distribution_with_anisotropy` 的 ξ=0 退化分支是否误用 p→E。
A:这里应该不是误用，因为动量各向异性时分布函数的值应该同时依赖p和x=cosθ，不能直接以能量E作为参数，因此和各向同性下的参数接口不同
3. 设计输运模块 API 与数据结构：在 src/relaxtime 新增 `TransportCoefficients` 模块，提供 `shear_viscosity`/`bulk_viscosity`/`electric_conductivity`/`transport_coefficients`，输入以 `quark_params`,`thermo_params`,`tau`(NamedTuple) 为主，显式传入 `degeneracy` 与 `charges`。
4. 复用现有 Gauss-Legendre 网格实现动量积分：沿用 [src/integration/GaussLegendre.jl](src/integration/GaussLegendre.jl) 的 `gauleg` 与默认节点/权重，封装通用积分核以复用到 η/ζ/σ（共享 p 网格、避免重复分配）。
5. 同步补齐文档与测试：新增 docs/api/TransportCoefficients.md（或拆分三份），并在 tests/unit 新增输运系数单测，主要断言：结果有限/非负、τ→0 时系数→0、u 与 d 在同参数下对称、σ 对电荷平方缩放关系成立。

### Further Considerations
1. 单位口径：输运模块对外接受 `T, μ` 用 MeV 还是 fm⁻¹？建议选一种并提供桥接（当前 τ 链多用 fm⁻¹，热力学导数链多用 MeV）。
使用fm⁻¹ 作为标准单位；PNJL中热力学量其实也应该使用fm⁻¹ 作为单位
2. 体粘滞里的 n：采用净重子数密度还是总粒子数密度？仓库里 [src/pnjl/analysis/ThermoDerivatives.jl](src/pnjl/analysis/ThermoDerivatives.jl) 目前返回 `rho_total=sum(rho_vec)/3`，需明确与公式中的 n 对应关系。
采用净重子数密度，仓库代码中的rho_total计算的就是净重子数密度，你可以为这个变量赋予更清晰的命名。
3. τ 的建模：先实现“每味常数 τ_i”版本（与现有 `relaxation_times` 一致），还是需要扩展到 τ(p)？两者 API 会不同，建议先做常数版本以匹配现有基础设施。
τ和p应该无关，直接使用现有代码即可
这份计划你希望按“最小可用（常数 τ_i + 各向同性 ξ=0）”先落地，还是一开始就把各向异性 ξ≠0 的分布也接进去？
