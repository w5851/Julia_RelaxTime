## Plan: 步骤三 - 介子传播子模块MesonPropagator.jl实施计划

根据公式文档第1-3节和第8节的要求，本计划构建介子传播子计算模块，支持一般介子（π、K、σ_π、σ_K）和混合介子（η/η'、σ/σ'）的完整计算。核心设计遵循第8节的性能优化原则：将极化函数和耦合矩阵预计算后作为参数传入。

### Steps

1. **创建MesonPropagator.jl核心模块** [src/relaxtime/MesonPropagator.jl]
   - 实现`meson_propagator_simple`函数：计算一般介子传播子（π、K、σ_π、σ_K），函数签名为`meson_propagator_simple(meson_type::Symbol, channel::Symbol, K_coeffs::NamedTuple, Π::ComplexF64) -> ComplexF64`
   - **参数说明**：
     - `meson_type`: 介子类型(`:pi`, `:K`, `:sigma_pi`, `:sigma_K`)
     - `channel`: 通道类型(`:P`赝标量或`:S`标量)
     - `K_coeffs`: 预计算的K系数NamedTuple(通过`EffectiveCouplings.calculate_effective_couplings`获取)
     - `Π`: 预计算的极化函数(ComplexF64)
   - **K系数自动选择**：函数内部根据`meson_type`和`channel`自动选择正确的K系数，降低调用者使用难度
   - **关键映射**：赝标量P（π、K）使用K⁺系数，标量S（σ_π、σ_K）使用K⁻系数
   - **K介子极化函数**：使用Π_{us}而非Π_{uu}
   - 实现`meson_propagator_mixed`函数：计算混合介子传播子（η/η'、σ/σ'），接受散射过程参数(q1::Symbol, q2::Symbol, q3::Symbol, q4::Symbol)表示散射过程q1+q2→q3+q4，以及散射道channel::Symbol(:t/:s/:u)，使用矩阵乘法形式`2det(K)/det(M) × J^T M J'`，返回ComplexF64
   - **散射过程约束**：只有q2和q4可能是反夸克(:ubar/:dbar/:sbar)，且若有反夸克则q2和q4同时为反夸克
   - **场算符映射**：函数内部根据channel自动将q1,q2,q3,q4映射到ψ,ψ̄,ψ',ψ̄'场算符，调用者无需处理复杂的场算符对应关系
   - **关键映射**：赝标量P（η/η'）使用K⁺系数，标量S（σ/σ'）使用K⁻系数
   - 实现`extract_flavor`辅助函数：从Symbol中提取味类型(去除bar标记)，如:ubar→:u
   - 实现`get_quark_wavefunction`辅助函数：根据味类型(:u/:d/:s)和是否为bar返回对应的列向量或行向量波函数
   - 在Constants_PNJL.jl中添加常量定义(使用ASCII bar命名，如ψbar_u而非ψ̄_u)：
     - Gell-Mann矩阵：`λ₀`和`λ₈`（3×3矩阵）
     - 夸克味波函数：`ψ_u`、`ψ_d`、`ψ_s`（列向量）和`ψbar_u`、`ψbar_d`、`ψbar_s`（行向量/1×3矩阵）
     - **所有新常量必须在Constants_PNJL.jl中export**
   - 实现`calculate_coupling_matrix`辅助函数:根据极化函数(Π_uu::ComplexF64, Π_ss::ComplexF64)、预计算的K系数(K_coeffs::NamedTuple)和通道类型channel::Symbol(:P或:S)计算复数耦合矩阵M(2×2 ComplexF64矩阵)
   - **采用方案A**:K系数通过EffectiveCouplings.calculate_effective_couplings预先计算并作为参数传入,函数内部根据channel自动选择对应的K系数(:P通道用K⁺系数,:S通道用K⁻系数),符合批量计算时复用K系数的性能优化原则
   - **M矩阵对称性**：利用M₀₈ = M₈₀减少计算
   - **M₀₈系数修正**：使用4√2/3（即4/3·√2），不是4/(3√2)
   - **det_M计算**：直接使用Julia标准库的det()函数计算M矩阵行列式，无需额外定义det_M函数(det_K已在EffectiveCouplings.jl中定义)
   - 添加完整文档字符串，说明参数单位（fm⁻¹、fm²）、各向异性修正（通过Π传入）和使用示例

2. **编写API文档** [api/MesonPropagator.md]
   - 遵循`api/EffectiveCouplings.md`格式：模块概述、依赖关系、单位约定、API参考（含参数表和物理意义）
   - **K系数预计算说明**：所有函数都假设K系数已通过`EffectiveCouplings.calculate_effective_couplings`预先计算,批量计算时可复用同一组K系数以提升性能
   - **K系数自动选择逻辑**：详细说明`meson_propagator_simple`函数如何根据`meson_type`和`channel`自动选择正确的K系数，减轻调用者负担
   - **明确介子类型映射表**（关键修正）：
     | 介子 | 通道 | 夸克组合 | 极化函数 | K系数 |
     |------|------|----------|---------|-------|
     | π | P | u̅u/d̅d | Π_{uu}^P | K₃⁺ |
     | K | P | u̅s/d̅s | Π_{us}^P | K₄⁺ |
     | σ_π | S | u̅u/d̅d | Π_{uu}^S | K₃⁻ |
     | σ_K | S | u̅s/d̅s | Π_{us}^S | K₄⁻ |
     | η/η' | P | (u̅u±s̅s)混合 | Π_{uu}^P, Π_{ss}^P | K₀⁺, K₈⁺, K₀₈⁺ |
     | σ/σ' | S | (u̅u±s̅s)混合 | Π_{uu}^S, Π_{ss}^S | K₀⁻, K₈⁻, K₀₈⁻ |
   - 说明Gell-Mann矩阵（λ₀、λ₈）和味波函数的预定义常量设计
   - 明确所有传播子函数返回ComplexF64类型
   - 提供完整使用示例：从预计算极化函数（ComplexF64）到调用传播子函数的端到端代码
   - 补充性能优化说明：批量计算时一次性预计算Π和M以减少重复开销

3. **开发单元测试** [test/test_meson_propagator.jl]
   - 基本功能测试：验证`meson_propagator_simple`和`meson_propagator_mixed`返回ComplexF64、单位正确（fm²）
   - **通道-K系数映射测试**：验证π和K使用K⁺，σ_π和σ_K使用K⁻，η/η'使用K⁺，σ/σ'使用K⁻
   - **K介子极化函数测试**：验证K介子使用Π_{us}而非Π_{uu}，对比两者差异
   - 物理约束测试：检查k→0极限下传播子实部发散（质量壳条件）、虚部为零（无衰变）
   - 混合介子测试：验证η/η'混合时矩阵乘法计算的正确性，与手动展开结果对比（相对误差<1e-10）
   - **M₀₈系数测试**：验证使用4√2/3而非4/(3√2)，对比两者数值差异
   - 各向异性修正测试：对比ξ=0和ξ=0.3时传播子的变化（预期修正量级~10-30%，此为经验值，无理论依据，通过测试确定实际修正量级即可）
   - 夸克味配置测试：测试u̅u、u̅s、s̅s不同夸克组合的传播子计算（涵盖所有波函数常量）
   - 参考`test/test_effective_couplings.jl`结构：使用@testset组织测试、预生成积分节点、输出测试总结

4. **创建测试总结文档** [test/test_meson_propagator_summary.md]
   - 记录测试通过率：xx/xx测试全部通过
   - 数值精度验证：矩阵乘法展开对比的相对误差结果
   - 物理合理性检验：传播子量级、各向异性修正幅度、不同介子类型的差异
   - 性能基准：预计算优化带来的速度提升（预期~3-5倍）
   - 已知限制：衰变宽度计算、Mott相变检测等高级功能留待后续步骤

### Further Considerations

1. **极化函数预计算策略** - ✅已确认：本模块接收ComplexF64类型的预计算Π值。网格插值实现推迟至步骤6性能优化阶段。
   
2. **Gell-Mann矩阵的完整性** - ✅已确认：当前只需实现λ₀和λ₈。其他矩阵（λ₃到λ₇）在公式文档中标记为"供参考，本模块暂不需要"。按需添加，避免过早优化。

3. **介子质量提取功能** - ✅已确认：本步骤专注于传播子计算。质量提取功能留待后续步骤或作为独立工具函数实现。

### 关键公式修正总结（已反映在上述步骤中）

根据公式文档修正，以下是实现时必须注意的5个关键点：

1. **通道-K系数映射**（与直觉相反）：
   - ❌ 错误：赝标量P用K⁻，标量S用K⁺
   - ✅ 正确：赝标量P用K⁺，标量S用K⁻

2. **K介子极化函数**：
   - ❌ 错误：使用Π_{uu}
   - ✅ 正确：使用Π_{us}

3. **M₀₈系数**：
   - ❌ 错误：4/(3√2) ≈ 0.9428
   - ✅ 正确：4√2/3 ≈ 1.8856

4. **数据类型**：
   - 所有传播子和M矩阵使用ComplexF64
   - 极化函数Π作为ComplexF64传入

5. **M矩阵对称性**：
   - M₀₈ = M₈₀（可利用此减少计算量）
