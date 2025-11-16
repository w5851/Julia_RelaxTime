## Plan: 步骤三 - 介子传播子模块MesonPropagator.jl实施计划

根据公式文档第1-3节和第8节的要求，本计划构建介子传播子计算模块，支持一般介子（π、K、σ_π、σ_K）和混合介子（η/η'、σ/σ'）的完整计算。核心设计遵循第8节的性能优化原则：将极化函数和耦合矩阵预计算后作为参数传入。

### Steps

1. **创建MesonPropagator.jl核心模块** [src/relaxtime/MesonPropagator.jl]
   - 实现`meson_propagator_simple`函数：计算一般介子传播子（π、K、σ_π、σ_K），函数签名为`meson_propagator_simple(meson_type::Symbol, K_coeffs::NamedTuple, Π::ComplexF64) -> ComplexF64`
   - **参数说明**：
     - `meson_type`: 介子类型(`:pi`, `:K`, `:sigma_pi`, `:sigma_K`)，介子类型已隐含通道信息
     - `K_coeffs`: 预计算的K系数NamedTuple(通过`EffectiveCouplings.calculate_effective_couplings`获取)
     - `Π`: 预计算的极化函数(ComplexF64)
   - **K系数自动选择**：函数内部根据`meson_type`自动选择正确的K系数，降低调用者使用难度
     - `:pi` → `K123_plus` (π是赝标量P通道，用K⁺)
     - `:K` → `K4567_plus` (K是赝标量P通道，用K⁺)
     - `:sigma_pi` → `K123_minus` (σ_π是标量S通道，用K⁻)
     - `:sigma_K` → `K4567_minus` (σ_K是标量S通道，用K⁻)
   - **关键映射**：赝标量P（π、K）使用K⁺系数，标量S（σ_π、σ_K）使用K⁻系数
   - **K介子极化函数**：使用Π_{us}而非Π_{uu}
   - 实现`meson_propagator_mixed`函数：计算混合介子传播子（η/η'、σ/σ'），函数签名为`meson_propagator_mixed(det_K::Float64, M_matrix::Matrix{ComplexF64}, q1::Symbol, q2::Symbol, q3::Symbol, q4::Symbol, channel::Symbol) -> ComplexF64`，接受预计算的det_K和M矩阵、散射过程参数(q1, q2, q3, q4)表示散射过程q1+q2→q3+q4，以及散射道channel::Symbol(:t/:s/:u)，使用矩阵乘法形式`2det(K)/det(M) × J^T M J'`，返回ComplexF64
   - **参数说明**：
     - `det_K`: 预计算的耦合矩阵行列式(通过`EffectiveCouplings.coupling_matrix_determinant`获取，也称为det_K)
     - `M_matrix`: 预计算的2×2复数耦合矩阵(通过`calculate_coupling_matrix`辅助函数获取)
     - `q1, q2, q3, q4`: 散射过程中的夸克/反夸克类型，可以是`:u/:d/:s`或`:ubar/:dbar/:sbar`
     - `channel`: 散射道类型，`:t`(t道)、`:s`(s道)或`:u`(u道)
   - **散射过程的物理意义**：这是2→2散射过程q1+q2→q3+q4，其中：
     - q1和q3**必定是夸克**(`:u`, `:d`, `:s`)
     - q2和q4**要么全为夸克，要么全为反夸克**(`:ubar`, `:dbar`, `:sbar`)
     - t道和u道：纯夸克散射过程
     - s道：夸克-反夸克湮灭过程
   - **散射道到场算符的映射表**：
     | 散射道 | 物理过程 | ψ | ψ̄ | ψ' | ψ̄' | 介子四动量 |
     |-------|---------|---|---|----|----|----------|
     | t道 | q₁+q₂→q₃+q₄(纯夸克) | q1 | q3 | q2 | q4 | q=p₁-p₃ |
     | s道 | q₁+q̄₂→q₃+q̄₄(湮灭) | q1 | q2 | q3 | q4 | P=p₁+p₂ |
     | u道 | q₁+q₂→q₃+q₄(交叉) | q1 | q4 | q2 | q3 | q=p₁-p₄ |
   - **介子四动量说明**：介子能量k₀和动量|k⃗|从四动量q或P中提取(k₀²-|k⃗|²=q²)，但在传播子函数中不直接使用，仅通过预计算的极化函数Π间接影响传播子
   - **场算符映射说明**：函数内部根据channel和散射粒子类型自动构建流算符向量J和J'，调用者只需提供散射过程的粒子标识即可，无需手动处理复杂的场算符对应关系
   - **场算符的产生/湮灭性质**(关键物理概念)：
     - 对于**夸克**：$\psi$是**湮灭算符**，$\bar{\psi}$是**产生算符**
     - 对于**反夸克**：$\psi$是**产生算符**，$\bar{\psi}$是**湮灭算符**
     - 判断依据：费曼图中箭头指向顶点为湮灭，离开顶点为产生
     - t道：每个顶点连接一个入射和一个出射粒子(直接传递)
     - s道：第一个顶点湮灭两个入射粒子，第二个顶点产生两个出射粒子(湮灭-产生)
     - u道：连接交叉的粒子对(交叉传递)
   - **旋量处理**：公式中的Dirac旋量指标通过矩阵乘法$\bar{\psi}\lambda\psi$自动收缩，只保留味空间的矩阵运算
   - **关键映射**：赝标量P（η/η'）使用K⁺系数，标量S（σ/σ'）使用K⁻系数
   - 实现`extract_flavor`辅助函数：从Symbol中提取味类型(去除bar标记)，如`:ubar→:u`，返回`(:u, true)`表示u味的反夸克
   - 实现`get_quark_wavefunction`辅助函数：根据味类型和is_bar标志返回对应波函数
     - `is_bar=false`: 返回列向量ψ（如ψ_u=[1,0,0]）
     - `is_bar=true`: 返回行向量ψ̄（如ψbar_u=[1 0 0]，即1×3矩阵）
   - 实现`calculate_current_vector`辅助函数：根据散射道映射表自动选择场算符，计算流算符向量
     - J = [ψ̄·λ₀·ψ, ψ̄·λ₈·ψ]（2×1列向量）
     - 矩阵乘法从左到右依次进行：`ψbar * λ * ψ`（1×3矩阵 × 3×3矩阵 × 3×1向量 = 标量）
   - **Constants_PNJL.jl中的常量**：以下常量已在`Constants_PNJL.jl`中定义并export(使用ASCII bar命名，如ψbar_u而非ψ̄_u)：
     - Gell-Mann矩阵：`λ₀`和`λ₈`（3×3矩阵）
     - 夸克味波函数：`ψ_u`、`ψ_d`、`ψ_s`（列向量）和`ψbar_u`、`ψbar_d`、`ψbar_s`（行向量/1×3矩阵）
   - 实现`calculate_coupling_matrix`辅助函数:根据极化函数(Π_uu::ComplexF64, Π_ss::ComplexF64)、预计算的K系数(K_coeffs::NamedTuple)和通道类型channel::Symbol(:P或:S)计算复数耦合矩阵M(2×2 ComplexF64矩阵)
   - **采用方案A**:K系数通过EffectiveCouplings.calculate_effective_couplings预先计算并作为参数传入,函数内部根据channel自动选择对应的K系数(:P通道用K⁺系数,:S通道用K⁻系数),符合批量计算时复用K系数的性能优化原则
   - **M矩阵对称性**：利用M₀₈ = M₈₀减少计算
   - **M₀₈系数修正**：使用`4*sqrt(2)/3`（数值约1.8856），不是4/(3√2)（数值约0.9428）
   - **det_M计算**：直接使用Julia标准库的det()函数计算M矩阵行列式，无需额外定义det_M函数(det_K已在EffectiveCouplings.jl中定义)
   - 添加完整文档字符串，说明参数单位（fm⁻¹、fm²）、各向异性修正（通过Π传入）和使用示例

2. **编写API文档** [api/MesonPropagator.md]
   - 遵循`api/EffectiveCouplings.md`格式：模块概述、依赖关系、单位约定、API参考（含参数表和物理意义）
   - **K系数预计算说明**：所有函数都假设K系数已通过`EffectiveCouplings.calculate_effective_couplings`预先计算,批量计算时可复用同一组K系数以提升性能
   - **det_K函数说明**：明确`EffectiveCouplings.coupling_matrix_determinant`函数用于预计算det_K（耦合矩阵行列式），文档中det_K和coupling_matrix_determinant指代同一概念
   - **K系数命名对照表**：添加物理符号与代码变量名的对应关系
     | 物理符号 | 代码变量名 | 说明 | 对应介子 |
     |---------|-----------|------|----------|
     | K₀± | `K0_plus/minus` | 单态通道 | η'/σ |
     | K₁±=K₂±=K₃± | `K123_plus/minus` | π通道(三个分量相等) | π, a₀ |
     | K₄±=K₅±=K₆±=K₇± | `K4567_plus/minus` | K通道(四个分量相等) | K, K₀* |
     | K₈± | `K8_plus/minus` | 八重态通道 | η₈, f₀ |
     | K₀₈± | `K08_plus/minus` | 混合通道 | η-η'混合 |
   - **K系数自动选择逻辑**：详细说明`meson_propagator_simple`函数如何根据`meson_type`自动选择正确的K系数，减轻调用者负担
   - **明确介子类型映射表**（关键修正，使用代码变量名）：
     | 介子 | 通道 | 夸克组合 | 极化函数 | K系数变量名 |
     |------|------|----------|---------|-------------|
     | π | P | u̅u/d̅d | Π_{uu}^P | `K123_plus` |
     | K | P | u̅s/d̅s | Π_{us}^P | `K4567_plus` |
     | σ_π | S | u̅u/d̅d | Π_{uu}^S | `K123_minus` |
     | σ_K | S | u̅s/d̅s | Π_{us}^S | `K4567_minus` |
     | η/η' | P | (u̅u±s̅s)混合 | Π_{uu}^P, Π_{ss}^P | `K0_plus`, `K8_plus`, `K08_plus` |
     | σ/σ' | S | (u̅u±s̅s)混合 | Π_{uu}^S, Π_{ss}^S | `K0_minus`, `K8_minus`, `K08_minus` |
   - 说明Gell-Mann矩阵（λ₀、λ₈）和味波函数的预定义常量设计
   - **旋量处理说明**：明确流算符向量中的Dirac旋量指标已通过矩阵乘法$\bar{\psi}\lambda\psi$自动收缩，本模块只处理味空间的3×3矩阵运算
   - 明确所有传播子函数返回ComplexF64类型
   - 提供完整使用示例：从预计算极化函数（ComplexF64）到调用传播子函数的端到端代码
   - 补充性能优化说明：批量计算时一次性预计算Π和M以减少重复开销
   - **各向异性修正的设计哲学**（重要设计决策）：
     - **传播子模块本身不直接依赖ξ参数**，这是有意的模块化设计
     - 各向异性通过预计算的极化函数值间接影响传播子
     - 调用链：`PolarizationAniso.polarization_aniso(ξ, ...)` → 返回Π(ComplexF64) → `MesonPropagator.meson_propagator(Π)`
     - 参考`src/relaxtime/PolarizationAniso.jl`：该模块实现了包含B0_correction项的各向异性极化函数
     - 使用时在调用链最前端选择`Polarization.jl`(各向同性)或`PolarizationAniso.jl`(各向异性)，传播子函数无需修改
     - 这种设计将各向异性的复杂度隔离在极化函数模块，保持传播子模块的简洁性

3. **开发单元测试** [test/test_meson_propagator.jl]
   - 基本功能测试：验证`meson_propagator_simple`和`meson_propagator_mixed`返回ComplexF64、单位正确（fm²）
   - **通道-K系数映射测试**：验证π和K使用`K123_plus/K4567_plus`，σ_π和σ_K使用`K123_minus/K4567_minus`，η/η'使用K⁺系列，σ/σ'使用K⁻系列
   - 物理约束测试：检查k→0极限下传播子实部发散（质量壳条件）、虚部为零（无衰变）
   - 混合介子测试：验证η/η'混合时矩阵乘法计算的正确性，与手动展开结果对比（相对误差<1e-10）
   - **散射道映射测试**：
     - **t道测试**：u+d→u+d（纯夸克散射，直接传递）
     - **s道测试**：u+d̅→u+d̅（夸克-反夸克湮灭，验证反夸克的正确处理）
     - **u道测试**：u+d→d+u（交叉散射）
     - **场算符验证**：检查不同散射道中ψ和ψ̄的产生/湮灭性质是否符合物理预期
   - 各向异性修正测试：分别用ξ=0和ξ=0.3计算极化函数Π，再传入传播子函数对比结果变化
   - 夸克味配置测试：测试u̅u、u̅s、s̅s不同夸克组合的传播子计算（涵盖所有波函数常量）
   - 参考`test/test_effective_couplings.jl`结构：使用@testset组织测试、预生成积分节点、输出测试总结

4. **创建测试总结文档** [test/test_meson_propagator_summary.md]
   - 记录测试通过率：xx/xx测试全部通过
   - 数值精度验证：矩阵乘法展开对比的相对误差结果
   - 物理合理性检验：传播子量级、各向异性修正幅度、不同介子类型的差异
   - 性能基准：预计算优化带来的速度提升（预期~3-5倍）
   - 已知限制：衰变宽度计算、Mott相变检测等高级功能留待后续步骤

### Further Considerations

1. **极化函数预计算策略** - ✅已确认：本模块接收ComplexF64类型的预计算Π值。网格插值实现推迟至步骤6性能优化阶段。各向异性通过`PolarizationAniso.jl`模块处理，传播子函数不直接依赖ξ参数。
   
2. **Gell-Mann矩阵的完整性** - ✅已确认：当前只需实现λ₀和λ₈。其他矩阵（λ₃到λ₇）在公式文档中标记为"供参考，本模块暂不需要"。按需添加，避免过早优化。

3. **介子质量提取功能** - ✅已确认：本步骤专注于传播子计算。质量提取功能留待后续步骤或作为独立工具函数实现。

4. **散射过程的粒子类型约束** - 当前实现支持任意夸克/反夸克组合，包括：
   - 纯夸克散射：q₁q₂→q₃q₄（t道、u道）
   - 夸克-反夸克散射：qq̄→qq̄（s道为主，t道也支持）
   - 未来可能需要的扩展：反夸克-反夸克散射（需验证场算符映射）

5. **K系数命名约定** - 代码中使用完整变量名（`K123_plus`），物理讨论中可简写为K₃⁺（代表K₁=K₂=K₃的简并值）。文档中两种写法应保持一致对应关系。

### 关键公式修正总结（已反映在上述步骤中）

根据公式文档修正，以下是实现时必须注意的关键点：

1. **通道-K系数映射**（与直觉相反）：
   - ❌ 错误：赝标量P用K⁻，标量S用K⁺
   - ✅ 正确：赝标量P用K⁺（`*_plus`），标量S用K⁻（`*_minus`）

2. **K系数命名约定**：
   - π介子：使用`K123_plus`（代表K₁⁺=K₂⁺=K₃⁺的简并值）
   - K介子：使用`K4567_plus`（代表K₄⁺=K₅⁺=K₆⁺=K₇⁺的简并值）
   - 混合介子：使用`K0_plus/minus`, `K8_plus/minus`, `K08_plus/minus`
   - `meson_propagator_simple`函数已移除冗余的`channel`参数，介子类型隐含通道信息

3. **K介子极化函数**：
   - ❌ 错误：使用Π_{uu}
   - ✅ 正确：使用Π_{us}（因为K介子是u̅s或d̅s组合）

4. **M₀₈系数**：
   - ❌ 错误：4/(3√2) ≈ 0.9428
   - ✅ 正确：`4*sqrt(2)/3` ≈ 1.8856

5. **数据类型**：
   - 所有传播子和M矩阵使用ComplexF64
   - 极化函数Π作为ComplexF64传入
   - K系数和det_K为Float64实数

6. **M矩阵对称性**：
   - M₀₈ = M₈₀（可利用此减少计算量）

7. **旋量与味空间处理**：
   - Dirac旋量指标通过矩阵乘法$\bar{\psi}\lambda\psi$自动收缩
   - 本模块只处理味空间的3×3矩阵运算（Gell-Mann矩阵和夸克波函数）
   - 矩阵乘法顺序：从左到右依次进行（1×3矩阵 × 3×3矩阵 × 3×1向量）
