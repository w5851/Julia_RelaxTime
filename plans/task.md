# 弛豫时间计算模块开发路线图

## 项目概述

基于PNJL模型的弛豫时间计算项目，从当前A_correction/A_aniso函数的完善出发，逐步构建有效耦合常数→介子传播子→散射矩阵元→弛豫时间的完整计算链。

## 开发步骤

### ✅ 步骤 1: 完善OneLoopIntegralsAniso模块的文档与测试 【已完成 2025-11-14】

**目标**: 为新添加的A_correction和A_aniso函数补充完整的API文档和单元测试

**完成情况**:
- ✅ **API文档** (`api/OneLoopIntegralsAniso.md`, 555行)
  - 详细说明了`A_correction`返回一阶修正项，需与`A`相加使用
  - 详细说明了`A_aniso`直接返回完整各向异性结果
  - 明确了两者的区别、适用场景（ξ<0.3 vs ξ≥0.3）和性能对比（1.5-2倍耗时差异）
  
- ✅ **单元测试** (`test/test_oneloopintegrals_aniso.jl`, 455行，74/74测试通过)
  - ξ=0边界测试：验证`A_correction(ξ=0) ≈ 0`和`A_aniso(ξ=0) ≈ A` ✓
  - 数值对比测试：扫描ξ=0.01~0.5，相对误差从0.01%到1.09% ✓
  - 线性近似验证：ratio检验通过，偏差<0.3% ✓
  - 系统性测试：完成T、μ、m参数扫描 ✓
  - 收敛性测试：确定最优节点配置（64×32） ✓
  - 性能基准测试：详细对比两种方法 ✓
  
- ✅ **测试总结文档** (`test/test_oneloopintegrals_aniso_summary.md`, 490行)
  - 相对误差分析：ξ≤0.2时误差<0.3%，ξ=0.5时误差1.09%
  - 适用范围建议：一阶近似在ξ<0.3时可靠
  - 收敛性结果：64动量节点和64×32双重积分配置
  - 性能对比：包含批量计算优化建议
  - 物理合理性检验：所有物理约束验证通过

**关键修复**:
- 🔧 修复了`QuarkDistribution_Aniso.jl`中df/dE符号错误（line 35, 63）
  - 问题：`df_dE = -β_fm * df_dx`（多余负号）
  - 修复：`df_dE = β_fm * df_dx`（正确链式法则）
  - 影响：修复后delta_f与correction符号一致，相对误差从~100%降至<0.3%

**测试结果摘要**:
| ξ值 | 相对误差 | 适用性 |
|-----|---------|--------|
| 0.10 | 0.06% | ⭐⭐⭐⭐⭐ 优秀 |
| 0.20 | 0.21% | ⭐⭐⭐⭐⭐ 优秀 |
| 0.30 | 0.43% | ⭐⭐⭐⭐ 良好 |
| 0.50 | 1.09% | ⭐⭐ 勉强可用 |

**依赖关系**:
- 输入依赖：`OneLoopIntegrals.jl`中的A函数、`QuarkDistribution_Aniso.jl`中的分布函数
- 文档参考：`api/OneLoopIntegrals.md`的格式、`test_b0_correction.jl`的测试结构

---

### ✅ 步骤 2: 构建有效耦合常数模块EffectiveCouplings.jl 【已完成 2025-11-14】

**目标**: 实现从单圈积分A到有效耦合常数K的计算

**完成情况**:
- ✅ **源代码模块** (`src/relaxtime/EffectiveCouplings.jl`, 220行)
  - 实现了`calculate_G_from_A`转换函数
  - 实现了`calculate_effective_couplings`主函数
  - 实现了`coupling_matrix_determinant`辅助函数
  - 完整的文档字符串和使用示例

- ✅ **API文档** (`api/EffectiveCouplings.md`, 详细记录)
  - 详细的参数说明表（包含典型值和取值范围）
  - 物理意义解释（味道通道、手征极限、SU(3)对称性）
  - 温度和化学势依赖性分析
  - 完整的使用示例和注意事项

- ✅ **单元测试** (`test/test_effective_couplings.jl`, 387行，96/96测试通过)
  - 基本功能测试：验证三个函数的返回值和数值正确性 ✓
  - 手征极限测试：K=0和G^f=0时K_α退化为G ✓
  - SU(3)对称性测试：G^u=G^s时味道简并和K₀₈=0 ✓
  - 物理约束测试：det K > 0（标量和赝标量通道） ✓
  - K系数量级检验：所有系数在合理范围内 ✓
  - 温度扫描测试：T=50-250 MeV，检验手征相变效应 ✓
  - 化学势扫描测试：μ=0-300 MeV，检验η-η'混合增强 ✓
  - 完整计算流程演示：从A函数到K系数的端到端计算 ✓

**测试结果摘要**:

| 测试类型 | 测试数量 | 通过率 | 关键结果 |
|---------|---------|-------|---------|
| 基本功能 | 31 | 100% | 所有函数正确返回预期值 |
| 手征极限 | 18 | 100% | K=0和G^f=0时严格退化为G |
| SU(3)对称性 | 6 | 100% | 味道简并关系验证通过 |
| 物理约束 | 12 | 100% | det K > 0，量级合理 |
| 数值稳定性 | 27 | 100% | 温度和化学势扫描无异常 |
| 完整流程 | 2 | 100% | 端到端计算成功 |

**典型物理结果**（T=150 MeV, μ=0）:

| 系数 | 数值 (fm²) | 物理意义 |
|------|-----------|---------|
| K₁₂₃⁻ | 1.776×10⁻¹ | π介子通道 |
| K₄₅₆₇⁻ | 1.731×10⁻¹ | K介子通道 |
| K₀₈⁻ | -2.125×10⁻³ | η-η'混合 |
| det K^P | 4.147×10⁻² fm⁴ | 赝标量行列式 |

**依赖关系验证**:
- 输入依赖：`OneLoopIntegrals.A`函数正常工作 ✓
- 输出接口：返回NamedTuple供下游模块使用 ✓
- 公式一致性：与`doc/formula/K_有效耦合常数byA.md`完全一致 ✓

**任务清单**（原计划）:
- ~~在`src/relaxtime/`下创建`EffectiveCouplings.jl`模块~~ ✓
  


---

### 步骤 3: 构建介子传播子模块MesonPropagator.jl

**目标**: 实现从极化函数Π和耦合常数K到介子传播子D的计算

**任务清单**:
- 在`src/relaxtime/`下创建`MesonPropagator.jl`模块
  
- 实现主函数`meson_propagator(meson_type, p0, p, m_quarks, μ_quarks, T, Φ, Φbar, ξ, K_coeffs)`：
  - 输入参数：
    - `meson_type`: Symbol类型，支持`:pi`, `:K`, `:eta`, `:eta_prime`, `:sigma`, `:sigma_prime`
    - `p0`: 能量分量（MeV或fm⁻¹）
    - `p`: 三维动量模长
    - `m_quarks`: 夸克质量数组`[m_u, m_d, m_s]`
    - `μ_quarks`: 化学势数组`[μ_u, μ_d, μ_s]`
    - `K_coeffs`: 来自`EffectiveCouplings.jl`的耦合系数
  - 返回：复数类型的传播子`D = real + im*imag`（单位fm⁻²）
  
- 实现不同介子类型的计算逻辑：
  - **一般介子**（π、K、σ等）：`D^{P(S)} = 2K_α^+ / (1 - 4K_α^+ Π^{P(S)})`
  - **混合介子**（η/η'、σ/σ'）：使用2×2矩阵求逆
    - `D = 2 det(K) / (M₀₀M₈₈ - M₀₈²)`
    - 其中`M_ij = δ_ij - 4K_ij Π_j`
  
- 调用极化函数：
  - 从`PolarizationAniso.polarization_aniso`获取Π值
  - 根据介子类型选择`:P`（赝标量）或`:S`（标量）通道
  - 预计算并复用A值以提升性能
  
- 实现辅助功能：
  - `extract_meson_mass(meson_type, ...)`: 通过求解`Re[D⁻¹(m, 0)] = 0`提取物理质量
  - `check_mott_transition(m_meson, m_q1, m_q2)`: 检测Mott相变（m_meson ≈ m_q1 + m_q2）
  
- 在`api/MesonPropagator.md`中记录：
  - 各介子类型对应的夸克味道组合
  - 复数返回值的物理意义（实部对应质量壳，虚部对应衰变宽度）
  - Mott相变的物理解释和数值处理
  - 典型传播子值和单位换算
  
- 在`test/test_meson_propagator.jl`中测试：
  - 基本功能：每种介子类型在标准参数下的传播子计算
  - 介子质量提取：与实验值对比（如m_π ≈ 140 MeV）
  - 简化情况：p0=0, p=0时的静态极限
  - 混合介子：η-η'混合角的正确性
  - 各向异性影响：ξ=0和ξ≠0的对比

**依赖关系**:
- 输入依赖：`EffectiveCouplings.jl`、`PolarizationAniso.jl`
- 输出用途：提供给`ScatteringAmplitude.jl`使用
- 公式来源：`doc/formula/Propagator_传播子byPolarization.md`

---

### 步骤 4: 从PDF论文中提取散射矩阵元公式

**目标**: 规范化记录散射矩阵元的理论公式，为代码实现做准备

**任务清单**:
- 手动从外部PDF论文中读取散射矩阵元公式（涉及介子交换的t通道贡献）
  
- 在`doc/formula/`下创建`ScatteringAmplitude_散射矩阵元.md`，按照`公式文档模板.md`格式：
  
  1. **公式标识**：
     - 程序实现：`scattering_amplitude`
     - 对应文件：`src/relaxtime/ScatteringAmplitude.jl`
  
  2. **物理意义**：
     - 物理背景：PNJL模型中夸克-夸克弹性散射q₁q₂→q₃q₄，通过介子传播子描述
     - 在计算链中的作用：连接介子传播子和散射截面
     - 相关物理量：Mandelstam变量s、t、u，介子传播子D，味因子
  
  3. **数学表达式**：
     - Mandelstam变量定义：
       - s = (p₁ + p₂)²
       - t = (p₁ - p₃)²
       - u = (p₁ - p₄)²
       - 约束：s + t + u = Σm_i²
     - 总振幅：`M_total = M_π + M_K + M_η + M_η' + ...`（不同介子通道求和）
     - 单通道形式：`M_α = g_α² D_α(t)`（其中g_α是顶点因子）
  
  4. **参数说明表**：
     - 输入：散射过程标识（如`:uu_to_uu`）、动力学变量s/t、热力学参数T/μ
     - 输出：复数振幅M（单位：无量纲或fm）
     - 味因子：不同夸克组合的耦合系数
  
  5. **依赖关系**：
     - 输入依赖：`MesonPropagator.meson_propagator`
     - 输出用途：计算微分散射截面`dσ/dt ∝ |M|²`
  
  6. **各向异性处理**：
     - 明确角度依赖：t = -2p²(1 - cosθ)，其中θ是散射角
     - 各向异性情况下，需对角度积分时考虑ξ的影响

**注意事项**:
- 记录公式来源：论文作者、年份、期刊、章节号、页码
- 标注近似条件：是否采用软动量交换近似、是否忽略u通道贡献等
- 明确同位旋因子和色因子的处理

---

### 步骤 5: 实现散射矩阵元模块ScatteringAmplitude.jl

**目标**: 基于步骤4的公式文档，实现散射振幅计算

**任务清单**:
- 在`src/relaxtime/`下创建`ScatteringAmplitude.jl`模块
  
- 实现主函数`scattering_amplitude(process, s, t, m_quarks, μ_quarks, T, Φ, Φbar, ξ)`：
  - 输入：
    - `process`: Symbol类型，如`:uu_to_uu`, `:ud_to_ud`, `:us_to_us`等
    - `s, t`: Mandelstam变量（fm⁻²）
    - 其他热力学参数
  - 输出：复数振幅M
  
- 实现通道分解逻辑：
  - 根据`process`确定参与的介子通道（如u-u散射包含π、η、η'交换）
  - 对每个通道调用`meson_propagator`获取D_α(t)
  - 应用味因子和顶点因子g_α
  - 求和得到总振幅：`M = Σ_α g_α² D_α(t)`
  
- 处理对称性和交叉对称性：
  - t通道和u通道的关系
  - 全同粒子的Bose/Fermi统计因子
  
- 在`api/ScatteringAmplitude.md`中记录：
  - 所有支持的散射过程列表
  - 单位约定（自然单位制，能量和动量单位为fm⁻¹）
  - 振幅的物理意义和典型量级
  - 使用示例：从入射动量计算振幅的完整代码
  
- 在`test/test_scattering_amplitude.jl`中测试：
  - 前向散射极限：t→0时振幅的行为
  - 后向散射：t→t_min的极限
  - 对称性检验：M(q₁q₂→q₃q₄) 与 M(q₁q₃→q₂q₄) 的关系
  - 不同温度下的振幅变化
  - 各向异性修正的量级（ξ=0 vs ξ=0.3）

**依赖关系**:
- 输入依赖：`MesonPropagator.jl`、`EffectiveCouplings.jl`
- 输出用途：提供给`CrossSection.jl`使用
- 公式来源：`doc/formula/ScatteringAmplitude_散射矩阵元.md`

---

### 步骤 6: 完成剩余模块（散射截面、驰豫时间）的构建循环

**目标**: 重复"公式提取→实现→测试"的流程，完成计算链的最后两步

#### 6.1 散射截面模块

**任务清单**:
- 从PDF论文提取公式，创建`doc/formula/CrossSection_散射截面.md`：
  - 微分截面：`dσ/dt = (1/16πs) |M(s,t)|²`
  - 总截面：`σ_total = ∫ dσ/dt dt`（积分范围：t_min到t_max）
  - Mandelstam变量的运动学边界
  
- 实现`src/relaxtime/CrossSection.jl`模块：
  - `differential_cross_section(process, s, t, ...)`: 返回dσ/dt（单位fm²）
  - `total_cross_section(process, s, ...)`: 对t积分得到σ_total
  - 使用自适应积分（QuadGK）处理可能的奇点
  
- 在`api/CrossSection.md`中记录：
  - 截面的物理意义（单位面积的散射概率）
  - 单位换算：fm² ↔ mb（毫靶恩）
  - 典型截面值和温度依赖性
  
- 在`test/test_cross_section.jl`中测试：
  - 单位正确性：验证dσ/dt积分后确实得到正确的σ_total
  - 低能/高能极限：与已知理论结果对比
  - 光学定理：Im[M(s, t=0)] 与 σ_total 的关系
  - 各向异性影响：ξ对截面的修正量级

#### 6.2 驰豫时间模块

**任务清单**:
- 从PDF论文提取公式，创建`doc/formula/RelaxationTime_驰豫时间.md`：
  - 平均散射率：`ω = ∫ dΓ σ(s) v_rel f₂(1±f₃)(1±f₄)`
  - 驰豫时间：`τ = 1/ω`
  - 相空间积分的处理（考虑Pauli阻塞/Bose增强）
  
- 实现`src/relaxtime/RelaxationTime.jl`模块：
  - `scattering_rate(p, m, μ, T, Φ, Φbar, ξ)`: 计算给定动量p的夸克散射率ω
  - `relaxation_time(p, m, μ, T, Φ, Φbar, ξ)`: 返回τ = 1/ω（单位fm）
  - 实现多重积分（对靶粒子动量p₂和出射粒子动量p₃、p₄积分）
  - 考虑各向异性分布函数的角度依赖
  
- 在`api/RelaxationTime.md`中记录：
  - 驰豫时间的物理意义（系统趋于平衡的特征时间）
  - 温度依赖性：高温下τ↓（散射更频繁）
  - 各向异性修正的效应：ξ如何影响τ
  - 与输运系数的关系：剪切粘滞η ∝ pτ
  
- 在`test/test_relaxation_time.jl`中测试：
  - 温度扫描：验证τ(T)的单调性和物理合理性
  - 动量依赖：τ(p)在低动量和高动量极限的行为
  - 化学势依赖：μ对τ的影响
  - 各向异性修正：扫描ξ=0~1，找到一阶近似失效的临界值
  - 与文献对比：在相同参数下与已发表结果比较

**依赖关系**:
- 输入依赖：`CrossSection.jl`、`QuarkDistribution_Aniso.jl`
- 输出用途：计算输运系数（剪切粘滞系数η/s、体粘滞系数ζ/s）
- 公式来源：`doc/formula/平均散射率_动量各向异性.md`（已存在）和新创建的文档

---

## 进一步考虑

### 1. 公式来源的规范化

**建议**:
- 在`doc/formula/`下创建`References.md`，统一记录所有公式的来源：
  - 论文标题、作者、期刊、年份、DOI
  - 每个公式文档对应的章节号和页码
  - 如有修改或简化，注明理由
  
**好处**:
- 便于公式验证和溯源
- 方便后续研究者理解理论基础
- 确保公式的准确性和一致性

### 2. 性能优化点

**潜在瓶颈**:
- 介子传播子计算需要多次调用极化函数Π，而Π又依赖于B0的复杂积分
- 散射截面和驰豫时间涉及多重积分（3维到9维）
  
**优化策略**:
- **预计算与插值**：在固定T、μ下，预计算Π(k0, k)在网格上的值，构建插值表
- **并行化**：利用Julia的多线程特性（`@threads`）并行计算不同动量点
- **自适应积分**：对平滑区域使用低精度（rtol=1e-3），奇点附近提高精度
- **缓存机制**：对于重复调用的函数（如分布函数），使用LRU缓存
  
**实施计划**:
- 在步骤3（MesonPropagator）实现Π的网格插值
- 在步骤6.2（RelaxationTime）实现多线程并行
- 在每个模块的测试中添加性能基准测试

### 3. 各向异性修正的完整性

**当前状态**:
- A_correction和B0_correction是ξ的一阶近似
- A_aniso是完整的Romatschke-Strickland分布

**待确认问题**:
1. 散射矩阵元是否需要考虑ξ的高阶修正？
2. 一阶近似在何种参数范围内有效（ξ < ?）？
3. 各向异性对介子质量的修正是否显著？

**验证计划**:
- 在最终的`test_relaxation_time.jl`中系统扫描ξ=0~1：
  - 对比一阶近似（A + A_correction）与完整计算（A_aniso）的相对误差
  - 绘制误差随ξ变化的曲线，确定一阶近似的适用范围
  - 检查ξ>0.5时是否需要引入二阶修正项
  
- 在各模块的文档中明确说明各向异性修正的阶数和适用条件

### 4. 代码质量保证

**遵循规范**:
- **命名规范**（`prompt/变量命名规范.md`）：严格使用`_inv_fm`后缀表示fm⁻¹单位
- **代码风格**（`prompt/代码风格指南.md`）：性能关键代码使用`@fastmath`、`@inline`、`@simd`
- **测试覆盖**：每个公开函数都有对应的单元测试

**文档一致性**:
- 所有API文档包含"参数说明表"和"使用示例"
- 所有测试总结文档包含"相对误差分析"和"性能对比"
- 公式文档严格按照`公式文档模板.md`的7个章节结构

**版本控制**:
- 每完成一个模块，提交一次git commit，附上详细的commit message
- 重要修复（如之前的B0_correction虚部修复）需记录在测试总结文档中

---

## 时间估算与里程碑

| 步骤 | 模块名称 | 预计工作量 | 里程碑 |
|------|---------|-----------|-------|
| 1 | OneLoopIntegralsAniso文档+测试 | 2-3天 | 完成A函数的各向异性计算链 |
| 2 | EffectiveCouplings.jl | 3-4天 | 实现A→K的转换 |
| 3 | MesonPropagator.jl | 5-7天 | 完成介子物理的核心计算 |
| 4 | ScatteringAmplitude公式提取 | 1-2天 | 理论基础准备完毕 |
| 5 | ScatteringAmplitude.jl | 4-5天 | 实现散射振幅计算 |
| 6.1 | CrossSection.jl | 3-4天 | 完成散射截面计算 |
| 6.2 | RelaxationTime.jl | 5-7天 | 完成弛豫时间主函数 |
| 优化与验证 | 性能优化+系统测试 | 7-10天 | 项目整体完成 |

**总计**: 约30-42天（6-8周）

**关键里程碑**:
- 第2周末：步骤1-2完成，验证K系数计算正确性
- 第4周末：步骤3完成，能够计算介子传播子和质量
- 第6周末：步骤4-5完成，实现散射矩阵元
- 第8周末：步骤6完成，输出第一版驰豫时间结果

---

## 成功标准

**技术指标**:
1. 所有模块通过单元测试，测试覆盖率>90%
2. 各向异性修正在ξ<0.3时与一阶近似误差<5%
3. 计算的介子质量与实验值误差<10%（m_π、m_K）
4. 驰豫时间的温度依赖性符合物理直觉（T↑则τ↓）

**文档完整性**:
1. 每个模块有对应的API文档、公式文档、测试总结
2. 所有公式可追溯到原始论文
3. 使用示例代码可直接运行

**性能要求**:
1. 单次驰豫时间计算<10秒（在标准参数下）
2. 温度扫描（10个点）<2分钟
3. 完整相图（100×100网格）<4小时

**物理合理性**:
1. 高温极限：τ → 0（理想气体极限）
2. 低温极限：τ → ∞（强耦合，无散射）
3. 各向异性效应：ξ↑导致某些方向τ↓（散射增强）

---

## 附录：现有资源总结

**已完成模块**（✅）:
- `Constants_PNJL.jl` - 物理常数和PNJL参数
- `QuarkDistribution.jl` - 各向同性夸克分布函数
- `QuarkDistribution_Aniso.jl` - 各向异性分布函数（RS分布）
- `OneLoopIntegrals.jl` - A和B0函数（各向同性）
- `OneLoopIntegralsAniso.jl` - A_correction、B0_correction、A_aniso
- `PolarizationAniso.jl` - 极化函数（支持各向异性）
- `GaussLegendre.jl` - 高斯-勒让德积分工具
- `CauchyPV.jl` - 柯西主值积分工具

**现有公式文档**（15个）:
- 基础积分：A.md、B0.md、B0_extra.md、B0_各向异性一般情况.md等
- 分布函数：PNJL_夸克有效分布函数.md、PNJL_夸克有效分布函数_动量各向异性.md
- 高层计算：K_有效耦合常数byA.md、Polarization_极化函数byB0.md、Propagator_传播子byPolarization.md
- 散射率：平均散射率_动量各向异性.md

**现有测试**（8个）:
- test_gausslegendre.jl、test_cauchypv.jl
- test_quark_distribution_antiderivative.jl、test_quark_distribution_aniso.jl
- test_oneloopintegrals.jl、test_b0_correction.jl、test_polarization_aniso.jl
- test_singularity_inequality.jl

**现有API文档**（4个）:
- OneLoopIntegrals.md、PolarizationAniso.md
- GaussLegendre.md、CauchyPV.md

**待创建资源**（下一阶段）:
- API文档：OneLoopIntegralsAniso.md、EffectiveCouplings.md、MesonPropagator.md、ScatteringAmplitude.md、CrossSection.md、RelaxationTime.md
- 公式文档：ScatteringAmplitude_散射矩阵元.md、CrossSection_散射截面.md、RelaxationTime_驰豫时间.md
- 源代码：EffectiveCouplings.jl、MesonPropagator.jl、ScatteringAmplitude.jl、CrossSection.jl、RelaxationTime.jl
- 测试代码：test_oneloopintegrals_aniso.jl、test_effective_couplings.jl、test_meson_propagator.jl、test_scattering_amplitude.jl、test_cross_section.jl、test_relaxation_time.jl

---

**计划创建日期**: 2025年11月14日  
**目标完成日期**: 2025年12月底至2026年1月初  
**项目状态**: 
- ✅ **第一阶段（步骤1）已完成** (2025-11-14): OneLoopIntegralsAniso模块文档与测试完备，df/dE符号错误已修复
- ✅ **第二阶段（步骤2）已完成** (2025-11-14): EffectiveCouplings模块实现、文档与测试完备，96个测试全部通过
- 📋 **准备开始第三阶段（步骤3）**: 构建介子传播子模块MesonPropagator.jl