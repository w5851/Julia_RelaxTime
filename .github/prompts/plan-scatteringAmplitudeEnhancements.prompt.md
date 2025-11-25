## Plan: 完善散射矩阵元模块

为散射矩阵元模块添加新散射过程、批量计算功能、性能测试，并修复极化函数ξ<0的问题。

### Steps

1. **诊断并修复ξ<0极化函数问题**
   - 在`src/relaxtime/OneLoopIntegralsAniso.jl`中检查`B0_correction`函数对ξ符号的处理逻辑
   - 创建`test_unit/test_polarization_xi_negative.jl`测试ξ正负值的差异
   - 检查是否错误使用`abs(ξ)`或条件判断跳过ξ<0的情况
   - 修复bug或在文档中说明物理原因

2. **添加缺失的电荷共轭散射过程**
   - 在`src/relaxtime/ScatteringAmplitude.jl`的`SCATTERING_MESON_MAP`中添加`:dubar_to_dubar` (d+ū→d+ū) 和`:subar_to_subar` (s+ū→s+ū)
   - `:dubar_to_dubar`介子映射：t道为`[:pi, :sigma_pi, :mixed_P => true, :mixed_S => true]`，s道为`[:pi, :sigma_pi]`
   - `:subar_to_subar`介子映射：t道为`[:mixed_P => true, :mixed_S => true]`（仅η/η'/σ/σ'），s道为`[:K, :sigma_K]`
   - 验证新过程与`:udbar_to_udbar`和`:usbar_to_usbar`在同位旋对称下物理等价（可通过数值对比散射矩阵元）

3. **实现批量计算函数**
   - 在`src/relaxtime/ScatteringAmplitude.jl`中添加`calculate_all_scattering_amplitudes_squared(s, t, quark_params, thermo_params, K_coeffs)`
   - 返回NamedTuple包含所有13种散射过程的|M|²值（包括新增的2种）
   - 添加到模块导出列表并编写API文档

4. **实现避免缓存的性能测试**
   - 在`test_unit/test_scattering_amplitude.jl`中添加批量计算性能测试
   - 每次迭代微调s和t参数（变化量>1e-11，超过缓存容差1e-12）
   - 测试1000次迭代，记录总耗时、平均用时和单过程平均用时
   - 验证缓存避免策略的有效性（对比重置缓存vs参数微调的性能差异）

### Further Considerations

1. **新散射过程的物理意义**：`:dubar_to_dubar`和`:subar_to_subar`与`:udbar_to_udbar`、`:usbar_to_usbar`在同位旋对称下矩阵元相同，但在弛豫时间计算中因粒子数密度不同而贡献不同。是否需要在文档中明确说明这种对称性？
A:是的，在文档中说明对称性

2. **介子映射的关键差异**：dubar过程的t道包含π介子（轻夸克对uu/dd都能形成），而subar过程的t道不包含π（需要ss成分）。代码实现时需注意这一物理约束。

3. **性能测试的缓存策略**：参数微调法（s×(1+i×1e-11)）模拟真实物理场景，重置缓存法测量纯计算时间。建议同时实现两种策略并对比结果？
A:是的，同时实现

4. **ξ<0问题的优先级**：当前测试显示ξ∈[-1,0]时极化函数值不变，这可能影响各向异性修正的正确性。建议优先诊断此问题，因为它可能影响所有散射过程的计算？
A:好的，诊断此问题