## Plan: 实现散射矩阵元模块ScatteringAmplitude.jl

在完成步骤3(介子传播子模块)的基础上，实现从总传播子到散射矩阵元平方的计算，支持所有11个散射过程(4个qq散射+7个qqbar散射)，直接计算|M|²而无需单独构建振幅M。

### Steps

1. **扩展TotalPropagator模块添加通道分离接口** - 在`TotalPropagator.jl`中实现`calculate_all_propagators_by_channel`函数，返回标量(S)和赝标量(P)通道分离的传播子：qq散射返回`(t_S, t_P, u_S, u_P)`，qqbar散射返回`(t_S, t_P, s_S, s_P)`。该函数内部需要调用修改后的`total_propagator_simple`和`total_propagator_mixed`，分别传入`:S`或`:P`通道参数。

2. **实现质心系动量计算辅助函数** - 创建`calculate_cms_momentum`函数，根据Mandelstam变量(s, t, u)和散射道类型(`:t`/`:u`/`:s`)计算传播子所需的介子四动量(k0, k)。质心系逻辑：先计算四个粒子能量E1=(s+m1²-m2²)/(2√s), E2=(s-m1²+m2²)/(2√s), E3=(s+m3²-m4²)/(2√s), E4=(s-m3²+m4²)/(2√s)；s道：k0=|E1+E2|, k=0；t道：k0=|E1-E3|, k=√(k0²-t)（需检查k0²-t≥0）；u道：k0=|E1-E4|, k=√(k0²-u)。

3. **创建ScatteringAmplitude.jl核心模块** - 在`src/relaxtime/`下创建新文件，实现`scattering_amplitude_squared`函数，根据公式第9节直接计算|M|²(包含1/(4N_c²)色自旋平均因子)。实现辅助函数`calculate_mandelstam_variables`(计算s_ij^±, t_ij^±, u_ij^±共18个辅助变量)和`get_quark_masses`(根据process从散射过程解析四个粒子质量，使用`quark_params`参数获取)。

4. **实现qq散射矩阵元平方** - 按公式9.1节实现夸克-夸克散射，计算三项：|M_u|²、|M_t|²、交叉项2Re(M_u M_t*)。调用`calculate_cms_momentum`分别计算t道和u道的(k0_t, k_t)和(k0_u, k_u)，然后调用`calculate_all_propagators_by_channel`获取4个传播子(D_t^S, D_t^P, D_u^S, D_u^P)。自动计算u变量(u = Σm² - s - t)。支持4个qq过程。

5. **实现qqbar散射矩阵元平方** - 按公式9.2节实现夸克-反夸克散射，计算三项：|M_s|²、|M_t|²、交叉项2Re(M_s M_t*)。调用`calculate_cms_momentum`分别计算s道和t道的(k0_s, k_s)和(k0_t, k_t)，获取4个传播子(D_s^S, D_s^P, D_t^S, D_t^P)。支持7个qqbar过程，正确处理不同过程的质量配置(如uū→ss中m3=m4=m_s)。

6. **创建API文档和单元测试** - 编写`api/ScatteringAmplitude.md`(包含参数表、Mandelstam变量单位fm^-2、返回值单位fm^-2、11个过程列表、质心系计算说明、使用示例)，编写`test_unit/test_scattering_amplitude.jl`(验证物理约束如|M|²≥0、质心系动量计算正确性、对称性检查、温度扫描、与零各向异性ξ=0对比)。

### Further Considerations

1. **质心系边界条件处理** - 当k0²-t或k0²-u为负但接近零(-1e-12 < Δ < 0)时，设k=0以处理数值误差。是否需要对更小的负值抛出错误或返回NaN？建议在测试中验证Mandelstam约束s+t+u=Σm²是否严格满足。

2. **TotalPropagator接口修改影响** - 添加通道分离接口后，原`calculate_all_propagators`函数保持不变(返回混合S+P传播子)以兼容现有代码。新接口命名为`calculate_all_propagators_by_channel`，需在文档中明确两者的使用场景差异。

3. **批量计算优化** - 可实现`scattering_amplitude_squared_batch(process, s_array, t_array, ...)`用于扫描参数空间。内部可复用K系数和A函数值，预计性能提升2-3倍。在步骤6完成后评估是否需要此功能。
