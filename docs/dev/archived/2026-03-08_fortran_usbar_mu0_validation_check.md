---
title: Fortran 旧代码在 `usbar_to_usbar, mu_B=0, T=210/250 MeV` 上的运行与验证记录
archived: true
original: docs/dev/active/2026-03-08_fortran_usbar_mu0_validation_check.md
archived_date: 2026-03-08
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Fortran 旧代码在 `usbar_to_usbar, mu_B=0, T=210/250 MeV` 上的运行与验证记录

**记录日期**: 2026-03-08  
**目标**: 先确认旧 Fortran 项目是否能在当前工作区跑通，再检查其 `usbar_to_usbar` 在 `mu_B=0, T=210/250 MeV` 两个温度上是否和验证集对得上。  
**结论先行**: 旧项目可以构建并运行；但默认主程序不直接产出与验证集同口径的 `sigma(sqrt(s))` 曲线，需要补一个最小专用驱动。补驱动后得到的 Fortran 曲线与当前 Julia 诊断结论一致，即 `T=210 MeV` 明显峰位偏右且峰高过大，`T=250 MeV` 峰位较接近但峰高仍偏大，**并不能对上验证集**。

## 1. 旧 Fortran 项目可运行性检查

### 1.1 构建检查

- 使用 VS Code 任务 `Fortran: build all (generate .mod)` 成功编译全部 36 个源文件，`.mod` 正常生成到 `Relaxtime_fortran/mod/`。
- 该任务只负责模块编译，不会链接出新的可执行文件。
- 工作区中原本已有 `Relaxtime_fortran/codes/main/AA_MAIN.exe`，说明旧项目此前已能完成主程序链接。

### 1.2 主程序运行检查

- 从项目根目录直接执行 `codes/main/AA_MAIN.exe`，程序能够启动，不会立刻报错。
- 运行期间 `results/w_ij.dat` 成功更新，说明主程序至少能完成首个状态点的求解和输出。
- 但主程序默认流程是：
  - 固定 `mu_initial = mu_final = 5 MeV`
  - 扫描 `T = 50:350 MeV`
  - 输出松弛时间、平均散射率和输运系数
- 这条默认路径**不直接给出**验证集所需的 `usbar_to_usbar` 单过程 `sigma(sqrt(s))` 曲线，因此不能直接拿来和 digitized CSV 对比。

## 2. 旧代码中 `usbar_to_usbar` 的口径核对

### 2.1 过程码确认

- `Relaxtime_fortran/codes/mesons/module mesons.f90` 中定义了
  - `usb_code = reshape([1,3,1,3, 1,-1,1,-1, 1,-1,-1,1], [4,3])`
  - 注释为 `usb-usb`
- `Relaxtime_fortran/codes/relax time/z1 relax_time.f90` 中
  - `w(8)` 明确对应 `call averaged_rate( usb_code(:,:), w(8), wnn(8) ) ! 8.usbar->usbar`

### 2.2 直接验证时遇到的口径问题

- 旧主程序默认写出的是 `w_ij.dat` 中的平均散射率 `w(8)`，而不是 `cross_section.dat` 这种 `sigma(sqrt(s))` 曲线。
- `cross_section.f90` 虽然存在，但默认主程序并没有直接调用它来做 `usbar_to_usbar` 的专项扫描。
- 同时，`coefficients_tmu` 内含依赖输运系数计算的分支，对 `mu_B = 0` 的专项验证并不必要，且不适合拿来做最小化对表。

## 3. 为验证新增的最小专项驱动

### 3.1 新增文件

- 新增 `Relaxtime_fortran/codes/main/AA_USBAR_MU0_SCAN.f90`

### 3.2 该驱动做的事

- 只做两组定点:
  - `mu_B = 0, T = 210 MeV`
  - `mu_B = 0, T = 250 MeV`
- 对每个点执行:
  1. `set_integral`
  2. `newt(funcv_tmu)` 求平衡态 `phi_u, phi_d, phi_s, Phi, Phibar`
  3. `quantity` 计算当前态质量与热力学量
  4. `cross_section(usb_code)` 输出 `usbar_to_usbar` 的 `sigma(sqrt(s))`
- 该驱动绕开了旧主程序的整段输运扫描，只保留验证所需的最小链路。

### 3.3 编译方式

- 用 `gfortran -fopenmp` 对 `codes/` 下全部 Fortran 源文件重新链接。
- 链接时排除了旧的 `AA_MAIN.f90`，避免双 `program main` 冲突。
- 生成的专项可执行文件为:
  - `Relaxtime_fortran/codes/main/AA_USBAR_MU0_SCAN.exe`

## 4. 本次运行产物

专项驱动运行后生成了以下结果，并额外复制为语义化文件名以防被后续覆盖:

- `Relaxtime_fortran/results/cross_section.dat`
- `Relaxtime_fortran/results/usbar_mu0_state_summary.dat`
- `Relaxtime_fortran/results/usbar_mu0_cross_section_T210_T250.dat`
- `Relaxtime_fortran/results/usbar_mu0_state_summary_T210_T250.dat`

其中:

- `cross_section.dat` 共 600 行
  - 对应 `T=210 MeV` 300 点
  - 对应 `T=250 MeV` 300 点
- `usbar_mu0_state_summary.dat` 记录了两组点的平衡态和质量摘要

## 5. 定点状态摘要

从 `usbar_mu0_state_summary.dat` 读取到:

### 5.1 `T = 210 MeV, mu_B = 0`

- `check = F`，说明 `newt` 正常收敛
- `phi_u = phi_d = -3.3313153516832067e-1`
- `phi_s = -1.6526932509479366`
- `Phi = Phibar = 6.5283979473454001e-1`
- `m_u = m_d = 67.43029375495614 MeV`
- `m_s = 399.68954211996657 MeV`

### 5.2 `T = 250 MeV, mu_B = 0`

- `check = F`，说明 `newt` 正常收敛
- `phi_u = phi_d = -5.9535453233912031e-2`
- `phi_s = -1.0947335466632679`
- `Phi = Phibar = 7.6344587139797671e-1`
- `m_u = m_d = 15.956168869130746 MeV`
- `m_s = 310.96469966831575 MeV`

## 6. 与验证集的直接对比

### 6.1 对比数据源

- 验证集: `Julia_RelaxTime/tests/validation/data/relaxtime_sigma_literature_digitized_longtable_v1.csv`
- 过滤条件:
  - `process = usbar_to_usbar`
  - `muB_MeV = 0`
  - `T_MeV = 210` 或 `250`

### 6.2 对比方法

- 对 Fortran 输出 `cross_section.dat` 按温度拆分为两条曲线。
- 在验证集每个 `sqrt(s)` 位置上，对 Fortran 曲线做线性插值。
- 统计:
  - 峰位 `sqrt(s)_peak`
  - 峰高 `sigma_peak`
  - `MAD`（平均绝对误差）
  - `RMSE`
  - `MRD`（平均相对误差）

### 6.3 `T = 210 MeV`

- Fortran 峰位: `sqrt(s) = 0.496611 GeV`
- Fortran 峰高: `sigma = 51.202552 mb`
- 验证集峰位: `sqrt(s) = 0.467625 GeV`
- 验证集峰高: `sigma = 25.627966 mb`
- 误差摘要:
  - `MAD = 4.588342 mb`
  - `RMSE = 7.987449 mb`
  - `MRD = 0.688811`
  - `N = 62`

最差的几个点集中在文献峰附近:

- `T210MeV_usbartousbar_7`: `x = 0.497448 GeV`, validation `19.614585 mb`, Fortran `50.679422 mb`
- `T210MeV_usbartousbar_8`: `x = 0.509376 GeV`, validation `18.183520 mb`, Fortran `43.221740 mb`
- `T210MeV_usbartousbar_6`: `x = 0.486844 GeV`, validation `17.523395 mb`, Fortran `38.101093 mb`

判断:

- 峰位比验证集向右偏移约 `0.029 GeV`
- 峰高约为验证集的 `2.0` 倍
- 这和 Julia 侧现有诊断文档中的结论完全同向

### 6.4 `T = 250 MeV`

- Fortran 峰位: `sqrt(s) = 0.575951 GeV`
- Fortran 峰高: `sigma = 9.640291 mb`
- 验证集峰位: `sqrt(s) = 0.582275 GeV`
- 验证集峰高: `sigma = 5.131584 mb`
- 误差摘要:
  - `MAD = 1.137992 mb`
  - `RMSE = 1.684674 mb`
  - `MRD = 0.485093`
  - `N = 49`

最差的几个点集中在主峰附近:

- `T250MeV_usbartousbar_28`: `x = 0.582275 GeV`, validation `5.131584 mb`, Fortran `9.575839 mb`
- `T250MeV_usbartousbar_27`: `x = 0.553116 GeV`, validation `5.011060 mb`, Fortran `9.343859 mb`
- `T250MeV_usbartousbar_29`: `x = 0.611434 GeV`, validation `4.780014 mb`, Fortran `8.827315 mb`

判断:

- 峰位已经比较接近，左偏约 `0.006 GeV`
- 但峰高依然明显偏大，约为验证集的 `1.88` 倍
- 这同样和 Julia 侧既有诊断一致

## 7. 本轮结论

### 7.1 关于“旧 Fortran 项目能不能跑通”

- 可以。
- 构建层面没有阻塞，主程序也能实际启动并写出结果。

### 7.2 关于“旧 Fortran 代码能不能对上验证集”

- 在 `usbar_to_usbar, mu_B = 0, T = 210/250 MeV` 这两个关键温度上，**不能对上**。
- 旧 Fortran 结果和当前 Julia 诊断指向同一现象:
  - `T = 210 MeV`: 峰位偏右且峰高显著过大
  - `T = 250 MeV`: 峰位大致接近，但峰高显著过大

### 7.3 对后续排查的含义

- 这说明当前异常**不是 Julia 重写独有的问题**。
- 旧 Fortran 和 Julia 在这两个点上给出的偏差趋势一致，因此更像是:
  - 旧模型实现本身就带有该偏差，或
  - Julia 目前忠实复现了旧实现中的这一偏差来源。

## 8. 建议的下一步

- 以本次生成的 `usbar_mu0_cross_section_T210_T250.dat` 为基准，把 Fortran 与 Julia 在相同 `sqrt(s)` 采样点上做逐点直接对比。
- 若两者几乎重合，则后续工作重心应从“Julia 实现错误”转向“旧模型/旧公式口径为何与文献不符”。
- 若两者仍存在系统性差异，则可继续围绕:
  - 质量解
  - `B0` / 极化函数
  - strange `s` 道 `P` 通道
  - blocking / distribution 权重
  做 Julia-Fortran 的逐层对账。

## 9. Julia 与 Fortran 的逐点直接对比

### 9.1 对比输入

- Fortran 专项曲线:
  - `Relaxtime_fortran/results/usbar_mu0_cross_section_T210_T250.dat`
- Julia 现成结果:
  - `Julia_RelaxTime/data/outputs/results/relaxtime/cross_section/xs_T210_muB0_xi0p0.csv`
  - `Julia_RelaxTime/data/outputs/results/relaxtime/cross_section/xs_T250_muB0_xi0p0.csv`

补充说明:

- Fortran 输出列顺序为 `(sqrt_s_GeV, muB_MeV, T_MeV, sigma_mb)`。
- Julia CSV 为长表；本轮只筛选 `process = usbar_to_usbar`。
- 逐点对比时以 Fortran 的 `sqrt(s)` 采样点为基准，在 Julia 曲线上做线性插值。

### 9.2 生成的对比产物

- `Julia_RelaxTime/data/outputs/results/relaxtime/cross_section/fortran_compare/usbar_to_usbar_T210_julia_fortran_pointwise.csv`
- `Julia_RelaxTime/data/outputs/results/relaxtime/cross_section/fortran_compare/usbar_to_usbar_T250_julia_fortran_pointwise.csv`
- `Julia_RelaxTime/data/outputs/results/relaxtime/cross_section/fortran_compare/usbar_to_usbar_T210_T250_julia_fortran_summary.json`

点对点 CSV 额外记录了两套误差:

- 原始误差: `Julia - Fortran`
- `scaled10` 误差: `10 * Julia - Fortran`

之所以额外算 `scaled10`，是因为原始逐点结果一眼就显示 Julia 整体约为 Fortran 的 `0.1` 倍，需要直接验证这是否是稳定比例因子而不是局部偶然。

### 9.3 `T = 210 MeV`

- 重叠采样点数: `44`
- 原始逐点比例 `Julia / Fortran`:
  - mean = `0.100293`
  - median = `0.100140`
  - min = `0.095235`
  - max = `0.107767`
- 原始 RMSE: `12.855180 mb`
- 若把 Julia 整体乘以 `10` 后再比较，RMSE 降到 `0.503834 mb`
- 原始最大绝对误差: `46.326280 mb` at `sqrt(s) = 496.611 MeV`
- `scaled10` 最大绝对误差: `2.439828 mb` at `sqrt(s) = 496.611 MeV`
- 峰位对比:
  - Julia 峰位 = `499.222 MeV`
  - Fortran 峰位 = `496.611 MeV`
  - 峰位差 = `+2.611 MeV`
- 峰高对比:
  - Julia 峰高 = `4.931702 mb`
  - Fortran 峰高 = `51.202552 mb`
  - 原始峰高差 = `-46.270850 mb`
  - `scaled10` 峰高差 = `-1.885528 mb`

判断:

- 形状和峰位已经相当接近。
- 主要差异不是峰位置错了，而是 Julia 曲线整体几乎缩小了一个数量级。

### 9.4 `T = 250 MeV`

- 重叠采样点数: `49`
- 原始逐点比例 `Julia / Fortran`:
  - mean = `0.100132`
  - median = `0.100582`
  - min = `0.094579`
  - max = `0.102096`
- 原始 RMSE: `4.480155 mb`
- 若把 Julia 整体乘以 `10` 后再比较，RMSE 降到 `0.055521 mb`
- 原始最大绝对误差: `8.664144 mb` at `sqrt(s) = 575.951 MeV`
- `scaled10` 最大绝对误差: `0.176846 mb` at `sqrt(s) = 535.546 MeV`
- 峰位对比:
  - Julia 峰位 = `570.426 MeV`
  - Fortran 峰位 = `575.951 MeV`
  - 峰位差 = `-5.525 MeV`
- 峰高对比:
  - Julia 峰高 = `0.980711 mb`
  - Fortran 峰高 = `9.640291 mb`
  - 原始峰高差 = `-8.659580 mb`
  - `scaled10` 峰高差 = `+0.166817 mb`

判断:

- `T = 250 MeV` 的一致性比 `T = 210 MeV` 更高。
- 在乘以 `10` 之后，峰位和峰高都已经进入很小的残差范围。

### 9.5 逐点对账的结论

这一步已经能把问题进一步收缩到“整体归一化/单位因子”层面:

- Julia 与 Fortran 并不存在此前担心的那种“峰位完全跑开”或“曲线形状完全不同”的实现分叉。
- 两个温度下，`Julia / Fortran` 都稳定在约 `0.1`。
- 将 Julia 曲线整体乘以 `10` 后:
  - `T = 210 MeV` 的 RMSE 从 `12.855180 mb` 降到 `0.503834 mb`
  - `T = 250 MeV` 的 RMSE 从 `4.480155 mb` 降到 `0.055521 mb`

因此，当前最可疑的不是散射振幅的能量依赖形状，而是某个统一的截面归一化因子，例如:

- 单位换算链路
- mb/fm^2 因子
- 自旋/颜色平均或求和的整体因子
- 与旧 Fortran `cross_section.f90` 中 `cs*10d0` 对应的口径差异

### 9.6 对后续排查的直接建议

下一步优先级应调整为“找这个稳定 `x10` 差异从哪里进入”而不是继续做更细的峰形诊断。最直接的排查入口是:

1. 对照 Julia 与 Fortran 的最终截面输出公式，逐项检查是否存在额外或缺失的整体常数。
2. 回看 `cross_section.f90` 中 `cs*10d0` 的物理单位意图，并核对 Julia 侧最终输出单位是否已经提前吸收了同类因子。
3. 若总公式表面一致，再继续检查是否有统一作用于全部过程的颜色/自旋平均约定差异。

## 10. 单位链路核对结论

### 10.1 Fortran 的“内部单位”和“输出单位”不是一回事

从 `Relaxtime_fortran/codes/relax time/cross section.f90` 可直接看到，最终写盘语句是:

- `write(305, "(4(4XE17.10))") sqs*hc/1000d0, mu_B*hc, T*hc, cs*10d0`

这行代码至少说明了两件事:

- `sqrt(s)` 在输出前被乘 `hc/1000`，即从内部自然单位换成 `GeV`
- 截面 `cs` 在输出前被乘了 `10`

结合常用换算关系:

- `1 fm^2 = 10 mb`

可以得到明确判断:

- **Fortran 内部计算的 `cs` 更符合 `fm^2` 口径**
- **Fortran 写到 `cross_section.dat` 的第 4 列是 `mb`，不是原始 `fm^2`**

因此，如果讨论的是“旧 Fortran 输出文件的单位”，那么它不是 `fm^2`，而是已经换成了 `mb`。

### 10.2 Julia 主扫描脚本当前直接写出内部截面值

Julia 侧 `src/relaxtime/TotalCrossSection.jl` 的接口文档明确写明内部单位体系为:

- `s [fm^-2]`
- 各质量/能量量纲为 `fm^-1`
- `differential_cross_section(...)` 返回 `dσ/dt [fm^2]`

这说明 Julia 的总截面链路本身是在内部 `fm` 单位制下工作的。

进一步看 `scripts/relaxtime/scan_cross_section_vs_s_by_process.jl`，主扫描落盘时是:

- `σ = safe_sigma(process, s, params; n_points=opts.n_points)`
- `row = (..., sqrt_s_mev, σ, invσ2)`

这里没有看到任何 `σ * 10` 或 `MB_PER_FM2` 的换算，因此当前这两个主结果文件:

- `data/outputs/results/relaxtime/cross_section/xs_T210_muB0_xi0p0.csv`
- `data/outputs/results/relaxtime/cross_section/xs_T250_muB0_xi0p0.csv`

其 `sigma` 列更符合“**直接写出内部 `fm^2` 数值**”这一口径，而不是已经转换后的 `mb`。

### 10.3 调试脚本中的 `MB_PER_FM2` 进一步构成旁证

在 `scripts/debug/diag_usbar_mu0_channel_decomp.jl` 中，写 CSV 时表头显式使用了 `_mb` 后缀，并且逐列执行:

- `row.sigma_s_S * MB_PER_FM2`
- `row.sigma_s_P * MB_PER_FM2`
- `row.sigma_total_clamped * MB_PER_FM2`

这说明开发过程中已经默认承认:

- 内部截面变量本身是 `fm^2`
- 如果要输出成 `mb`，需要显式乘 `MB_PER_FM2 = 10`

这与上面的 Fortran 写盘逻辑完全同向。

### 10.4 因此，本轮单位结论是

用户提出的“可能是单位问题”这个方向是对的，但若严格区分“内部量”和“落盘文件”，更准确的说法应当是:

- **Fortran 内部截面是 `fm^2`，但 Fortran 输出文件已经转换成了 `mb`**
- **Julia 当前主扫描 CSV 很可能仍然直接写的是内部 `fm^2`，尚未转换成 `mb`**

这正好解释了前面观察到的稳定比例关系:

- `Julia / Fortran ≈ 0.1`
- `10 * Julia ≈ Fortran`

因为这里更像是:

- Julia CSV: `fm^2`
- Fortran `cross_section.dat`: `mb = 10 * fm^2`

### 10.5 这对当前问题的实际含义

到这一步，前面那条“整体差一个数量级”的结论已经可以从“数值现象”升级为“代码层面的单位链路结论”:

- Julia 和 Fortran 在 `usbar_to_usbar, mu_B=0, T=210/250 MeV` 上的曲线形状本身是对得上的
- 当前主差异首先是**输出口径不一致**，而不是散射振幅主体公式已经分叉
- 因而，若要继续做文献对比或 Julia-Fortran 自动回归，必须先统一输出单位，再谈残差阈值

## 11. 为后续 C++ 参考实现预备的复用流程

用户下一轮准备再引入一份 C++ 版本参考实现。为了避免重复摸索，本轮 Fortran 验证流程可以直接抽成以下检查单。

### 11.1 第一阶段: 先确认工程能否真实跑通

1. 确认构建方式，记录编译器、构建脚本、是否依赖 OpenMP / MKL / 第三方库。
2. 先运行默认主程序，确认它确实能启动、迭代并写出结果文件，而不是只做到“能编译”。
3. 记录默认主程序输出的物理量类型，判断它输出的是松弛时间、平均散射率，还是单过程 `sigma(sqrt(s))`。

### 11.2 第二阶段: 核对目标过程映射

1. 在 C++ 代码里定位 `usbar_to_usbar` 的过程编码、粒子编号、通道定义。
2. 确认它与 Fortran `usb_code`、Julia `:usbar_to_usbar` 是否是同一物理过程，而不是共轭过程或 flavor 排列不同。
3. 记录 `T=210/250 MeV, mu_B=0` 下使用的平衡态输入量和质量解来源。

### 11.3 第三阶段: 若默认主程序不输出同口径曲线，就补最小专项驱动

1. 只保留验证必需链路: 求平衡态、更新质量、调用单过程截面扫描。
2. 只跑 `T=210/250 MeV, mu_B=0` 两个点，避免把整套主程序扫描逻辑带进来。
3. 输出单独命名的结果文件，避免被默认结果覆盖。

### 11.4 第四阶段: 单位链路必须先查清楚

1. 分清内部单位和输出单位，分别记录。
2. 明确 `sqrt(s)`、`T`、`mu_B`、`sigma` 写盘前是否乘了 `hc`、`1000`、`10` 等换算因子。
3. 不接受只凭文件名或经验猜单位，必须找到最终写盘语句。

### 11.5 第五阶段: 做三层比较，而不是只做一层

1. 先和文献 digitized 数据比较，判断参考实现自身是否对得上文献。
2. 再和 Fortran 专项输出做逐点比较，判断是否与旧实现一致。
3. 最后和 Julia 在同一采样点逐点比较，区分“模型本身偏差”和“重写实现偏差”。

### 11.6 第六阶段: 固定产物格式，便于后续自动化

建议后续 C++ 验证也统一产出以下文件:

1. `state_summary`，记录 `Phi/Phibar`、夸克质量、收敛状态。
2. `cross_section` 曲线文件，至少包含 `sqrt(s)`, `mu_B`, `T`, `sigma`。
3. 与文献的误差摘要。
4. 与 Fortran / Julia 的逐点比较表。

这样下一轮引入 C++ 项目后，可以直接按本次 Fortran 的流程复刻，不需要重新设计验证口径。