# Fortran_Julia Legacy Reuse Audit

## Scope and provenance

本文档审计 `Fortran_Julia_full.zip`，判断早期 Fortran 到 Julia 的 NJL 输运移植是否仍有内容可以安全进入 `Julia_RelaxTime` 主线。本文记录的是历史复核，不把归档中的报告、输出或“完成度”表述当作当前生产证据。

- 来源归档：`D:\Desktop\_cleanup_quarantine\2026-05-27\_safety_exports\Fortran_Julia_full.zip`
- SHA-256：`15898DC6381E5EC869E390ABC6565765D9D2D885CEAAF096335F5BE72EC88E53`
- 归档规模：152 个条目；解包后 147 个文件、1,007,381 bytes
- 审计解包目录：`D:\Desktop\_cleanup_quarantine\2026-05-27\_audit_tmp\Fortran_Julia_full`
- 审计日期：2026-08-17

## 内容概况

归档包含一个早期 include 驱动的 NJL-Julia 工程，主要覆盖：

- Gauss-Legendre、Cauchy 主值和其它数值积分工具；
- 极化函数、传播子、散射振幅和夸克方程求解；
- 平均散射率、弛豫时间、剪切粘滞和电导率的实验性入口；
- 大量阶段报告、比较报告、临时脚本和测试输出。

这些模块名称与当前主项目有明显重叠，但重叠本身不构成可迁移性。旧工程依赖本地 include 顺序、早期参数结构和隐式全局状态，不能直接接入当前 `Models` 与 `src/relaxtime/` 合同。

## 可靠性证据

### 报告之间不一致

`Julia/PROJECT_COMPLETION_REPORT.md` 同时声称“100% 完成度”和“测试通过率 60%”。报告列出的 30/30 或 9/9 收敛点只能说明某些试验扫描完成，不能证明散射截面、散射率和输运系数已经通过物理验证。

`Julia/output/calculation_summary.txt` 也只记录了 30 个计算点全部收敛；它没有给出可复现的 Git 提交、输入指纹、公式版本或独立参考值。

### 散射率仍没有生产级证据

`Julia/Stage_I_D_Scattering_Rate_COMPLETION_REPORT.md` 记录了从 `sigma_ij = 1.0e-3` 占位值切换到 `amplitude2` 的尝试，但同一报告的示例输出仍为：

```text
w_ij (归一化): 0.0 GeV
w_ij_n (未归一化): 0.0 GeV
```

报告把该零值解释为“技术修复完成”，并把参数优化和 Fortran 对照列为后续工作。因此这不是可直接迁移的数值结果。

代码中仍可见未完成或仅用于诊断的路径：

- `Julia/src/types.jl` 的默认物理常数带有“占位，需要替换”注释；
- `Julia/src/main.jl` 将输运系数标为占位实现；
- `Julia/src/numerics/transport.jl` 在异常/无效路径返回零散射结果；
- `Julia/src/numerics/scattering_rate_optimized.jl` 的散射截面缓存键只由 `process_code` 构成，未编码质量、化学势、Polyakov 参数、温度和网格设置，不能作为当前生产缓存契约。

### 测试和输出不能直接成为 baseline

归档包含许多名称为 `test_*` 的脚本和成功日志，但报告本身已经暴露出通过率不一致、占位实现和零值回退。输出目录没有当前主项目要求的 manifest、配置指纹、失败语义和回归基线契约。因此不迁移这些测试结果、DAT/CSV 输出或性能数字。

## 当前主项目对应能力

当前主项目已经有更近且职责清晰的实现边界：

- `src/relaxtime/RelaxTime.jl`：弛豫时间入口；
- `src/relaxtime/MesonPropagator.jl`：介子传播子相关逻辑；
- `src/relaxtime/ScatteringAmplitude.jl`：散射振幅；
- `src/relaxtime/TotalCrossSection.jl`：总截面；
- `src/relaxtime/AverageScatteringRate.jl`：平均散射率；
- `src/relaxtime/TransportCoefficients.jl`：输运系数；
- `Models` 和分层 unit/integration/regression/validation 测试：入口、配置和验证边界。

旧工程没有提供一个可证明优于这些当前模块的独立算法或经过当前参数契约验证的实现。直接复制反而会引入旧单位、旧接口、旧缓存语义和错误回退。

## 可保留的启发

只有以下方法层想法值得在未来任务中重新设计，而不是复制代码：

1. 用显式的 2x4 过程矩阵表达散射通道；
2. 将截面预计算与主散射率积分分开；
3. 用 Fortran 对照点驱动分阶段验证。

如果未来重新使用这些想法，应在当前参数类型、缓存指纹、单位约束、失败语义和 regression/validation 分层下重写，并先建立可复现参考值。

## 决定

- 直接迁移代码：不执行。
- 迁移旧输出或测试结果：不执行。
- 立即可做的主项目贡献：保留本审计记录，避免将旧“完成报告”误当作生产证据。
- 未来若需要增强输运验证：新建当前主线任务，使用现有 `src/relaxtime/` 合同和独立 Fortran/reference 数据，而不是恢复该旧工程。

本项已完成彻底清理：原解包目录、`Fortran_Julia_full.zip` 和本项没有其它主线依赖的残留均已删除；主项目只保留本审计记录。
