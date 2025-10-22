# Julia RelaxTime

## 【重要单位约定】
```
=============================================================================
本模块统一使用自然单位制 (ℏ = c = 1)，所有物理量单位为 fm⁻¹：
- 温度 T: fm⁻¹ 单位 (转换：T[MeV] = T[fm⁻¹] × ℏc = T[fm⁻¹] × 197.327 MeV·fm)
- 化学势 μ: fm⁻¹ 单位 (μ[MeV] = μ[fm⁻¹] × 197.327)  
- 质量 m: fm⁻¹ 单位 (m[MeV] = m[fm⁻¹] × 197.327)
- 四动量 k₀, k: fm⁻¹ 单位 (p[MeV] = p[fm⁻¹] × 197.327)
- 能量变量 z: fm⁻¹ 单位 (E[MeV] = E[fm⁻¹] × 197.327)
- 偏振函数 Π: fm² 单位 (真空极化振幅，量纲 = 1/k²)
- BPM积分 B: fm² 单位 (泡泡图积分振幅)
- Polyakov环参数 A: 无量纲
- 截断参数 Λf: fm⁻¹ 单位 (由常数模块提供，Λf ≈ 3.05 fm⁻¹)
=============================================================================
```

## 项目说明

本项目用于计算弛豫时间相关的物理量。

## 安装与使用

1. 激活项目环境：
```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

2. 使用模块：
```julia
include("src/relaxtime/relaxtime.jl")
```

## 项目结构

下面列出了仓库中主要文件和文件夹及其用途，帮助快速定位代码和文档：

- `Project.toml`, `Manifest.toml`：Julia 项目环境与依赖清单，用于激活和重现环境。
- `src/`：主要源代码目录。
	- `integration/`：数值积分相关实现，如 `CauchyPV.jl`, `GaussLegendre.jl`。
	- `relaxtime/`：包含弛豫时间计算实际实现代码（模块内部实现）。
- `api/`：接口与公式的说明文档（Markdown），例如 `CauchyPV.md`, `GaussLegendre.md`。
- `doc/`：项目文档（草稿与模板），如 `formula/` 目录下的公式说明。
- `prompt/`：编码规范和变量命名指南（Markdown 与 YAML），用于团队协作时保证风格一致。
- `results/`：运行或测试生成的结果文件（数据、图表等）。
- `test/`：单元测试文件，用来验证 `integration` 与 `relaxtime` 的正确性，如 `test_cauchypv.jl`, `test_gausslegendre.jl`。
- `README.md`：本文件，项目介绍与快速上手说明。
