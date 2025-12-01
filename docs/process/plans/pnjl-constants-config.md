# PNJL 常量与参数管理草案

## 目标
- 将物理常数与可调模型参数解耦，方便在 RelaxTime 与 PNJL 求解器之间复用。
- 允许根据不同文献参数集快速切换，同时保留默认值。
- 在 CI/服务端可验证参数范围，避免误写导致的数值发散。

## 目录结构建议
```
config/
  pnlj/
    default.toml          # 官方默认参数
    lattice2020.toml      # 示例：来自特定文献的参数集
    custom.local.toml     # 开发者个人覆盖（可加入 .gitignore）
src/
  Constants_PNJL.jl       # 仅保留字段声明 + 读取逻辑
```

## 加载流程
1. `Constants_PNJL.jl` 定义 `struct PNJLParams` 以及 `load_pnjl_params(; profile::String="default")`。
2. 读取 `config/pnjl/{profile}.toml`，若未找到则 fallback 到 `default.toml` 内置字典。
3. 在模块末尾导出 `const DEFAULT_PARAMS = load_pnjl_params()`，同时提供 `set_global_params!(params::PNJLParams)` 供测试或服务器注入。
4. 服务器/脚本可以通过环境变量 `PNJL_PARAM_PROFILE` 指定配置文件，实现部署时零代码切换。

## TOML 示例
```toml
[name]
label = "PNJL-standard"
version = "2021-09"

[physical]
hbarc = 197.327         # MeV·fm
alpha_em = 1/137.035999

[model]
Nc = 3
Nf = 3
Lambda_MeV = 602.3
G_over_Lambda2 = 1.835
K_over_Lambda5 = 12.36
m_u_MeV = 5.5
m_s_MeV = 140.7

[polyakov]
T0_MeV = 210.0
a0 = 3.51
a1 = -2.47
a2 = 15.2
b3 = -1.75
b4 = 7.555
```

## 校验与日志
- `load_pnjl_params` 在解析后运行 `validate(params)`：检查耦合常数、Λ、T0 是否落在预期区间；若超界则抛出 `ArgumentError`。
- 将选用的参数集写入日志或 `/health` 响应，方便排查。

## 与 PNJL_Simulation 的接口
- `Function_PNJL_aniso.jl`、`Function_PNJL.jl` 先调用 `load_pnjl_params()` 拿到 `PNJLParams`，再将 `params.model.m_u_fm = params.model.m_u_MeV / params.physical.hbarc` 等值缓存到局部常量。
- 当 RelaxTime 作为上层 HTTP 服务时，可在启动脚本里先 `include("config/pnjl/loader.jl")`，再把 `PNJLParams` 传入求解器（避免两个项目各自复制常量）。

## 后续迭代
- 将不同模块（气液、旋转、各向异性）使用的参数字段标注在同一文件，必要时拆分 `pnjl-gas-liquid.toml`。
- 如果需要在运行时修改个别参数，可在配置中支持 `overrides` 小节或通过 API 提供校验过的补丁。
