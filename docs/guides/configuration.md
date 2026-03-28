# 配置管理指南

## 概述

本文档定义QCD模型库的配置管理策略，包括配置文件格式、参数管理、环境配置和最佳实践。

**目标**：提供灵活、可维护的配置系统，支持多种使用场景和环境。

## 阅读指引：实况层 vs 概念层

- **仓库实况（当前支持）**：以本仓 `config/` 真实文件与当前脚本入口为准，可直接复现。
- **概念示例（非现行实现）**：用于说明可选范式（如 JSON 配置、通用模板渲染、版本迁移脚手架），不代表仓库已默认启用。
- 若两者有冲突，**一律以仓库实况层为准**。

## 当前项目 Profile 口径（2026-02）

当前模型参数文件以 `config/models/` 为主：

- `config/models/njl/default.toml` + `config/models/njl/unittest.toml`
- `config/models/pnjl/default.toml` + `config/models/pnjl/unittest.toml`
- `config/models/pnjl/magnetic_default.toml`（外磁场 PNJL 配置模板）
- `config/models/rpnjl/default.toml` + `config/models/rpnjl/unittest.toml`

覆盖关系统一为：

1. 先加载 `default.toml`（完整基线）
2. 再把 `unittest.toml` 作为覆盖层 deep-merge 到 default
3. `unittest` 中未出现的字段，继承 `default`（即“只写差异”）

对应环境变量：

- `NJL_PARAM_PROFILE` / `PNJL_PARAM_PROFILE` / `RPNJL_PARAM_PROFILE`
- `PHYSICS_PARAM_PROFILE`（物理常量层）

推荐漂移检查命令：

```powershell
julia --project=. scripts/dev/check_model_profile_matrix.jl
```

该检查会验证：

- default+unittest 合并后关键字段完整；
- unittest 覆盖键必须是 default 已存在键（防止误拼写/静默新增漂移）；
- 输出各模型覆盖键列表，便于审阅 profile 差异。

## 最小可复现配置矩阵（models / physics / workflows）

| 配置域 | 最小文件集 | 主要消费入口 | 说明 |
|---|---|---|---|
| models | `config/models/pnjl/default.toml` + `config/models/pnjl/unittest.toml` | `src/models/Models.jl` | `default` 为基线，`unittest` 为差异覆盖 |
| physics | `config/physics/default.toml` + `config/physics/unittest.toml` | `src/models/Models.jl` / 相关 workflow | 跨模型共享常量层 |
| workflows | `config/workflows/relaxtime/default.toml` | `scripts/relaxtime/*` / workflow 入口 | 工作流级参数与 schema/alias 协同 |

推荐快速校验：

```powershell
julia --project=. scripts/dev/check_model_profile_matrix.jl
```

---

## 1. 配置层次

### 1.1 配置类型

QCD模型库使用三层配置结构：

```
1. 默认配置（代码中硬编码）
   ↓
2. 配置文件（TOML/JSON）
   ↓
3. 运行时参数（函数参数、环境变量）
```

**优先级**：运行时参数 > 配置文件 > 默认配置

### 1.2 配置范围

| 配置类型 | 位置 | 用途 | 示例 |
|---------|------|------|------|
| **模型参数** | `config/models/<model>/default.toml` | 物理参数、耦合常数 | G, K, Lambda |
| **物理常量参数** | `config/physics/default.toml` | 跨模型共享常量 | hbarc, rho0 |
| **环境配置** | `.env` 或环境变量 | 路径、日志级别 | DATA_DIR, LOG_LEVEL |
| **运行时配置** | 函数参数 | 特定计算的参数 | T, mu, xi |

---

## 2. 配置文件格式

### 2.1 推荐格式：TOML

**优点**：
- 人类可读
- 支持注释
- 类型明确
- Julia原生支持

**示例**：`config/models/pnjl/default.toml`

```toml
# PNJL模型默认参数配置
# 单位：自然单位 (fm^-1)，除非另有说明

[model]
name = "PNJL"
description = "Polyakov-loop extended Nambu-Jona-Lasinio model"

[parameters.coupling]
# 四夸克耦合常数 (fm^2)
G = 10.08
# 六夸克耦合常数 (fm^5)
K = -39.0
# 紫外截断 (fm^-1)
Lambda = 0.6

[parameters.masses]
# 流夸克质量 (fm^-1)
m_u0 = 0.0056
m_d0 = 0.0056
m_s0 = 0.135

[parameters.polyakov]
# Polyakov势参数
T0 = 0.19  # fm^-1
a0 = 3.51
a1 = -2.47
a2 = 15.2
b3 = -1.75

[numerical]
# 积分节点数
thermal_nodes = 64
vacuum_nodes = 128
# 求解器容差
solver_tol = 1e-10
solver_max_iter = 100

[physical_constraints]
# 参数有效范围
T_min = 0.0
T_max = 1.0
mu_min = 0.0
mu_max = 2.0
xi_min = -1.0
xi_max = 1.0
```

### 2.2 备选格式：JSON（概念示例，非现行实现）

**用途**：与外部工具交互、Web API

**示例**：`config/models/pnjl/default.json`（概念示例）

```json
{
  "model": {
    "name": "PNJL",
    "description": "Polyakov-loop extended Nambu-Jona-Lasinio model"
  },
  "parameters": {
    "coupling": {
      "G": 10.08,
      "K": -39.0,
      "Lambda": 0.6
    },
    "masses": {
      "m_u0": 0.0056,
      "m_d0": 0.0056,
      "m_s0": 0.135
    }
  }
}
```

---

## 3. 配置加载（概念示例，非现行实现）

> 以下代码块用于说明设计范式；当前仓库以 `src/models` 与既有脚本中的实际加载路径为准。

### 3.1 加载配置文件

```julia
using TOML

"""加载PNJL模型配置"""
function load_pnjl_config(config_file::String="config/models/pnjl/default.toml")
    if !isfile(config_file)
        @warn "Config file not found, using defaults" config_file
        return get_default_pnjl_config()
    end
    
    config = TOML.parsefile(config_file)
    return config
end

"""获取默认PNJL配置"""
function get_default_pnjl_config()
    return Dict(
        "parameters" => Dict(
            "coupling" => Dict(
                "G" => 10.08,
                "K" => -39.0,
                "Lambda" => 0.6
            ),
            "masses" => Dict(
                "m_u0" => 0.0056,
                "m_d0" => 0.0056,
                "m_s0" => 0.135
            )
        )
    )
end
```

### 3.2 配置合并

```julia
"""合并配置（运行时覆盖文件配置）"""
function merge_configs(base::Dict, override::Dict)
    result = deepcopy(base)
    for (key, value) in override
        if haskey(result, key) && isa(result[key], Dict) && isa(value, Dict)
            result[key] = merge_configs(result[key], value)
        else
            result[key] = value
        end
    end
    return result
end

# 使用示例
file_config = load_pnjl_config()
runtime_config = Dict("parameters" => Dict("coupling" => Dict("G" => 11.0)))
final_config = merge_configs(file_config, runtime_config)
```

### 3.3 配置验证

```julia
"""验证配置完整性和有效性"""
function validate_config(config::Dict, schema::Dict)
    # 检查必需字段
    for (key, spec) in schema
        if spec["required"] && !haskey(config, key)
            throw(ConfigurationError("Missing required field: $key"))
        end
    end
    
    # 检查值的范围
    for (key, value) in config
        if haskey(schema, key)
            spec = schema[key]
            if haskey(spec, "min") && value < spec["min"]
                throw(ConfigurationError("$key=$value below minimum $(spec["min"])"))
            end
            if haskey(spec, "max") && value > spec["max"]
                throw(ConfigurationError("$key=$value above maximum $(spec["max"])"))
            end
        end
    end
    
    return true
end

# 配置schema示例
const PNJL_CONFIG_SCHEMA = Dict(
    "G" => Dict("required" => true, "min" => 0.0, "max" => 100.0),
    "K" => Dict("required" => true, "min" => -1000.0, "max" => 0.0),
    "Lambda" => Dict("required" => true, "min" => 0.0, "max" => 10.0)
)
```

---

## 4. 模型配置

### 4.1 创建模型实例

**方式1：使用配置文件**

```julia
# 从默认配置创建
model = create_model(:PNJL)

# 从指定配置文件创建
model = create_model(:PNJL, config_file="config/models/pnjl/custom.toml")

# 从配置字典创建
config = load_pnjl_config()
model = create_model(:PNJL, params=config["parameters"])
```

**方式2：直接指定参数**

```julia
# 使用关键字参数
model = create_model(:PNJL, 
    G=10.08,
    K=-39.0,
    Lambda=0.6,
    m_u0=0.0056,
    m_d0=0.0056,
    m_s0=0.135
)

# 使用参数字典
params = Dict(
    :G => 10.08,
    :K => -39.0,
    :Lambda => 0.6
)
model = create_model(:PNJL, params=params)
```

### 4.2 配置文件组织

```
config/
├── models/
│   ├── pnjl/
│   │   ├── default.toml      # 默认参数
│   │   └── unittest.toml     # 单测覆盖参数
│   ├── rpnjl/
│   │   ├── default.toml
│   │   └── unittest.toml
│   ├── njl/
│   └── njl2/
├── physics/
│   ├── default.toml          # 共享物理常量
│   └── unittest.toml
└── workflows/
    └── relaxtime/
```

### 4.3 场景配置

**高密度场景**：`config/models/pnjl/high_density.toml`

```toml
[parameters.coupling]
G = 10.08
K = -39.0
Lambda = 0.6

[numerical]
# 高密度需要更多节点
thermal_nodes = 128
solver_tol = 1e-12

[scan_range]
T_min = 0.05
T_max = 0.3
mu_min = 0.5
mu_max = 1.5
```

---

## 5. 环境配置

### 5.1 环境变量

**支持的环境变量**：

| 变量名 | 用途 | 默认值 | 示例 |
|--------|------|--------|------|
| `QCD_DATA_DIR` | 数据目录 | `./data` | `/path/to/data` |
| `QCD_CONFIG_DIR` | 配置目录 | `./config` | `/path/to/config` |
| `QCD_LOG_LEVEL` | 日志级别 | `INFO` | `DEBUG`, `WARN` |
| `QCD_LOG_FILE` | 日志文件 | `stderr` | `/path/to/log.txt` |
| `QCD_CACHE_DIR` | 缓存目录 | `./cache` | `/tmp/qcd_cache` |

**使用示例**：

```julia
# 读取环境变量
data_dir = get(ENV, "QCD_DATA_DIR", "./data")
log_level = get(ENV, "QCD_LOG_LEVEL", "INFO")

# 设置环境变量（Julia中）
ENV["QCD_LOG_LEVEL"] = "DEBUG"

# 设置环境变量（Shell中）
# Linux/Mac
export QCD_DATA_DIR=/path/to/data

# Windows
set QCD_DATA_DIR=C:\path\to\data
```

### 5.2 .env文件

**创建**：`.env`（不提交到Git）

```bash
# QCD模型库环境配置
QCD_DATA_DIR=/home/user/qcd_data
QCD_CONFIG_DIR=/home/user/qcd_config
QCD_LOG_LEVEL=DEBUG
QCD_LOG_FILE=/home/user/qcd.log
```

**加载**：

```julia
using DotEnv

# 加载.env文件
DotEnv.config()

# 现在可以使用环境变量
data_dir = ENV["QCD_DATA_DIR"]
```

---

## 6. 配置最佳实践

### 6.1 配置文件管理

**DO（应该做）**：
- ✅ 提供默认配置文件作为模板
- ✅ 在配置文件中添加注释说明
- ✅ 使用版本控制管理配置模板
- ✅ 为不同场景提供预设配置
- ✅ 验证配置的完整性和有效性

**DON'T（不应该做）**：
- ❌ 在配置文件中硬编码路径
- ❌ 提交包含敏感信息的配置
- ❌ 使用过于复杂的配置结构
- ❌ 在代码中散布配置逻辑

### 6.2 参数命名

**命名约定**：
- 使用描述性名称：`thermal_nodes`而非`n`
- 包含单位信息：`T_MeV`或在注释中说明
- 使用一致的命名风格：`snake_case`
- 分组相关参数：`parameters.coupling.*`

### 6.3 默认值

**设置默认值的原则**：
- 默认值应该是"安全"的（不会导致错误）
- 默认值应该适用于最常见的场景
- 在文档中说明默认值的来源和适用范围
- 允许用户轻松覆盖默认值

### 6.4 配置文档

**每个配置文件应包含**：
- 文件用途说明
- 参数含义和单位
- 有效范围
- 默认值来源（如参考文献）
- 使用示例

---

## 7. 高级配置（概念示例，非现行实现）

### 7.1 配置继承

```toml
# config/models/pnjl/high_density.toml
[inherit]
base = "config/models/pnjl/default.toml"

[parameters.coupling]
# 只覆盖需要改变的参数
G = 11.0

[numerical]
thermal_nodes = 128
```

**实现**：

```julia
function load_config_with_inheritance(config_file::String)
    config = TOML.parsefile(config_file)
    
    if haskey(config, "inherit") && haskey(config["inherit"], "base")
        base_file = config["inherit"]["base"]
        base_config = load_config_with_inheritance(base_file)
        delete!(config, "inherit")
        return merge_configs(base_config, config)
    end
    
    return config
end
```

### 7.2 配置模板

```toml
# config/workflows/templates/scan.toml
[scan]
type = "{{SCAN_TYPE}}"  # TmuScan, TrhoScan
T_min = {{T_MIN}}
T_max = {{T_MAX}}
T_steps = {{T_STEPS}}
```

**使用**：

```julia
function render_config_template(template_file::String, vars::Dict)
    template = read(template_file, String)
    for (key, value) in vars
        template = replace(template, "{{$key}}" => string(value))
    end
    return TOML.parse(template)
end

# 使用示例
vars = Dict(
    "SCAN_TYPE" => "TmuScan",
    "T_MIN" => 0.1,
    "T_MAX" => 0.3,
    "T_STEPS" => 50
)
config = render_config_template("config/workflows/templates/scan.toml", vars)
```

### 7.3 配置验证Schema

```julia
using JSONSchema

# 定义JSON Schema
const CONFIG_SCHEMA = Schema("""
{
  "type": "object",
  "properties": {
    "parameters": {
      "type": "object",
      "properties": {
        "coupling": {
          "type": "object",
          "properties": {
            "G": {"type": "number", "minimum": 0, "maximum": 100},
            "K": {"type": "number", "minimum": -1000, "maximum": 0},
            "Lambda": {"type": "number", "minimum": 0, "maximum": 10}
          },
          "required": ["G", "K", "Lambda"]
        }
      },
      "required": ["coupling"]
    }
  },
  "required": ["parameters"]
}
""")

# 验证配置
function validate_config_schema(config::Dict)
    result = validate(CONFIG_SCHEMA, config)
    if !isvalid(result)
        throw(ConfigurationError("Invalid configuration: $(result.errors)"))
    end
    return true
end
```

---

## 8. 配置迁移（概念示例，非现行实现）

### 8.1 版本兼容性

**配置版本标记**：

```toml
[meta]
version = "1.0"
created = "2026-02-01"
description = "PNJL default configuration"

[parameters]
# ...
```

**版本升级**：

```julia
function upgrade_config(config::Dict, from_version::String, to_version::String)
    if from_version == "1.0" && to_version == "2.0"
        # v1.0 → v2.0: 添加新字段
        if !haskey(config["parameters"], "new_field")
            config["parameters"]["new_field"] = default_value
        end
    end
    
    config["meta"]["version"] = to_version
    return config
end
```

### 8.2 配置迁移工具

```julia
"""迁移旧格式配置到新格式"""
function migrate_config(old_config_file::String, new_config_file::String)
    old_config = load_old_format(old_config_file)
    new_config = convert_to_new_format(old_config)
    
    # 验证新配置
    validate_config_schema(new_config)
    
    # 保存新配置
    open(new_config_file, "w") do io
        TOML.print(io, new_config)
    end
    
    @info "Configuration migrated" old_file=old_config_file new_file=new_config_file
end
```

---

## 9. 配置示例（含概念示例）

> 本节含“生产模板/扫描模板”等说明性示例，不等价于仓库默认自动加载策略。

### 9.1 完整的PNJL配置

```toml
# config/models/pnjl/production.toml
# PNJL模型生产环境配置

[meta]
version = "1.0"
created = "2026-02-01"
author = "QCD Team"
description = "Production configuration for PNJL model"

[model]
name = "PNJL"
type = "isotropic"

[parameters.coupling]
G = 10.08      # fm^2
K = -39.0      # fm^5
Lambda = 0.6   # fm^-1

[parameters.masses]
m_u0 = 0.0056  # fm^-1
m_d0 = 0.0056  # fm^-1
m_s0 = 0.135   # fm^-1

[parameters.polyakov]
T0 = 0.19      # fm^-1
a0 = 3.51
a1 = -2.47
a2 = 15.2
b3 = -1.75

[numerical.integration]
thermal_nodes = 64
vacuum_nodes = 128
angular_nodes = 32

[numerical.solver]
method = "nlsolve"
tol = 1e-10
max_iter = 100
show_trace = false

[physical_constraints]
T_min = 0.0
T_max = 1.0
mu_min = 0.0
mu_max = 2.0
xi_min = -1.0
xi_max = 1.0

[logging]
level = "INFO"
file = "logs/pnjl_production.log"
format = "structured"

[output]
data_dir = "data/outputs"
format = "csv"
precision = 6
```

### 9.2 扫描配置

```toml
# config/scans/tmu_scan.toml
# T-μ扫描配置

[scan]
type = "TmuScan"
description = "Temperature-chemical potential phase diagram scan"

[scan.range]
T_min = 0.05
T_max = 0.3
T_steps = 50
mu_min = 0.0
mu_max = 1.5
mu_steps = 50

[scan.options]
xi = 0.0
adaptive_refinement = true
save_intermediate = true
parallel = true
num_workers = 4

[scan.output]
file = "data/outputs/results/pnjl/scan/tmu/tmu_scan_results.csv"
save_convergence_info = true
```

### 9.3 输出目录迁移说明（root `outputs/` -> `data/outputs/`）

- 当前默认口径：所有新脚本输出应写入 `data/outputs/`（按 `results/`、`figures/` 分层）。
- 推荐细分：
    - 性能基准产物写入 `data/outputs/perf/<module>/`；
    - 参数扫描产物写入 `data/outputs/results/<module>/scan/<feature>/`；
    - 导数量化产物写入 `data/outputs/results/<module>/derivatives/`；
    - 回归校验产物写入 `data/outputs/results/<module>/regression/`；
    - 调试比对产物写入 `data/outputs/results/debug/<topic>/`；
    - 图像产物写入 `data/outputs/figures/<module>/`。
- 历史兼容口径：根目录 `outputs/` 仅用于读取历史产物，不再作为默认落盘路径。
- 迁移建议：
    - 旧命令若显式写 `--output outputs/...`，请改为 `--output data/outputs/...`。
    - 若历史脚本需复现实验，可保留旧入参路径，但新产物应回写到 `data/outputs/`。
- 空目录说明：`data/outputs/results/debug/baseline_compare/` 作为对比产物标准目录，若当前没有执行过 Fortran/Julia 对齐导出任务则可为空，不视为迁移遗漏。

---

## 10. 故障排查

### 10.1 常见问题

**问题1：配置文件未找到**

```julia
# 错误
ERROR: Config file not found: config/models/pnjl/custom.toml

# 解决
# 1. 检查文件路径是否正确
# 2. 使用绝对路径
# 3. 检查工作目录
pwd()  # 查看当前目录
```

**问题2：配置参数无效**

```julia
# 错误
ERROR: ConfigurationError: G=150.0 above maximum 100.0

# 解决
# 1. 检查参数范围
# 2. 查看配置文档
# 3. 使用默认值
```

**问题3：配置合并冲突**

```julia
# 问题：运行时参数未生效
# 原因：配置合并顺序错误

# 正确做法
final_config = merge_configs(file_config, runtime_config)
# 而非
final_config = merge_configs(runtime_config, file_config)
```

### 10.2 调试配置

```julia
# 打印当前配置
function print_config(config::Dict, indent=0)
    for (key, value) in sort(collect(config))
        if isa(value, Dict)
            println("  "^indent, key, ":")
            print_config(value, indent+1)
        else
            println("  "^indent, key, " = ", value)
        end
    end
end

# 使用
config = load_pnjl_config()
print_config(config)
```

---

## 11. 检查清单

配置新模型时：
- [ ] 创建默认配置文件
- [ ] 添加参数注释和单位
- [ ] 定义参数约束
- [ ] 实现配置加载函数
- [ ] 实现配置验证
- [ ] 编写配置文档
- [ ] 提供使用示例
- [ ] 添加配置测试

---

## 12. 参考资料

- TOML规范：https://toml.io/
- Julia TOML包：https://github.com/JuliaLang/TOML.jl
- 数据契约：`docs/api/data_contracts.md`
- 错误处理：`docs/guides/error_handling.md`
- 参数类型API：`docs/api/PARAMETER_TYPES_API.md`

---

**最后更新**：2026-02-01  
**维护者**：开发团队
