# Julia 模块加载与依赖

本仓库采用 include 驱动的组织方式。用户入口为 `Models` 与 `src/models/entrypoints.jl`；模块职责和依赖方向见 [架构规则](../architecture/dependency_rules.md)。

## include、using 与 import

- `include` 在指定模块的作用域中执行文件，可以加载同一模块的代码片段，也可以由聚合模块加载子模块。
- `using` 引入已有模块及所选名字；扩展其他模块的函数时用 `import`。
- 同一模块的代码片段共享作用域，无需互相导入。加载顺序要保证顶层表达式用到的类型、常量和方法已经定义。
- 子模块有独立作用域，使用相对导入表达依赖；不需要为了减少文件数把所有子模块合并。

```julia
module Example
include("types.jl")
include("kernels.jl")
include("workflow.jl")
end
```

## 本仓库的加载边界

聚合入口负责加载实现，调用方优先使用聚合后的接口。库文件内部的路径使用 `@__DIR__` 定位；稳定 CLI 可以从仓库根目录启动。

现有共享模块的 `Main` 单例和 guarded include 需要保留，以免产生两份类型或重复定义。沿用现有边界时可采用：

```julia
if !isdefined(Main, :ParameterTypes)
    Base.include(Main, parameter_types_path)
end
using Main.ParameterTypes: QuarkParams, ThermoParams
```

其中 `parameter_types_path` 应由文件位置构造为确定路径。这是 include 驱动仓库的加载约定，不要求所有新 helper 都经 `Main` 注册。普通子模块优先由其所属聚合模块加载，并使用相对导入。

`module` 声明必须位于顶层，不能直接放进普通 `if` 块。需要入口幂等性时，优先在调用方对既有模块做加载检查。

## 拆分与依赖

- 按职责拆分；同一个物理内核保留一份实现。目录或 API 外形一致不要求为每个模型建立空的占位模块。
- 稳定入口与内部实现分开。内部文件迁移时更新调用者和当前文档，不复制整份实现维持旧路径。
- 遇到循环依赖，先检查模块归属；确实共享的类型或纯函数可以下沉。仅为消除一个循环而增加通用框架通常没有必要。
- `Base.invokelatest` 只用于已确认的动态加载边界，普通调用直接派发；登记规则见 [架构规则](../architecture/dependency_rules.md)。

## 验证

修改加载路径或模块归属后，检查统一入口能够加载、导出函数仍指向预期实现，并运行受影响的接口及工作流测试。数值行为可能变化时补充对应回归，不通过放宽容差处理差异。

常用入口：

```sh
julia --project=. scripts/dev/check_models_entry_contract.jl
julia --project=. scripts/dev/check_pnjl_migration_guard.jl
```

具体测试选择见 [命令参考](agent_command_reference.md) 与 [测试治理](testing_governance.md)。
