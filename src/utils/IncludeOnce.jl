"""IncludeOnce

极轻量工具：在 include 驱动（非 package 化）项目中，提供“只 include 一次”的样板封装。

动机：
- 同一 Julia session 内，多处脚本/测试可能重复 `Base.include(Main, path)`。
- 重复 include 会触发 `WARNING: replacing module ...`、world-age 噪声、以及类型分裂。

约定：
- `target` 通常为 `Main`（作为模块单例锚点）。
- `sym` 是目标模块名，例如 `:PNJL`、`:ThermoFacade`。
- `path` 应当是 `normpath(joinpath(...))` 得到的稳定绝对路径。
"""
module IncludeOnce

export include_once!

"""include_once!(target, sym, path) -> Module

若 `target.sym` 尚未定义，则 `Base.include(target, path)`，然后返回 `target.sym`。

注意：
- 本函数只负责“幂等加载”。模块内部是否 include-once（`@eval module ... end`）由被 include 的文件自行保证。
"""
function include_once!(target::Module, sym::Symbol, path::AbstractString)::Module
    if !isdefined(target, sym)
        Base.include(target, path)
    end
    return getfield(target, sym)::Module
end

end # module IncludeOnce
