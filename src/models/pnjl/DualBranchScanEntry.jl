"""\
    DualBranchScanEntry

独立入口（按需加载）：把 `PNJL.DualBranchScan` 从主模块默认导出中移出后，
如果你需要做一阶相变区域的双分支扫描，可直接 include 本文件来启用。

用法：

```julia
include("src/models/pnjl/DualBranchScanEntry.jl")

# 之后 DualBranchScan 位于 PNJL 子模块内
using PNJL
PNJL.load_dual_branch_scan!()

res = PNJL.DualBranchScan.run_dual_branch_scan(; T_mev=100.0, mu_range=0.0:10.0:400.0)
info = PNJL.DualBranchScan.find_phase_transition(res)
```

说明：
- 本文件会确保 `PNJL` 已加载，然后调用 `PNJL.load_dual_branch_scan!()`。
- 不会把 DualBranchScan 的符号 re-export 到 `PNJL` 顶层命名空间。
"""

const _INCLUDE_ONCE_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "utils", "IncludeOnce.jl"))
if !isdefined(Main, :IncludeOnce)
    Base.include(Main, _INCLUDE_ONCE_PATH)
end
const IncludeOnce = Main.IncludeOnce

const _PNJL_PATH = normpath(joinpath(@__DIR__, "PNJL.jl"))
IncludeOnce.include_once!(Main, :PNJL, _PNJL_PATH)

Main.PNJL.load_dual_branch_scan!()

