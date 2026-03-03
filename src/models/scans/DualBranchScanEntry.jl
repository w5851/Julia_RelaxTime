"""
    load_dual_branch_scan!() -> Module

Models 域双分支扫描入口：保留按需加载语义，统一返回 `DualBranchScan` 模块。
"""
function load_dual_branch_scan!()
    return DualBranchScan
end
