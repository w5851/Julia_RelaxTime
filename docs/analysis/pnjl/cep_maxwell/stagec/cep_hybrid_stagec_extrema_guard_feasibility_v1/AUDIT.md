# Stage-C extrema-guard feasibility audit

输入 manifest、curve、source run、calculation SHA 和 replay producer SHA 均写入
`manifest.json`。所有 S-shape、Maxwell 和 geometry 结果由当前 Julia PhaseCore
重新计算；历史 oracle 标签只在 replay 完成后用于 gate。Stage-C 分类曲线始终是
完整 Stage-B 全域曲线与选定局部点的并集。

成本按 Stage-A、完整 Stage-B 和选定 Stage-C rho key 的并集计算，Stage-B 网格不
免费消费。离线 replay 不能证明 residual/Jacobian、fallback/retry 或 runner-minutes；
这些必须由后续 Actions targeted/full shadow 验证。source Actions 的 aggregate
physical gate 失败被保留为诊断事实，但所有输入曲线均需 finite/converged。
