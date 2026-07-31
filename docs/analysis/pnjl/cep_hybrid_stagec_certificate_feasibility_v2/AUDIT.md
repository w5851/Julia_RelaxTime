# Stage-C certificate feasibility v2 audit

输入 manifest SHA、source run ID、calculation SHA 和本地 replay producer SHA 均写入
`manifest.json`。所有最终 curve 证书均由 Julia 核心重算；本审计不把历史 oracle row 的
geometry 直接当作 hybrid 证书。完整 raw `curve_points.csv` 继续留在 Actions/local artifact，
仓库只保留 selected-point index 与代表性曲线。

成本按 Stage-A、完整 Stage-B 和 selected Stage-C rho key 的并集计算，Stage-B semantic
grid 不免费计入。离线 replay 不能证明 residual/Jacobian/runner-minutes，须由后续 Actions
targeted/full shadow 验证。
