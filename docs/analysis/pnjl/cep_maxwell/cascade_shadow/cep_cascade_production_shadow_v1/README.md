# PNJL CEP cascade production shadow v1

本目录是 `main@d1e386deb950169b6d5db9cc3dd0617845ee60e8` 上 source run
`30527097410` 的只读 replay evidence。聚合/绘图代码来自
`d6ffa6b172d022157bea1ebb5118d5c500528dfc`；它不覆盖 reference，不启动
C0/C1/C2、C3/O1 或 transport。

verdict: `hybrid_required`。

- oracle anchors：稳定；24/24 anchor artifact、hash、finite/converged 合同通过。
- Actions critical path：608 s；source runner-minutes：87。
- cascade solver work：3927 fixed-rho requests、39592 residual calls、31582 Jacobian calls、139.36 s。
- memoized dense：15384 fixed-rho requests、130660 residual calls、99796 Jacobian calls、232.13 s。
- cascade fallback/retry：0；dense fallback/retry：0。
- 严格 gate 原因：`xi=0.5,T=5` 为 oracle ambiguous 但 cascade confirmed；`xi=-0.5` 没有 cascade confirmed-first-order 双端证据。

`method_costs.csv` 是由全部 anchor 的切片 cache 与 telemetry 重建的总量；原始 job
artifact 只通过 `curve_index.csv`、source hash 和 manifest 追溯，完整
`curve_points.csv` 保留在 Actions/local artifact。`manifest.json` 同时记录
source calculation SHA、postprocess SHA、source run ID 和 aggregation correction。

这份结果只支持“cascade 具有明显 solver-work 优势但需要 hybrid 收口”，不支持全域
cascade 直接晋升 production/reference。
