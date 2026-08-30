# Audit boundary

输入包括 base final aggregate 与五点 tolerance-contract aggregate。两套 manifest、文件
哈希、曲线键、finite/converged 状态和双 SHA provenance 在执行时验证；共同曲线数值字段
另有 `curve_identity_audit.csv` 记录。rho=0 的化学势端点不是唯一可识别的
thermodynamic root，重放保留该差异为 provenance warning，并以 tolerance artifact 覆盖
重叠端点；正密度曲线才用于数值 identity gate。所有 Maxwell/geometry 证书由当前 Julia PhaseCore
与 outer geometry comparison 重算，内部 solver tol 与外层 gate 分离。

本次是 solver-free replay：不估算 residual/Jacobian/runner-minute 成本，不复制 raw
`curve_points.csv` 到仓库。若 verdict 为 `feasible_candidate`，仍必须由 targeted/full
Actions shadow 验证，再由作者审核物理曲线和成本。
