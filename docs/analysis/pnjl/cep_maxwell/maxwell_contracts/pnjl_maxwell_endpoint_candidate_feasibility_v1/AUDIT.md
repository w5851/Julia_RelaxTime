# Audit boundary

输入曲线的 SHA、行数、calculation SHA 和 producer head 写入 `manifest.json`。所有派生
表只使用已有 `(method,xi,T,rho)` 记录；没有插值、补点或 equilibrium 重算。`legacy_*`
字段只记录当前 public Maxwell 的对照结果，不作为新候选 gate。

`endpoint_dependent=true` 表示左外交点落在首个 rho 单元内，提示后续 Actions 需要正密度
局部细化；它不是零密度平台的硬编码，也不是 production promotion 证书。
