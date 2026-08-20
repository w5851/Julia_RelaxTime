# C2 phase surfaces v6: v5 baseline + crossover endpoint overlay

本包明确以 v5 后处理结果为基线，而不是直接从 C2 原始 CSV 重建。v5 的 Maxwell boundary、spinodal、CEP bracket、物理 crossover 筛选、unresolved 诊断和无三角化规则均保留；本包只追加 endpoint expansion 的 crossover 结果。

固定 calculation SHA 为 `3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`。endpoint expansion 来自 numerical run `32240898122` 的 aggregate replay，replay run provenance 位于 `tables/endpoint_replay_manifest.json`；本地物化过程 `solver_called=false`、`reference_write=false`、`oracle_labels_consumed=false`。

v5 基线包含 Maxwell `6886` 行、spinodal `6886` 行、保留 crossover `1157` 行和 `93` 个 CEP bracket。新增 overlay 为 `186` 行，覆盖 `93` 个非均匀 xi 切片；它不是完整均匀 xi 网格，也不代表 Maxwell 区已补全。

图中不使用插值或三角化；超过显式 gap 上限的 native support 保持空白。endpoint 点以空心圆标识，便于审核其相对 v5 端点的延伸，但该标识不表示新的物理证书。物理上 CEP 是 crossover 与 first-order 的共同端点，因此后续 `interpolated_noncertified` 派生层可以把选定的 CEP boundary estimate 作为 crossover 的终止边界；v6 本身尚未执行这一步，`mu_q > mu_CEP` 的响应峰仍不绘制为 crossover。

本包 verdict 为 diagnostic-only，不能触发 phase-reference 晋升、正式 production 或 transport。Maxwell 区数值扩展仍是独立任务。
