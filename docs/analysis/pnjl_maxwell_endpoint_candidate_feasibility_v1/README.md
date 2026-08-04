# PNJL Maxwell endpoint/candidate feasibility v1

verdict: `feasible_candidate`。本包只对已有 Actions 曲线做 solver-free replay；没有调用
equilibrium/Newton solver，没有修改 production、reference 或历史 evidence。

严格 candidate 要求：

- 面积函数只在恰有三个去重交点、且切片具有 `+→−→+` S-shape topology 时有效；
- 无效拓扑会重置 sign-change 状态，不能跨两个交点/三交点边界二分；
- 同一曲线所有有效 sign-change 都会枚举，多个稳定根保持 inconclusive；
- `rho=0` 只记录端点代表值，不复制重复点，也不插入 `(0,0)`。

当前低温证据中的旧 Maxwell 残差约 `36--38` 是端点/两交点伪变号；严格 replay
应在约 `331.55 MeV` 找到唯一三交点根。该结论仍是离散曲线 candidate，不能替代
后续正密度近零补点的跨分辨率证书。

输入：D:\Desktop\Julia_RelaxTime_issue130_artifacts\cep_hybrid_extrema_guard_20260803\required_three_deep_run_30809754119\cep-deep-oracle-cep_hybrid_extrema_guard_required_three_deep_20260803-aggregate, D:\Desktop\Julia_RelaxTime_issue130_artifacts\cep_hybrid_extrema_guard_20260803\targeted_replay_run_30808473818\cep-hybrid-production-shadow-cep_hybrid_extrema_guard_targeted_replay_20260803-aggregate-replay。
