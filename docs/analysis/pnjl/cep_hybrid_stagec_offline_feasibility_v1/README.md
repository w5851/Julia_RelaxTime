# PNJL Hybrid Stage-C offline feasibility v1

verdict: `maxwell_candidate_inconclusive`。本目录只使用 aggregate replay CSV 做离线重放，不调用 equilibrium solver，不改变 production、reference 或历史 evidence。

- source run: `30601857594`
- source calculation SHA: `fde2b929b60575f1daacb84a1b9b8ff6e3b0a0cc`
- tested caps: `48, 96, 160, 224`
- classification mismatches across cap evaluations: `8`
- selected cap: `None`

Stage-C classification uses the complete oracle-level-0 0.00625 curve as an offline semantic reference, merges selected local 0.003125 points, enumerates all stable +→−→+ candidates, and uses the existing oracle geometry row only as a non-solver geometry gate. The oracle curve is not charged to simulated hybrid cost.

若 verdict 不是 `feasible_candidate`，不得创建 production PR；当前结果只用于定位弱 S、Maxwell 候选或成本问题。
