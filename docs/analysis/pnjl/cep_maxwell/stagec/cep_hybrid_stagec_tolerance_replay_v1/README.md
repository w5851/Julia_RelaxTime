# PNJL Hybrid Stage-C tolerance-contract replay v1

verdict: `oracle_inconclusive`。本目录将已完成的 24-anchor hybrid shadow final aggregate 与
Maxwell tolerance-contract 五点 revalidation 做 solver-free 合并重放；没有调用
equilibrium solver，不修改 production、reference 或历史 v1/v2 evidence。

- base source run: `30601857594` / `fde2b929b60575f1daacb84a1b9b8ff6e3b0a0cc`
- tolerance revalidation run: `30730990835` / `467be1fce847a9c991ec362c3335be07fccbe604`
- internal Maxwell solver tolerance: `5.0e-6`
- outer Maxwell geometry gate: `5.0e-5`
- tested caps: `12, 24, 48, 96, 160`
- simulated dense reference: `15384.0` unique solves
- solver called: `false`

五个 revalidated anchor 的局部 certificate 与 base Stage-A/B 证据分开记录；新深层
`0.0015625` 曲线只作为独立收敛证据，不被免费计入 hybrid cost。该 replay 仍是
candidate/diagnostic evidence，不能直接晋升 production/reference。
