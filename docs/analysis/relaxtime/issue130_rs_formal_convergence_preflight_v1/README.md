# Issue #130 RS formal numerical convergence preflight v1

生成日期：2026-08-30

工作区：`D:/Temp/julia_relaxtime_issue130_shared`

分支：`codex/issue130-rs-formal-convergence`

当前 `origin/main`：`cc4f71a47f41d3eba20509af0f8409c75d7b9a6d`

## 结论

本包是 solver-free、Action-free 的输入与 provenance 预检，不是数值收敛 verdict，也没有在本机调用 PNJL equilibrium solver。

现有两套 p128 `prod_v2` artifact 的 calculation SHA 都是
`3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`，但 source workflow head 都是
`22874505877491754eed27519ad8a7b871c82571`。该 head 的 phase-reference CLI 默认是
`strict`，并通过旧 adapter 保留 candidate-strict → legacy fallback；当时仓库还没有
accepted v2 runtime。故：

* p128 可以作为 `strict` 显式口径的高阶端复用；
* p128 不能直接作为当前 `accepted-primary` 默认 runtime 的高阶端；
* 从当前 main 只改 p104 的网格输入，会生成 accepted-primary p104，和旧 strict p128
  混配，不能宣称是同一 reference 下的 formal convergence；
* 因此本预检没有触发 p104 Action。

## 两条可行路线

### 推荐：accepted-primary matched pair

1. 固定 calculation SHA `3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`、同一 mode-A/mode-B
   15+15 shard topology、同一 direct-coexistence `xi=±0.003` 语义。
2. 在当前 accepted-primary runtime 下先重跑 p128，确保高阶端与目标默认 reference 完全一致。
3. 在同一 workflow contract 和同一 accepted-primary reference 下补 p104：
   `104/18/30/20/440`。
4. 重用现有 convergence evaluator 的字段集合和状态/finite/negative/failure gates，
   但把 phase-reference layer/mode、source SHA、workflow SHA 写入新的 result-side
   convergence manifest。

这是唯一能证明“当前默认 accepted-primary runtime 的 p104→p128 数值收敛”的路线。

### 可选：strict-only historical pair

若作者明确只需要验证 strict 显式口径，可在历史 workflow ref
`22874505877491754eed27519ad8a7b871c82571` 上，仅改变 p104 的五个积分输入，保持
strict 默认与旧 source provenance；现有 p128 可以复用。这条路线成本较低，但结论只适用于
strict-explicit，不覆盖当前 accepted-primary 默认。

## “是否需要新增 Action 文件”的回答

对 p104 网格本身，不需要新增 workflow 文件：现有
`.github/workflows/relaxtime-phase-guided-transport-production.yml` 已有
`tau_p_nodes`、`tau_angle_nodes`、`tau_phi_nodes`、`tau_n_sigma`、`sigma_grid_n` 输入，
可表达 `104/18/30/20/440`，并能按 mode 和 xi 列表分 shard。

但“只改输入”不能修复 p128 与 accepted-primary 的 reference provenance 不一致。若选择
推荐路线，应先让 p128/p104 都在 accepted-primary 下生成（必要时新增一个显式记录
`phase_reference_layer=accepted`、`phase_reference_mode=runtime` 的小型 workflow contract
变更），再运行数值 gate。

mode-A 的 direct-coexistence 特殊点不得删掉：`xi=0` 在严格共存 anchor 上由已认证的
`xi=-0.003/+0.003` 两侧点表示；p104 和 p128 必须使用相同的有效 key 拓扑。

## 输入证据

| 证据 | 结论 |
| --- | --- |
| p128 mode-A/B manifest | source solver 已运行；aggregate/import 均未调用 solver；`numerical_status=diagnostic_only`；两套 manifest 与 effective config 哈希见 `input_inventory.csv`。 |
| source workflow head 228 | CLI 默认 `phase_reference_layer=:strict`、`phase_reference_mode=:runtime`；adapter 只有 v1 strict/derived/render 与 legacy fallback。 |
| current main cc4 | CLI 默认 `phase_reference_layer=:accepted`、`phase_reference_mode=:runtime`；accepted v2 manifest 已标记 `accepted_for_downstream`。 |
| 现有 workflow | p104 五个积分输入可由 dispatch 参数表达；workflow 当前未显式暴露 phase-reference layer/mode。 |
| historical p104↔p128 gate | 仅是旧 `first_canonical_v1` case 的 production-grade 历史对照，不是当前 prod_v2 accepted-primary 证明。 |

## 运行边界

本包未修改任何 result/reference/figure，未调用 solver，未 dispatch Action，未改变容差、Maxwell、equilibrium solver、transport kernel 或 direct-coexistence 路由。下一步只有在 formal scope 明确后才可产生数值 Action。
