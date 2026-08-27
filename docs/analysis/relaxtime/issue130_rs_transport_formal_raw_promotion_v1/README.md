# Issue #130 RS transport formal raw promotion v1

本包记录作者已审核后的两套 RS `prod_v2` 结果（mode-A 与 mode-B）登记为
`approved` formal raw，以及分析入口默认迁移到 `prod_v2` 的不可变审计边界。
它不是新的数值重跑，也不是 RS 全域 numerical convergence 证明。

## 结论

- 两套结果分别位于：
  - `data/outputs/results/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2/`
  - `data/outputs/results/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2/`
- 作者已接受两套结果及同构审核图的 raw/formal layout；registry 中状态为 `approved`，但
  `manuscript_eligible=false`。
- 分析脚本 `build_phase_guided_transport_xi001_jump_analysis.py` 与
  `phase_guided_p128_mechanism_scan.jl` 的默认 case 已切到 `prod_v2`，并使用版本化的
  `_analysis_v2` 输出目录；旧 `prod_v1` 分析和结果树保持不变。
- 本次 promotion/default migration 没有调用 solver、没有写 production CSV、没有切换
  runtime workflow 默认，也没有删除 legacy fallback 或旧 `prod_v1`。

## 证据文件

- `manifest.json`：输入哈希、两套结果/图像 manifest、registry 状态和不变边界。
- `input_inventory.csv`：逐 mode 的 row、hash、质量警告和状态摘要。
- `claim_ledger.csv`：可支持的 raw promotion claim、范围和未完成事项。

数值来源为 calculation SHA
`3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`、workflow head
`22874505877491754eed27519ad8a7b871c82571`。source artifact 的 numerical audit 仍为
`diagnostic_only`；其质量警告和历史 sidecar 缺陷在原始 result manifest 中保留。

## 后续边界

分母/传播子近零点和非一阶相变突变点的文献显示清理必须另建 versioned
`publication_clean` 派生层，不能从本包所指向的 raw CSV 中删除或覆盖点。legacy fallback、
显式 `--phase-reference-mode legacy` rollback 以及旧 `prod_v1` 结果继续保留，直到另行审查
并授权 retirement。
