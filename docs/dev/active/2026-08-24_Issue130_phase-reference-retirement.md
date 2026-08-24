# Issue #130：legacy phase-reference canonical 路径退役

状态：active。runtime switch、consumer adapter、solver-free parity 和限定 numerical pilot 已完成，
作者已接受新 reference 的限定结果并授权旧 reference retirement 与后续 RS transport 重验。

## 范围与语义

- canonical runtime source 保持
  `data/reference/pnjl/issue130_phase_reference_v1/strict`；只消费 certified 行。
- 历史 `boundary.csv`、`cep.csv`、`crossover_dense.csv` 和 `spinodals.csv` 不再占用
  `data/reference/pnjl/` 根路径，迁入 versioned
  `data/reference/pnjl/legacy_phase_reference_v1/` snapshot。
- snapshot 保持原始字节和 SHA-256，只用于 candidate 缺键/未认证键的逐键 fallback、显式
  `--phase-reference-mode legacy` rollback，以及历史结果复现。
- 本任务不是物理删除：strict candidate 仍需要 boundary/crossover/CEP/spinodal fallback，只有后续
  reference 覆盖完整且所有 consumer 不再需要 rollback 时，才可另行评估删除 snapshot。
- `data/reference/pnjl/crossover.csv` 是另一组历史 fixed-point 输入，不属于本次 dense legacy
  retirement，保持原路径不变。
- 已完成且 hash-bound 的 import/runtime-audit evidence builder 保留原始根路径字符串；其历史包
  应在对应 source SHA 上重放，不用当前 HEAD 伪造同名历史 manifest。

## 实现与验收

- [x] 以 `git mv` 迁移六个 legacy 文件，新增 byte/sha256 retirement manifest 和说明。
- [x] Julia transport、phase-guided、Paper/Python、TaylorDiff 与 workflow 消费者改读 snapshot；
  PNJL phase-diagram workflow 对当次生成的根路径文件使用显式参数，避免误读 snapshot。
- [x] 增加 snapshot hash、canonical-root absence、fallback/rollback 与消费者路径测试。
- [x] 运行 focused Julia/Python、workflow schema、task ledger、active docs、docs consistency、
  script/data-path governance 和 `git diff --check`。
- [ ] PR 通过 CI 后合并；从合并 SHA 准备新的 candidate-runtime RS transport production，
  不覆盖旧 production。

## 已有放行证据

- runtime switch：PR #260，merge SHA
  `1ccf29310fb20c30bcd154f0b4966e25a7565225`。
- solver-free parity：PR #261。
- paired numerical pilot：run `32684074876`，candidate/legacy 各完成 19/19 点；共同的五个
  `tau_u_ubar_ratio_high` 质量警告不构成 candidate-specific regression。
- quality audit 修复：PR #263，merge SHA
  `d4193bcf77bb23d740787a103782e3e6fc96dbba`。
- 作者结论：新 reference 效果符合预期，授权 retirement 和 RS transport 后续运行。

## 非目标与停止条件

- 不修改 equilibrium solver、transport kernel、Maxwell、phase-reference 数值或任何容差。
- 不把限定 pilot 称为全域 RS convergence；全量重验必须使用新 case、不可变 provenance 和独立
  production audit。
- 若 snapshot hash、fallback/rollback 或任何消费者契约失败，停止 retirement，不触发 RS production。
