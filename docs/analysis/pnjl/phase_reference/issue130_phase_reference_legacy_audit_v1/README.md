# Issue #130 PNJL legacy phase-reference retirement audit v1

这是 solver-free 的阶段 A 审计，不调用 Julia/PNJL，不修改 reference，不删除文件。
本次 repo HEAD：`ead55ece4a799948e54799288177f45908ff9d88`。

## Verdict

- verdict：`retirement_inconclusive`
- legacy fallback keys：`382`
- active consumer/rollback blockers：`24`
- candidate uncertified rows：`4074`
- physical deletion eligible：`False`
- next action：retain snapshot; migrate consumers and/or close candidate certified-key gaps before path retirement

## 关键语义

- 覆盖使用 adapter 的真实语义 key：boundary/spinodals 为 `(xi,T_MeV)`，
  crossover 为 `(xi,mu_MeV)`，CEP 为 `(xi)`；浮点 key 统一四舍五入到 12 位。
- candidate 的 unresolved/non-certified 行不会被视为 certified；因此 legacy fallback
  只要仍有缺键或未认证键就不能物理删除。
- 静态 consumer 扫描只覆盖当前 Git tracked source/config/test；不能证明外部生成器、
  本地脚本副本或未来代码不存在引用。

## 四表起点

| table | candidate rows | candidate certified | candidate uncertified | legacy rows | fallback keys |
|---|---:|---:|---:|---:|---:|
| boundary | 7162 | 3091 | 4071 | 48 | 35 |
| crossover | 1343 | 1343 | 0 | 336 | 315 |
| cep | 93 | 90 | 3 | 11 | 1 |
| spinodals | 6886 | 6886 | 0 | 57 | 31 |

## Consumer 处置

`tables/consumer_matrix.csv` 给出每个 tracked occurrence 的角色和是否阻塞物理清理；
`tables/claim_ledger.json` 保留机器可读的 claim/boundary 对象，CSV 供表格审阅。
当前共记录 `70` 个含 legacy token 的路径，其中 `24` 个仍需迁移或保留 rollback。

## 下一阶段

1. 先保留 `data/reference/pnjl/legacy_phase_reference_v1/`，不创建删除 PR。
2. 对 blocker consumer 设计 candidate-only/explicit-rollback contract，并增加 focused tests。
3. 迁移完成后重跑本 audit；只有 fallback=0 且未知 active 引用为 0 才能生成 physical-deletion proposal。
