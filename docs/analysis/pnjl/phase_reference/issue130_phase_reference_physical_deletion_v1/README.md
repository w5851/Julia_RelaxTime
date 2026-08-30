# Issue #130 PNJL legacy phase-reference physical deletion v1

这是一个只删除工作树副本的独立 physical-deletion proposal。作者已授权创建
本 PR；合并仍需单独审核。删除对象是已经完成 canonical-path retirement、且不再被
accepted-primary runtime 使用的 `legacy_phase_reference_v1` snapshot。

## Scope

精确删除 allowlist 中的 8 个文件及其空目录：

`data/reference/pnjl/legacy_phase_reference_v1/`

不删除 `issue130_phase_reference_v1`、`issue130_phase_reference_v2` 的 strict/accepted/render
层、独立的 `data/reference/pnjl/crossover.csv`、任何 `docs/analysis` 历史 evidence、
raw rho-mu 外部归档指针、solver、Actions artifact 或 RS `prod_v2` 结果。

## Gates and provenance

- accepted 是 phase-reference 下游默认层；strict 只在显式模式使用；legacy 不再是
  runtime fallback 或 rollback。
- path-retirement merge：`9aa4c313901ca0c91e851f58514e3df9aa124df4`；该提交包含
  删除前的完整 legacy snapshot，可作为 Git 恢复边界。
- calculation SHA：`3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`。
- 删除前 tree SHA：`c128ef6358a5813533fc5a9726047585a750a4421f162293564bea8e363764e6`；
  8 files、50749 bytes。
- `solver_called=false`、`production_write=false`；本 PR 不重跑任何数值流程。

删除后的历史 hash、字节数、恢复命令和排除项均见 `deletion_manifest.json`、
`deletion_allowlist.csv` 与 `restore_commands.md`。合并后 fallback/rollback 不再可用；
恢复只允许在新的审核分支从 recovery ref 显式取回，不能作为静默 runtime 路径。
