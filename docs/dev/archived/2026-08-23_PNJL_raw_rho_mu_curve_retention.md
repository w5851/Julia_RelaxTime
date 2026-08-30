---
title: PNJL raw rho-mu curve retention
archived: true
original: docs/dev/active/2026-08-23_PNJL_raw_rho_mu_curve_retention.md
archived_date: 2026-08-23
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# PNJL raw rho-mu curve retention

状态：已完成，按三层策略归档。本任务只负责保留和恢复 C2 raw rho-mu 证据，不改变
equilibrium solver、Maxwell、hybrid policy、phase-reference promotion 或 transport。

## 三层保留策略

1. **仓库合同层**：`docs/analysis/pnjl/raw_curve_archive_v1/README.md` 定义字段、单位、
   精确 `(xi,T_MeV)` coverage、1281 个 rho 点、source manifest、hash 和恢复命令。
2. **不可变外部数据层**：Zenodo record `21980679` / DOI `10.5281/zenodo.21980679` 保存
   10458 条完整 byte-preserved 曲线；仓库只保留 pointer，不把约 1.88 GB ZIP 提交进 Git。
3. **恢复验证层**：`pnjl-raw-curve-archive-zenodo-restore.yml` 执行外层 ZIP hash、内部
   manifest、每条曲线/source-manifest hash、完整 coverage 和 sample restore 校验。

## 固定证据

- implementation merge：`94588784ab33f5936d1eae9eded6b8a7c1310e80`
- ledger registration：`16a52e323bd0c82b96ed663fe6a83a6cf20b43b5`
- audit run：`31941614867`
- source run：`32013771445`
- restore verification：`32039692479`
- calculation SHA：`3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`
- postprocess SHA：`67d73f871578e35759c08b3c75200c51646cf6cd`
- archive SHA-256：`467be7fb1075d1a5f0de3dd0d8afe29d9206a156c0ca7135a1e50967a4f18ccc`
- inner `archive_manifest.json` SHA-256：`514d9a7dd4cf537e8b209ed7df1cb996f52da48ab0b3672f27c3437d0cba4e52`

## DoD

- [x] 全量 exact-coordinate curve set 已验证：93 个 xi、10458 个 `(xi,T_MeV)`、每曲线 1281 个 rho 点。
- [x] 所有曲线保留完整 28 列 source schema；无插值、删点或重新求解的 archive-side 改写。
- [x] representative index 只引用 full curve index，不产生第二份代表性数据副本。
- [x] 外部 archive pointer、hash、provenance 和 restore path 已进入 Git。
- [x] restore verification 已完成；未来恢复必须重新执行 full-domain validation。

## 后续边界

任务已完成并归档；若未来需要新的 raw schema、额外坐标或更新计算 SHA，应创建新的
versioned retention task，不在本记录上覆盖历史证据。
