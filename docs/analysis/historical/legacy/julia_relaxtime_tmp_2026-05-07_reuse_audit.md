# Julia_RelaxTime_tmp_2026-05-07 Legacy Reuse Audit

## Scope and provenance

本文档审计一份 2026-05-07 的临时输出归档，判断其中是否还有未被主项目或外部文献库覆盖的内容。

- 来源归档：`D:\Desktop\_cleanup_quarantine\2026-05-27\_safety_exports\Julia_RelaxTime_tmp_2026-05-07_full.zip`
- SHA-256：`34DB170472C62E481B705FC94E683FA570EDDEA92C4C5B019733E778D72116AC`
- 归档规模：12 个条目，包含 4 篇 PDF、4 个文本头和 4 个 CSV
- 审计解包目录：`D:\Desktop\_cleanup_quarantine\2026-05-27\_audit_tmp\Julia_RelaxTime_tmp_2026-05-07_full`
- 审计日期：2026-08-17

## 文献副本

四篇 PDF 与 `D:\Desktop\paper` 中的文献文件逐字节一致，哈希如下：

| arXiv 标识 | SHA-256 | 外部文献库中的对应文件 |
|---|---|---|
| `1608.05383` | `8CBA7BE55E2299E033A923C1C36D8318F02301D35BFB0568ADC54AEBF9FAEC43` | `D:\Desktop\paper\docs\papers\PNJL\Blaschke 等 - 2017 - Mott dissociation of pions and kaons in hot, dense quark matter.pdf` |
| `1305.3907` | `196FBD1C9B07705F6F5642F725E0D09A89B18D8617D8241FC8328FB8621BF09D` | `D:\Desktop\paper\docs\papers\介子数密度(及其比值)\Blaschke 等 - 2014 - Generalized Beth--Uhlenbeck approach to mesons and diquarks in hot, dense quark matter.pdf` |
| `2301.09882` | `B54777DA443D241B2C208CD8A1BE0587DAB6FC530B15658743B6E2BC60A2DF9A` | `D:\Desktop\paper\docs\papers\PNJL\Maslov和Blaschke - 2023 - Effect of mesonic off-shell correlations in the PNJL equation of state.pdf` |
| `1310.4960` | `BD02699EC1BC7ADA5B055FB04B3E59250C1C83B5A4BAE80D26DE0846697E94A1` | `D:\Desktop\paper\docs\papers\介子压强引入EOS状态方程\Yamazaki和Matsui - 2015 - Quark-hadron phase transition in a three flavor PNJL model for interacting quarks.pdf` |

因此删除临时 PDF 副本不会丢失文献材料，也不需要将这些二进制文件复制进主项目。

## CSV 结果审计

`freezeout_phase_shift_gbu_baseline_extract.csv` 是旧的 freezeout 摘要表，哈希为
`43EB32838A71F736A25232657C2CDD82C63C667E891A1D9340DA2F3850EBBE89`。它没有 command、Git SHA、参数指纹或 manifest；同一主题的当前主项目 baseline 位于
`tests/baselines/relaxtime/baseline_meson_density_freezeout_phase_shift_gbu_path_v1.csv`，并由
`tests/regression/relaxtime/test_meson_density_freezeout_phase_shift_gbu_path_regression.jl` 约束。两者哈希不同，故旧表不能直接覆盖当前 baseline。

三个 `r3_probe*.csv` 是诊断 probe，而不是可晋升的生产结果。其字段中可见 `m_pi_MeV`、`m_K_MeV` 为 `NaN`，`gamma_pi_MeV`、`gamma_K_MeV` 为零，且 `equilibrium_iterations` 为 `-1`。这类状态不能作为介子质量、宽度或 freezeout 物理结论。

当前主项目已经提供更完整的对应边界：

- `docs/guides/sop/workflows/meson_density.md` 规定 command、remote manifest、收敛摘要和 figure manifest；
- `docs/api/models/scans/` 与 `docs/api/relaxtime/` 记录稳定入口；
- `tests/unit/`、`tests/integration/`、`tests/regression/` 和 `tests/validation/` 覆盖介子密度、介子质量和 Mott 相关 workflow；
- `data/outputs/` 中的当前结果按 case、manifest 和 provenance 组织。

## 决定

- PDF：不迁移，外部文献库已有逐字节相同副本。
- CSV：不迁移，不把无 manifest 的旧 probe 或摘要表写入当前输出和 baseline。
- 立即可做的主项目贡献：保留本审计记录；当前没有需要继续运行才能得到的可靠结果。

本项已完成彻底清理：原 ZIP、审计解包目录及本项没有其它主线依赖的残留均已删除；外部文献库和主项目当前 workflow 不受影响。
