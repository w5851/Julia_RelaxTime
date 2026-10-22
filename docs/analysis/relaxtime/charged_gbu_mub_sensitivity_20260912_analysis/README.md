# Charged GBU μB 灵敏度：历史诊断与来源留存

留存日期：2026-10-22。状态：`historical_diagnostic_retention`；原始数值与历史说明保持字节不变。
本文收口 `D:/w/charged-gbu-mub-sensitivity` 的既有研究，不启动新扫描。

固定 `T` 的 μB 扰动曾建议进一步研究联合 `T–μB` 冻结路径。该建议现由
`meson-freezeout-experiment / meson-freezeout-parameter-envelope` 承接：
[冻结线研究链](../meson_freezeout_defaults_20261017_analysis/README.md) 已随
[PR #327](https://github.com/w5851/Julia_RelaxTime/pull/327) 合入，
merge commit 为 `9c727d106abb9e931ba22dd26c206dc38d2efdff`。
这里保留 q0/有限 q 对照、失败记录和来源，不再重复安排物理路径研究。
冻结线参数带和图像仍按原研究的作者审阅状态处理。

## q0 与有限 q 的既有对照

[逐点对照表](q0_finite_q_comparison.csv) 使用有限 q 的
`muB_sensitivity_20260911_v3`，以及 q0 的
`muB_sensitivity_q0_20260913_probe3`、`muB_sensitivity_q0_20260913_5_77`。
九个点的 `T`、μB、夸克质量、flavor chemical potentials、介子化学势、
Polyakov 变量及 BQS 残差逐字段数值相同。背景均为 quark-only
`FixedMuBConservedCharges`、`rhoQ/rhoB=0.4`、`rhoS=0`，无介子反馈。

| √sNN (GeV) | 有限 q 中心 K+/π+ | q0 中心 K+/π+ | q0 相对有限 q |
| ---: | ---: | ---: | ---: |
| 3.0 | 0.509649 | 0.273127 | −46.41% |
| 5.0 | 0.538711 | 0.372656 | −30.82% |
| 7.7 | 0.545570 | 0.415422 | −23.86% |

这是相同已记录背景下两条计算路线的差异。有限 q 使用无限热 kernel 的稀疏诊断，
`q_nodes=8`、`qmax=8 fm⁻¹`、`mesh=48`、`cut_nodes=24`、
`omega_nodes=32`、`upper=32 fm⁻¹`；q0 使用 `q0_lambda_reference`，
`nq=24`、`qmax=8 fm⁻¹`、`mesh=512`、`nw=4800`、
`lower=1e-5 fm⁻¹`、`upper=56 fm⁻¹`、`thermal=24 fm⁻¹`。
完整设置以原 manifest 为准。路线与积分窗口不同，表中的差异不能直接解释成
离散化误差或 kernel 正确性的判定；没有全 q 尾部或零端点收敛认证。

两条路线的 ±50 MeV 记录均显示正 μB 响应。有限 q 的 3 GeV 加密结果
`0.376790 / 0.507245 / 0.686801` 保留在 `muB_sensitivity_20260912_refine3_v2`；
这只是原诊断的分辨率比较，不提升为 production baseline。

## 原始运行与失败

九个目录的全部 64 个原文件保存在
[raw_outputs.zip](../../../../data/outputs/results/relaxtime/analysis/charged_rpa_phase_backend/muB_sensitivity_retained_20261022/raw_outputs.zip)。
解压后使用原仓库相对路径；CSV、manifest、progress、PNG/PDF 和失败记录均未改写。

| 原运行目录 | 留存内容与限制 |
| --- | --- |
| `muB_sensitivity_20260911_probe3` | 首轮 3 GeV 诊断；记录的旧扫描脚本字节未恢复 |
| `muB_sensitivity_20260911_v1` | 九条有效 ratio；仅有 `status=running` 的 progress，无最终 manifest，不补造完成状态 |
| `muB_sensitivity_20260911_v2` | 九点结果与图；记录的旧扫描脚本字节未恢复 |
| `muB_sensitivity_20260911_v3` | 原报告引用的九点有限 q 结果；记录的源码和中心输入 hash 已核对 |
| `muB_sensitivity_20260912_refine3` | `upper=48 fm⁻¹` 时六个通道失败，原错误为 `ArgumentError: Bose-safe normal-gap geometry required`；三条 ratio 均为 NaN |
| `muB_sensitivity_20260912_refine3_v2` | 改用 `upper=32 fm⁻¹` 的三点加密诊断；记录的源码和输入 hash 已核对 |
| `muB_sensitivity_q0_20260913_probe3` | 3 GeV、±50 MeV q0 对照；历史 manifest 未记录扫描脚本 hash |
| `muB_sensitivity_q0_20260913_5_77` | 5/7.7 GeV、±50 MeV q0 对照；同样保留来源限制 |
| `muB_sensitivity_q0_20260913_target3` | 3 GeV 的 −270/0/+50 MeV 历史探针；非对称偏移缺少中心三元组，导数仍为 NaN、`missing_central_triplet`，不当作拟合或有效中心导数 |

原失败、未完成元数据、NaN 和 provenance 缺口都是证据的一部分。
有限 q 早期 `probe3/v1/v2` 中记录的三个旧扫描脚本 hash 与现存脚本不同；
[manifest](manifest.json) 逐项列出 recorded/captured hash，不以现存版本冒充旧版本。

## 源码、输入与恢复

[source_snapshot.zip](../../../../data/outputs/results/relaxtime/analysis/charged_rpa_phase_backend/muB_sensitivity_retained_20261022/source_snapshot.zip)
保存 capture HEAD `66413d6dacd3f23e458a0e7f70f6d2765861703c` 的现存工作树字节：
Julia 环境文件、`src/`、完整配置、`scripts/relaxtime/` 桥接依赖、相关 `causal_gbu_*.jl`、
有限 q 扫描脚本、原始分析说明和历史测试，以及三组已复用中心输入。
未提交的 q0 扫描脚本也在快照中。413 个文件的逐项 hash 与压缩包 hash 见 manifest。

历史 Git checkout 由已推送的 annotated tag
[`archive/charged-gbu-mub-sensitivity-20260912`](https://github.com/w5851/Julia_RelaxTime/tree/archive/charged-gbu-mub-sensitivity-20260912)
保存，指向上述 capture HEAD。tag 记录提交历史；未提交的 q0 脚本和原始结果仍以两个留存包为准。

有限 q 原 manifest 的 `source_git_head=b3922e5fc2e1f70150c812859c2464a8685a26f6`
是当时基点，扫描脚本当时尚未提交；已记录的源码 hash 才用于确认 v3/refinement
对应字节。q0 manifest 记录 capture HEAD，但缺少脚本 hash。
本次保存 q0 脚本可以恢复现存实现，不能追认其在历史运行时完全相同。

在一个空目录中解压两个包即可恢复各自的仓库相对路径，例如在仓库根执行：

```powershell
python -m zipfile -t data/outputs/results/relaxtime/analysis/charged_rpa_phase_backend/muB_sensitivity_retained_20261022/raw_outputs.zip
python -m zipfile -t data/outputs/results/relaxtime/analysis/charged_rpa_phase_backend/muB_sensitivity_retained_20261022/source_snapshot.zip
python -m zipfile -e data/outputs/results/relaxtime/analysis/charged_rpa_phase_backend/muB_sensitivity_retained_20261022/source_snapshot.zip D:/Temp/charged-gbu-mub-restored
python -m zipfile -e data/outputs/results/relaxtime/analysis/charged_rpa_phase_backend/muB_sensitivity_retained_20261022/raw_outputs.zip D:/Temp/charged-gbu-mub-restored
```

使用 manifest 的 SHA-256 核对包和解压文件；数值重放应使用恢复的历史环境，
不能替换为当前主干源码。原命令与参数位于恢复后的
`docs/analysis/relaxtime/charged_gbu_mub_sensitivity_20260912.md` 及各原 manifest。
快照用于历史恢复，不注册新的稳定 CLI 或测试。

解压目录不含 `.git`，原扫描脚本在生成运行 manifest 时会读取 Git HEAD；
若后续授权重放，可使用上述历史 tag 建立独立 checkout（例如
`git worktree add --detach D:/Temp/charged-gbu-mub-replay archive/charged-gbu-mub-sensitivity-20260912`），
核对 HEAD 等于 manifest 的 capture HEAD，再叠加两个留存包并校验文件 hash。
不能把单纯解压目录当作完整 Git checkout。

[validation.json](validation.json) 记录本次压缩包 CRC、477 个成员的字节/hash、
原报告三个 hash、九点背景对应和失败保留检查。已在空目录解压并实际加载 q0/有限 q
模块，确认桥接依赖可用；未调用数值求解器。

## 收尾边界

[PR #325](https://github.com/w5851/Julia_RelaxTime/pull/325) 已于 2026-10-08 合并，
merge commit 为 `7e5c7e896d4a4e66b73535ae11250ddfc74f1052`；合并 PR 的十项检查均成功。
[等高图任务](../../../dev/active/2026-09-30_charged-gbu-phase-contours.md) 的代码集成已收口，
M4 的作者研究产物审阅仍单独保留；不因代码合并或本次留存而标记 accepted/promoted。

留存包已随 [PR #328](https://github.com/w5851/Julia_RelaxTime/pull/328) 合入主干。
2026-10-22 经作者授权，历史 tag 已推送，原 `D:/w/charged-gbu-mub-sensitivity`
worktree 和 `codex/charged-gbu-mub-sensitivity` 分支已回收。
回收前，64 个原结果文件和未提交的 q0 脚本均与主干归档逐字节/hash 一致；
额外的 620 字节 `test/SuppositionDB/README.md` 为通用示例说明，无科研结果。
历史实现按 tag 与留存包恢复，M4 作者审阅和原有来源缺口维持上述状态。
