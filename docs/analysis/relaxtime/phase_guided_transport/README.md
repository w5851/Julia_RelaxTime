# Phase-guided transport

## 当前结果

当前正式图层为 **publication_clean_v15**，唯一选择入口是
[current](publication_clean_current.json)。
[v15 正式交付](../../../../data/outputs/results/relaxtime/transport/phase_guided/publication_clean_v15_pdf/README.md)
包含 76 张 PDF，直接引用已接受的 76 张彩色和 76 张灰度 PNG。
论文已接入三张输运组合图；正文的科学解释与审阅状态由论文的 `writing_package/00_status.md` 管理。

数值源为[η 加密正式网格](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_eta_augmented_20261029_v1/README.md)：
1867 点，1810 个前驱原值、9 个原坐标重算值和 48 个新增 ξ 点。
图中使用实际计算值，不平滑、不插值。精度结论覆盖已经执行的分层、联合代表点、
完整网格有效性和局部积分检查，不等于所有点的统一误差上界。

- 图件与复现规则：[绘图合同](plotting_case_contract.md)。
- 计算、验收与研究边界：[统一状态](../threshold_cache_upgrade_20261021_analysis/production_plan_20261023.md)。
- π、K、η、σ 的局部结构：[Mott 分析](../threshold_cache_upgrade_20261021_analysis/mott_local_bumps_20261028.md)。

## 旧资产退役

2026-11-06 按作者意见，让被取代的完整旧包退出当前分支：v1–v13 的图像、
显示点表、孤立 manifest、重复源码 ZIP、旧 raw/收敛表及只服务这些版本的生成器和测试一并退役。
只证明旧版本自身仍完整的交叉引用不构成保留理由；退役不新增全量备份或逐文件清单。
v14 和最初 full-RS 试排的重复图件及审阅副本也已退役。

**v14 是已验收的 1819 点前驱，不能称为未收敛结果。**
以下内容仍直接服务当前结果，继续保留：

- v15 的数值包、晋升与验收依据、PNG/PDF、冻结源码和实际字体定位记录。
- v15 使用的 v14/raw 绘图模块，以及新采样 1819 点前驱和独立收敛资格。
- v14 的版式 manifest：v15 的 PNG 生成器仍读取它，v15 自身也保存其冻结副本。
  该 manifest 记录的旧图需从 Git 恢复后才能复核；它不是当前图层入口。
- v15 先前指针明确绑定的少量历史元数据，以及支撑现有结论的 Mott 结果和回归基线。

旧 v1–v13 数值不用于当前物理预言。需要调查旧算法误差时，从 Git 恢复对应版本，
在独立 checkout 中复现；当前基线中的小规模受控比较已记录旧 κ/λ 漂移的收敛原因。

本次退役前的完整版本为
`6a3d08663c845587ae457fbe1c90173e0e5dcc25`；例如可用
`git show 6a3d08663c845587ae457fbe1c90173e0e5dcc25:<repository-relative-path>`
查阅某份旧记录，或在独立 checkout 恢复整套历史入口。
更早已退役的 557 张 PNG 见提交 `bad5f98ce99ffd81a9ff0f60336a6c3c32a3d44a`。
Git 历史保留原始证据，不承诺当前分支仍可运行每个历史生成器。
