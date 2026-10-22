# 默认化学冻结线比较的实验输入

`experimental_ratios.csv` 对应用户投影片引用的 Universe **9**, 106 (2023),
Fig. 1 左图。共 22 个中心重离子点；不含 pp、PHSD 或 SMASH。

- 14 个 Au+Au 发布值作为主比较：AGS 5 点、STAR BES 5 点、STAR 14.5 GeV 1 点、STAR 高能 3 点。
- NA49 的 5 个 Pb+Pb 中快度点辅助比较 horn，和 Au+Au 分开评分。
- ALICE 的 2 个点是 `(K+ + K-)/(pi+ + pi-)`，`|y|<0.5`、0–5%，仅辅助显示高能平台。
- HADES 的 2.4 GeV 点从综述原图读取，约 0.010；取点范围只表示读图误差，图中不把它画成实验误差棒，也不用于选择默认方案。

原始 HEPData JSON、STAR Table IX 的表格转录与复核图片在 `raw/`。
URL、DOI、SHA-256 和筛选细节见 `sources.json`。后处理入口为
`scripts/analysis/relaxtime/prepare_freezeout_experiment_data.py`；重建时指定新 `--output`，
避免覆盖已冻结的实验 CSV。发布表中的数字保持来源精度；HADES 读图值只应解读为约 0.010。

## 元数据勘误与核验边界（2026-10-17 问答复核）

本 CSV 的 `NA49_17.27` 中心度被统一写成了 `0-7.2%`，这是元数据错误。
原始 NA49 2002 文献（PRC 66, 054902，[nucl-ex/0205002](https://arxiv.org/abs/nucl-ex/0205002)，
PDF p.20 Table I）给出：40A、80A GeV 为 0–7.2%，158A GeV 为 **0–5%**；
158A GeV 对应此处的 17.27 GeV 点。该点发布比值 0.179 和绘图纵坐标不因这项勘误改变，
它也不属于 Au+Au 主评分集合。

该原文 PDF 已补存为文献目录中的 `S16_NA49_2002_PbPb_pion_kaon.pdf`，SHA-256：
`af24ec3a2ccad39fb337aa7226204d7ec2d525d5f5fb1e9152e111dda1ce0a69`。
此前对旧 SPS 点沿用的快度字段来自汇编的全局 qualifier；原始 2002 文献 p.8 对其
中快度 `dn/dy` 使用 `|y|<0.6` 描述。后续修订实验表时应按能量回溯原文，不能把
汇编的全局 qualifier 当作所有原始实验均已满足的统一接受度。

本次重新从冻结来源重建的 22 行与 CSV 全字段一致，并复核了 STAR Table IX 的三个
转录值和 HADES 的坐标校准。这证明提取过程可复现；不证明编纂元数据已无误。
旧 CSV、原始 JSON 和图像保持原字节，生成器也仍对应旧版元数据；这份勘误是当前
解读入口。修正后的机器可读表及生成器需要作为后续版本处理，不覆盖旧图的输入。

注意以下口径：

1. NA49 Table 10 是中快度数据，Table 9 是全相空间数据；这里取前者。表内低于 6 GeV
   的点来自 Au+Au AGS，不能沿用表的全局 `Pb Pb` 标签。保留表中 2.35 GeV 的能量值，
   没有重新换算或把它移到另一个碰撞能量。
2. STAR 2009 的 HEPData Fig. 30a 记录实际只有负电荷比值；高能正电荷数据来自原文
   PDF p.27 Table IX：62.4、130、200 GeV 分别为 `0.162±0.012`、`0.167±0.011`、`0.159±0.019`。
   130 GeV 的中心度为 0–6%，其余为 0–5%。
3. STAR 14.5 GeV 从同一中心度的 `26.7/141` 得到比值。误差传播假设两产额协方差为零，
   明确标为派生误差；不据此给出正式拟合优度。
4. ALICE 2020 的 `syst. uncorr.` 是跨中心度不相关的系统误差部分，不再与完整 `syst.`
   重复相加（arXiv:1910.07678，Sec. 4.1，积分产额的误差说明）。原分量仍在 JSON 中。
5. 原综述低能 AGS 点及不同实验有各自的接受度/系统误差；本轮等权 RMSE 只用于相对选择，
   不把混合实验集合当作具有已知协方差的精密拟合。

`hades_digitized/` 保留原图 crop、取点 overlay、坐标锚点和读图范围。
从本目录的持久 crop 重放可用 `hades_digitize_replay.json`（输出到新 `hades_replay/`），
不依赖初次准备图像时的临时路径。原始提取 metadata 保持不变。
