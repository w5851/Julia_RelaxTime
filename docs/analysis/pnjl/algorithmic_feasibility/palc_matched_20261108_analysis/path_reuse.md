# PALC 路径复用与校正成本

状态：算法与数值验收完成（2026-11-12），交付 PR [#346](https://github.com/w5851/Julia_RelaxTime/pull/346)。
最终算法源码 `36bce656e286d1bfc3195b4d23c4172ccf94b550`；基准为 PR #345 的端点 PALC 与自然密度连续扫描。
结果集中保留在 [path_reuse_optimization.json](../../../../../data/outputs/results/analysis/palc_matched_20261108/path_reuse_optimization.json)，约 1.8 MiB。原始逐次日志和重复状态数组留在 Actions 临时 artifact；本文件及紧凑结果足以恢复配置、根数、误差、计时和结论。

## 最终行为

- PNJL 的 Tμ、对称/非对称 Tρ、FixedSigma（固定 s/nB）及 BQS 路径默认使用 Newton 目标校正。先检查残差，最多 8 步、每步最多 5 次减半回溯；失败从原初值恢复 trust-region，最终仍检查原残差与硬约束。
- NLsolve 的函数、Jacobian 准备及求解器工作区在单条路径内复用。参数改变时清空值/Jacobian 缓存；结果复制后返回。目录种子搜索仍采用 trust-region。
- PALC 默认在同一单调段内使用四点三次目标预测，邻接回弯或不规则间距时使用线性预测；每个交点独立校正。
- PALC Jacobian 默认中心有限差分；ForwardDiff 为显式选项。两者均复用缓冲。AD 在不同机器上的性能排序不同，不统一设为默认。
- 对称 Tρ 扫描和密度会话默认复用邻温最高密度端点初值；Tμ/非对称 Tρ 保留显式选项。每个新温度仍重新求解整条路径。
- FixedEntropy 和 RPNJL 的默认预测、校正及 Jacobian 选择保持原设置。独立 residual-first 开关已并入 Newton。积分公式、节点、权重及非 fixedmu 八维联立约束不变。

选项见 [solver API](../../../../api/models/solver/Overview.md)、[Tμ API](../../../../api/models/scans/TmuScan.md)、[Tρ API](../../../../api/models/scans/TrhoScan.md)。

## 性能与工作量

Julia 1.12.5、单线程/单 BLAS 线程、24×8 tensor 积分、iterations=300。每方法每轮 5 次路径求解，3 轮轮换方法顺序；首轮含编译，仅后两轮用于 warm 中位数。时间包括全部候选求解和相同的选中结果重构；验证和序列化在外。对照在同一份新代码中显式恢复 PR #345 的数值控制，双方共享新的残差/FD 缓冲；因此倍率反映控制器与预测改进，未单独计入共享缓冲相对旧提交的收益。最终矩阵每批计时前统一 GC，批内 GC 保留在总时间中。

| 场景 | 原设置 / 新默认 | 分配下降 |
| --- | ---: | ---: |
| Tμ 升降温序列，ξ=0/0.4 | 1.361–1.376× | 15.1–15.4% |
| 对称 Tρ 升降温序列 | 1.279–1.283× | 22.7–22.8% |
| 非对称 Tρ 升降温序列 | 1.159–1.163× | 16.8–17.3% |
| T=120 MeV 细 Tμ，89 个目标/111 根 | 1.321× | 约 17% |
| FixedSigma 三个窗口（最终 FD） | 1.147–1.281× | 约 7–12% |
| BQS（最终 FD） | 1.316× | 约 13% |
| FixedEntropy（原算法） | 0.998× | 基本相同 |

温度序列为 T=100,110,120,130,140,160,180,200 MeV，升/降各一遍。独立粗 Tμ 为 μq=0:15:450 MeV；细窗口 296:0.25:318 MeV，含 11 个三根目标。密度网格为 0 加 0.2:0.15:3.5；FixedSigma 为 s/nB=12.0759，T=100:2:140（ξ=0/0.4）及 120:2:160；BQS 为 μB=900 MeV、T=100:3:160。

细 Tμ 的根 Jacobian 次数由 387 降至 253，迭代由 275 降至 252；最终 Actions 根阶段约 159.6→112.7 ms，PALC 约 19.1→19.0 ms，减少终止点的额外 Jacobian 是主要收益。密度路径可能触发一次局部 Newton 恢复，全部目标仍收敛。端点复用单独贡献较小：首批序列 Tμ 0.989–1.009×，对称 Tρ 1.023–1.070×，非对称 Tρ 0.972–1.055×，所以仅对称密度默认启用。

### Jacobian、积分与内存

Tμ 微基准中，就地残差从每次 96 bytes 降至 0，耗时约 29.7 μs；FD/AD Jacobian 都为 0 bytes，约 298/353 μs。FixedSigma 就地残差约 1216 bytes，FD/AD Jacobian 约 19456/6432 bytes；最终 runner 的 AD 耗时 1.87 ms，FD 1.62 ms。AD 在远机和首批 Actions 的 FixedSigma 更快，表明不能从残差调用次数或单机结果推断所有机器的收益。

FixedSigma 的 AD/FD 最大绝对差约 2.95e-5，Jacobian 最大元素约 85.3，相对范数差 2.93e-7；差分步长减半后误差约缩小四倍，符合中心差分截断误差。Tμ 相对差约 1e-10。微基准与步长检查由同一比较入口生成。

远机 T=120 MeV、μq=310 MeV 的抽样剖析中，620 个样本的 617 个经过 thermal omega，522 个经过 calculate_log_sum；主要成本是嵌套 Dual 的乘法及 exp/log，节点和权重已经缓存。此处仅定位热点，未改变热积分公式。下一层积分/AD 算术优化应由具体负载另行驱动。

## 数值与简并选态

原容差保持 state=1e-4、thermodynamics=1e-6、residual=1e-6。完整矩阵有 120 条路径/切片、720 条方法/轮次记录；最后 FD 补充矩阵 5 条路径、30 条记录。每条记录核对全部候选根、各根热力学量及最大压力选态。同机原设置/新默认均严格一致，包含稳定、亚稳和中间分支。

完整矩阵跨方法最大差：状态 2.74e-8、压力 4.86e-10、密度 3.64e-8、熵 2.25e-8、能量 4.09e-8；参考最大残差 7.29e-8。单独路径与 PR #345 冻结根逐一匹配。FixedEntropy 有 21 个目标、42 个电荷共轭根，默认仍不晋升为生产 PALC。

历史 FixedEntropy 选中密度符号在 7 个目标与 PR #345 不同：两根压力只差约 1e-14，根集合状态差不超过 4.45e-15。固定熵约束不规定密度符号，μ/ρ 变号同时交换 Phi/Phi_bar 后给出等压根，舍入会改变 argmax。报告保留 `strict historical selected=false`；只有两根完整、满足共轭关系且压力差在 64 ulp 内，才接受历史简并说明。其它模式、丢根或非简并压力差仍失败。同机新旧比较未作该例外，求解器选择规则和数值容差也未修改。

实际消费者测试在远机通过 3047 个断言，覆盖默认 Tμ/Tρ、密度会话、T120 三根、S 形/Maxwell 结构、BQS、isentropic 介子入口、CSV/callback 顺序与 resume。Python 比较器 30 项通过。最终源码的 13 项常规 PR 检查通过。

## 运行来源与复现

| 批次 | 源码 | 结果 |
| --- | --- | --- |
| [扩展分项 43155017475](https://github.com/w5851/Julia_RelaxTime/actions/runs/43155017475) | 84e192e0 | 1800 s 预算退出，完成 119 路径/切片、1239 记录；已完成数值全部一致，FixedEntropy 未完成 |
| [完整默认矩阵 43165548215](https://github.com/w5851/Julia_RelaxTime/actions/runs/43165548215) | e1a9c6a7 | 120 路径/切片、720 记录全部完成；原汇总因历史固定熵简并报失败，修正汇总后通过 |
| [最终 FD 43169410396](https://github.com/w5851/Julia_RelaxTime/actions/runs/43169410396) | 36bce656 | 5 路径、30 记录全部通过，workflow success |

首批 growing JSON 每方法重复写盘造成额外开销；现在按方法组写盘并记录预期矩阵。缺失、重复或未计划的记录不能通过完整验收。早期 T200 ξ=0.4 与 sigma_t160 的组合总耗时还受到前序工作触发的 GC 影响，最终矩阵统一了批前 GC；历史总耗时保留，未事后扣掉 GC 制造提速。

```text
PALC_TIMING_SAMPLES=5 julia --startup-file=no --project=. scripts/analysis/pnjl_bifurcation_spike/run_path_reuse_comparison.jl output.json 3 all baseline,default
python scripts/analysis/pnjl_bifurcation_spike/compare_path_reuse.py --source output.json --output reviewed --summarize-only
```

Windows 在运行前设置 `$env:PALC_TIMING_SAMPLES="5"`。`--methods baseline,default` 可通过 Python CLI 或 Actions 输入设置；最终 FD 小矩阵的 cases 为 `sigma,sigma_xi4,sigma_t160,bqs,entropy`。首批分项严格复现需使用对应源码版本。约束温度微基准可同时输出 AD/FD 与步长差异。

本轮未重算完整生产。磁场验证由独立对话负责；PALC 收尾后进入 PR3 的实际低温消费者及目标导数收敛核验。
