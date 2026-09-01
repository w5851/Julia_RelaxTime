# 定向相位路线检索补充（2026-09-01）

本补充不是新的全文综述，而是为下一 PR 的“相位归一化--S-matrix--测度”闭合做的三组
定项核对。查询、日期、数据库和限制已登记在 `tables/search_log.csv` 的 Q11--Q13；
这里只记录能改变代码设计的原始来源结论。

## 1. Retarded propagator 与 on-shell S-matrix

### 已能唯一确定的部分

Dashen--Ma--Bernstein 的 S-matrix 统计力学给出连通 on-shell S-matrix 的能量导数；
现代两体重写明确写出 `S(E)=exp(2 i delta(E))`，从而

```math
\frac{1}{4\pi i}S^{-1}\overleftrightarrow{\partial_E}S
\longrightarrow \frac{1}{\pi}\partial_E\delta(E).
```

等价地，如果代码存的是完整 `arg(S)=2 delta`，测度必须是
`d(arg S)/(2 pi)`。这一点可以由 [Dashen--Ma--Bernstein, Phys. Rev. 187,
345](https://doi.org/10.1103/PhysRev.187.345) 和 [S-matrix formulation of
thermodynamics with N-body scatterings](https://link.springer.com/article/10.1140/epjc/s10052-017-5106-0)
逐式核对。

### 不能由来源唯一确定的部分

现有 PNJL/BU 论文把相位写成带 `i eta` 的传播子或逆传播子复比值/反正切，并
要求阈值束缚态跳变和高能相位消失。例如 [Dubinin et al., Phys. Rev. D 96,
094008](https://arxiv.org/html/1608.05383) 的 Eq. (14)--(17) 和 [Blaschke et al.,
Particles 3 (2020)](https://arxiv.org/html/1912.13162) 的 Eq. (14)--(20) 支持这一
物理对象和端点要求，但没有把 `-arg(Delta^R)`、主值 `unwrap` 的方向、有限窗口
`high_energy_zero` 平移规定为唯一数值算法。

因此下一 PR 将把三层分开：

1. `delta` 是物理散射相位，`S=exp(2 i delta)`；
2. `arg(S)` 是完整 S-matrix 相位，等于 `2 delta`；
3. `-arg(Delta^R_inverse)` 是当前项目的显式 propagator-to-phase diagnostic
   mapping，只有在独立 profile/端点测试通过后才可宣称与 `delta` 对齐。

## 2. Unequal-mass charged channel

[Dubinin et al., Phys. Rev. D 96, 094008](https://arxiv.org/html/1608.05383) 明确说明：
非等质量组分会在正电荷 kaon 通道产生 anomalous low-energy mode；Mott 转变时阈值
相位从 `pi` 跳到 `0`，而连续散射相位在高能端补偿。更直接的 K/pi 论文
[Blaschke et al., arXiv:1912.13162](https://arxiv.org/html/1912.13162) 给出：`K+`
和 `K-` 使用不同有序组分，Eq. (18)--(20) 分别定义标准 BU 密度、GBU 替换
`delta -> delta - sin(2 delta)/2`，并把 anomalous mode 的保留/去除作为独立比较。

这些来源支持“保留 ordered `Pi_us/Pi_su`、显式报告 anomalous component、分别做
Levinson/Mott gate”，但不支持有限网格简单符号变号计数作为唯一绑定态算法。后者
仍是本项目的 diagnostic gate，后续应与独立 pole/阈值计数对照。

## 3. `eta -> 0+`、PV 与端点

PNJL 论文使用 `q0+i eta` 或传播子在上下边界的组合来取得 cut 上的复相位，并把
阈值跳变、高能消失和连续谱补偿作为 Levinson/Mott 约束；[Wergieluk et al.,
arXiv:1212.5245](https://arxiv.org/html/1212.5245) 与上述 S08/S10 是本项目已有
的原始依据。它们支持以下实现分工：

- 实部：在 `eta -> 0+` 时取 Cauchy 主值；
- 虚部：由解析 cut 的 `i0+` 跳变给出，有限 `eta` 只作数值探针；
- 相位：先固定 `delta`/`arg(S)` 变量，再做连续分支和端点检查；
- 端点：阈值束缚态计数、Mott 前后相位差、高能尾部和积分测度必须一起验收。

来源没有规定项目所用 Gauss--Legendre 节点、`eta` 序列、`omega_max` 外推或误差阈值。
因此本项目会把 PV oracle、端点和收敛参数写入后续 PR 的测试契约，而不将任一有限
`eta` 结果视为解析极限。

## 4. 对下一步的影响

当前可以闭合并测试的是纯代数层：`delta <-> S`、`d(delta)/pi <->
d(arg S)/(2 pi)`、标量/对角矩阵的 `Im tr(S^-1 dS)`。当前不能提前宣称闭合的是
传播子到物理 S-matrix 的具体映射，以及真实 charged profile 的 PV、端点、Levinson、
Mott、节点和截断数值 gate。后者按用户批准的第三阶段单独实施，且不修改 production
默认、PNJLCore、`Omega_M` 或 `K03/K38` 路线。
