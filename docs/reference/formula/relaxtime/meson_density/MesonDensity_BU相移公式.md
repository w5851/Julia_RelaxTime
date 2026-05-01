# 介子数密度：BU 相移公式

本文固定当前介子数密度主线中的 BU 相移数密度表达，以及它与当前项目传播子定义的对接方式。

相关文档：

- [2026-05-01_介子数密度BU主线公式草案.md](../../../../dev/archived/2026-05-01_%E4%BB%8B%E5%AD%90%E6%95%B0%E5%AF%86%E5%BA%A6BU%E4%B8%BB%E7%BA%BF%E5%85%AC%E5%BC%8F%E8%8D%89%E6%A1%88.md)
- [2026-05-01_PhaseB_介子数密度记号对照与最小物理口径.md](../../../../dev/archived/2026-05-01_PhaseB_%E4%BB%8B%E5%AD%90%E6%95%B0%E5%AF%86%E5%BA%A6%E8%AE%B0%E5%8F%B7%E5%AF%B9%E7%85%A7%E4%B8%8E%E6%9C%80%E5%B0%8F%E7%89%A9%E7%90%86%E5%8F%A3%E5%BE%84.md)
- [Propagator_传播子byPolarization.md](../propagator/Propagator_%E4%BC%A0%E6%92%AD%E5%AD%90byPolarization.md)

## 1. BU 主公式

广义 Beth-Uhlenbeck 形式写为：

```math
n_M(T)
= d_M \int \frac{d^3q}{(2\pi)^3}
\int \frac{d\omega}{2\pi}
g_M(\omega)\frac{d\delta_M(\omega,q)}{d\omega}.
```

### 来源

- `Blaschke:2020bzh` Eq. (24)
- `Blaschke:2021yml` Eq. (8)

## 2. 分部积分后的等价形式

对上式分部积分后，得到：

```math
n_M(T)
= \frac{d_M}{T}
\int \frac{dq\,q^2}{2\pi^2}
\int_0^\infty \frac{d\omega}{2\pi}
g_M(\omega)\bigl[1+g_M(\omega)\bigr]\delta_M(\omega).
```

### 来源

- `Blaschke:2020bzh` Eq. (25)

## 3. 当前优先实现口径

当前主线优先采用**分部积分后的 `\delta_M` 本体形式**，原因是：

1. 当前项目后续相移将由传播子数值相位提取；
2. 直接对 `\delta_M` 求 `\omega` 导数更容易放大数值噪声；
3. 相移本体形式更利于先检查 Levinson 约束和相位分支选择。

导数形式保留为理论原式和交叉验证口径。

## 4. 当前项目中的相移定义

文献常通过传播子的极坐标形式定义相移：

```math
S^M_{ij}(\omega,q)=|S^M_{ij}(\omega,q)|e^{i\delta_M(\omega,q)}.
```

对当前项目简单介子道，优先采用项目化解释：

```math
\delta_M(\omega,q;\xi)=\arg \mathcal{D}_M(\omega+i\eta,q;\xi).
```

其中当前项目主传播子为：

```math
\mathcal{D}_M(\omega,q;\xi)
= \frac{2K_M}{1-4K_M\Pi_M(\omega,q;\xi)}.
```

这意味着 BU 相移主线必须服从当前项目的传播子与极点条件口径，而不是直接切换成文献里另一套归一化记号。

## 5. Levinson 约束与物理检查

在 BU 实现中，应满足下列物理检查：

1. 束缚态对应相移 `+\pi` 跳变；
2. 高能端相移应回到 `0`；
3. Mott 解离后，原束缚态跳变应转为连续谱/共振结构。

这些检查先于数值积分结果本身，是相移实现是否可信的基本门禁。

## 6. 当前边界

本文只固定：

1. BU 数密度表达；
2. 相移与当前项目传播子之间的关系；
3. 分部积分形式作为当前优先实现口径。

它不负责：

1. BW 近似；
2. 稳定粒子极限；
3. 冻结线现象学拟合；
4. 各向异性下的完整角向积分实现细节。
