# 介子数密度：charged 化学势 profile

本文固定 charged / freeze-out 子链里，`meson chemical profile` 这一层的最小口径。

## 1. 角色

当前路径分层已经拆成三层：

1. freeze-out profile  
   负责 `\sqrt{s_{NN}} \to (T,\mu_B)`
2. flavor chemical profile  
   负责 `\mu_q \to (\mu_u,\mu_d,\mu_s)`
3. meson chemical profile  
   负责：
   - `\mu_\pi`
   - `\mu_K`
   - `d_\pi`
   - `d_K`
   - charge label

因此当前 workflow 中，charged / freeze-out 数密度扫描的最小结构是：

```text
sqrt(s_NN)
-> freeze-out profile
-> flavor chemical profile
-> meson chemical profile
-> meson-density workflow
```

## 2. 当前实现的物理边界

当前 meson chemical profile 是**介子层**对象，而不是 flavor-level 化学势对象。

也就是说，当前实现直接固定的是：

```math
\mu_\pi,\quad \mu_K,\quad d_\pi,\quad d_K
```

而不是直接固定：

```math
\mu_u,\quad \mu_d,\quad \mu_s.
```

这样做的原因是：

1. 当前数密度后处理层本来就直接消费 `\mu_\pi / \mu_K`；
2. 若现在直接把 `\mu_s` 注入平衡态主链，会把问题扩成 flavor-`\mu` solver 子链；
3. 对 charged / freeze-out 最小工作流而言，先把 meson-level 口径显式化更稳。

## 3. 当前配置口径

当前已配置：

1. `default`
   - 总 `\pi / K` 聚合通道
   - `d_\pi = 3`
   - `d_K = 4`
   - `\mu_\pi = 0`
   - `\mu_K = 0`
2. `charge_resolved_neutral`
   - charge-resolved 最小口径
   - `d_\pi = d_K = 1`
   - `\mu_\pi = \mu_K = 0`
3. `blaschke2019_mu_pi_100`
   - 面向 `Blaschke:2019col` Figure 4 right panel 的最小对照 profile
   - `d_\pi = d_K = 1`
   - `\mu_\pi = 100~\mathrm{MeV}`
   - `\mu_K = 0`
4. `blaschke2019_mu_pi_134p5`
   - 同上
   - `\mu_\pi = 134.5~\mathrm{MeV}`

## 4. 与 `\mu_s` 的边界

当前还不能把 `\mu_s` 直接等同于 `\mu_K`，原因是：

1. `\mu_s` 是 flavor-level strange chemical potential；
2. `\mu_K` 是介子层化学势；
3. 两者的关系依赖：
   - flavor 约束
   - freeze-out / critical path 规则
   - 文献中采用的近似

因此当前仓库把这件事拆成三层：

1. freeze-out profile
2. flavor chemical profile
3. meson chemical profile

也就是说，`mu_s` 现已升格为显式 profile 对象，但仍与 `mu_pi / mu_K` 分层治理。

## 5. 与文献对照的关系

当前这一层足够支持：

1. charge-resolved `K^-/\pi^-`
2. 固定 `\mu_\pi`
3. freeze-out baseline 路径

但还**不足以**严格复现文献中所有 horn / stitched-path 结果，因为它还没有包含：

1. `\mu_s(x)`
2. pseudo-critical / stitched path
3. anomaly mode 或更高层现象学规则

因此当前这层的角色是：

- 让 charged / freeze-out 子链先进入统一 workflow
- 不过早把 flavor-level 现象学和路径代理混进底层公式层

### 5.1 与旧 `K / sigma_K` 标签的边界

需要特别注意，charged meson-density 子链里的：

- `K_plus`
- `K_minus`

是后加入的 charged / ordered strange-light label；它们不同于更早主链里的 generic：

- `K`
- `sigma_K`

当前代码语义下：

- `K_plus` 走 ordered `(u,s)` 赝标量分支
- `K_minus` 走 ordered `(s,u)` 赝标量分支
- 旧的 generic `K` 在已支持 charged label 的实现里，通常与 `(u,s)` 分支同向
- 旧的 generic `sigma_K` 仍不是 charge-resolved scalar label

因此，charged meson-density workflow 不应把旧 `sigma_K` 自动解释成某个固定电荷态；若需要这一步，必须在更低层传播子 / 极化函数语义上显式拆分。
