# 介子数密度：flavor 化学势 profile

本文固定 charged / freeze-out 子链里，flavor-level 化学势现象学的最小对象。

## 1. 当前定位

当前链路已拆成三层：

```text
sqrt(s_NN)
-> freeze-out profile
-> flavor chemical profile
-> meson chemical profile
-> meson-density workflow
```

其中：

1. freeze-out profile  
   负责 `sqrt(s_NN) -> (T, mu_B)`
2. flavor chemical profile  
   负责 `mu_q -> (mu_u, mu_d, mu_s)`
3. meson chemical profile  
   负责 `mu_pi, mu_K, d_pi, d_K`

## 2. 为什么这一层要单独存在

`mu_s` 不是 `mu_K`。

更准确地说：

1. `mu_s` 属于 flavor-level quark chemical potential
2. `mu_K` 属于 meson-level effective chemical potential
3. 两者不应继续混在同一个 profile 对象里治理

因此当前最小实现把它显式拆出来，避免后续把 strange phenomenology 继续塞进文档注释或脚本常量。

## 3. 当前最小规则族

当前先固定一族最保守、最容易治理的规则：

```math
\mu_u = c_u \mu_q,\quad
\mu_d = c_d \mu_q,\quad
\mu_s = c_s \mu_q,
```

其中：

```math
\mu_q = \mu_B / 3.
```

这对应当前实现中的：

- `family = flavor_ratio_to_muq_v1`

## 4. 当前已配置 profile

### 4.1 `default`

```math
c_u = c_d = c_s = 1.
```

这等价于当前项目原本的对称 `FixedMu` 口径：

```math
\mu_u = \mu_d = \mu_s = \mu_q.
```

### 4.2 `blaschke2017_mu_s_0p2`

```math
c_u = 1,\quad
c_d = 1,\quad
c_s = 0.2.
```

即：

```math
\mu_s = 0.2 \mu_q.
```

当前把它作为 strange phenomenology 的第一条显式 profile，用来承接文献中常见的最小 `mu_s` 代理口径。

## 5. 当前工程判断

截至当前版本：

1. freeze-out 不升格为底层 `ConstraintMode`
2. `mu_s` 先升格为 path/workflow 层的显式 profile
3. profile 结果直接透传到平衡态主链的 flavor `mu` 求解

因此当前推荐理解是：

- `freeze-out` 是路径层对象
- `flavor chemical profile` 是路径上的 flavor 规则
- 底层仍复用已有 equilibrium / meson workflow，而不是另造一套脚本拼链
