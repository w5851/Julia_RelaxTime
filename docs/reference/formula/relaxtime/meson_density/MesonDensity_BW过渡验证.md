# 介子数密度：BW 过渡验证

本文固定当前介子数密度主线中的 BW 过渡验证层，用于在稳定粒子极限和 BU 全相移公式之间建立可检查的中间台阶。

相关文档：

- [MesonDensity_稳定粒子与KPi比值.md](./MesonDensity_%E7%A8%B3%E5%AE%9A%E7%B2%92%E5%AD%90%E4%B8%8EKPi%E6%AF%94%E5%80%BC.md)
- [MesonDensity_BU相移公式.md](./MesonDensity_BU%E7%9B%B8%E7%A7%BB%E5%85%AC%E5%BC%8F.md)
- [MesonMass_RPA_Pole.md](../propagator/MesonMass_RPA_Pole.md)

## 1. BW 极点近似

文献中对窄宽度共振采用：

```math
z_M = E_M - i\Gamma_M/2.
```

对应的相移写为：

```math
\delta_M(\omega)
= -\arctan\frac{\Gamma_M/2}{\omega-E_M}.
```

### 来源

- `Blaschke:2020bzh` Eq. (26) 周边

## 2. 相移导数

由上式可得：

```math
\frac{d\delta_M}{d\omega}
= \frac{1}{2}\frac{\Gamma_M}{(\omega-E_M)^2+\Gamma_M^2/4}.
```

### 来源

- `Blaschke:2020bzh` Eq. (26)

## 3. BW 数密度公式

把 BW 相移导数代入 BU 原式，可得：

```math
n_M(T)
= d_M \int \frac{d^3q}{(2\pi)^3}
\int \frac{d\omega}{2\pi}
g_M(\omega)\frac{\Gamma_M}{(\omega-E_M)^2+\Gamma_M^2/4}.
```

### 来源

- `Blaschke:2020bzh` Eq. (27)

## 4. 与当前项目的对接方式

当前项目不直接把文献中的 BW 参数当作独立输入。更合适的做法是：

1. 先从当前项目主极点条件
   ```math
   1-4K_M\Pi_M(p_0,q;\xi)=0
   ```
   求得 `M` 与 `\Gamma`
2. 再把它们代回 BW 近似表达式
3. 用作稳定极限与 BU 之间的过渡验证层

这样做可以保证 BW 近似仍服从当前仓库的传播子和极点口径。

## 5. 稳定极限退化关系

当

```math
\Gamma_M \to 0
```

时，BW 数密度应退化回稳定粒子极限。

这条退化关系是当前 `Phase D` 的关键一致性检查之一：

1. 若退化不成立，则实现、单位或归一化存在问题；
2. 若退化成立，则 BW 层可作为过渡验证台阶继续保留。

## 6. 当前用途

本文只服务于：

1. `Phase D3`：是否需要 BW 过渡验证的判断
2. `Phase D4`：把 BW 明确标记为“验证台阶而非最终目标”
3. 后续与 BU 结果的最小差异对照

本文不把 BW 视为当前主线最终物理结果。
