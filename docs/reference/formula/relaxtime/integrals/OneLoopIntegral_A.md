# 有限温度有限密度场论中的单圈积分 A（one-line integral）- 实现文档

## 1. 公式标识
- **程序实现**: `A`
- **对应文件**: [src/relaxtime/OneLoopIntegrals.jl](/C:/Users/Wmzx/.codex/worktrees/346c/Julia_RelaxTime/src/relaxtime/OneLoopIntegrals.jl:536)

## 2. 数学表达式

### 2.1 Matsubara 求和前
```math
A(m, \mu, \beta, \Lambda)
= \frac{16\pi^2}{\beta}\sum_n e^{i\omega_n\eta}
\int \frac{d^3p}{(2\pi)^3}
\frac{1}{(i\omega_n+\mu)^2-E^2},
\qquad E=\sqrt{p^2+m^2}.
```

这里项目实现最终会把角积分做掉，仅保留 `p\in[0,\Lambda]` 的一维积分。

先对被积函数做部分分式分解：
```math
\frac{1}{(i\omega_n+\mu)^2-E^2}
= \frac{1}{2E}
\left[
\frac{1}{i\omega_n-(E-\mu)}
-\frac{1}{i\omega_n-(-E-\mu)}
\right].
```

在当前文档采用的收敛因子 `e^{i\omega_n\eta}`、`\eta\to0^+` 约定下，有
```math
\frac{1}{\beta}\sum_n \frac{e^{i\omega_n\eta}}{i\omega_n-a}
= n_F(a)-1,
\qquad
n_F(x)=\frac{1}{e^{\beta x}+1}.
```

因此
```math
\begin{aligned}
\frac{1}{\beta}\sum_n e^{i\omega_n\eta}
\frac{1}{(i\omega_n+\mu)^2-E^2}
&=
\frac{1}{2E}
\Bigl[
n_F(E-\mu)-1
-\bigl(n_F(-E-\mu)-1\bigr)
\Bigr] \\
&=
-\frac{1}{2E}
\Bigl[
1-n_F(E-\mu)-n_F(E+\mu)
\Bigr].
\end{aligned}
```

这里最后一步用了
```math
n_F(-x)=1-n_F(x).
```

### 2.2 Matsubara 求和后常见 shorthand
```math
A(m,\mu,\beta,\Lambda)
= -16\pi^2 \int \frac{d^3p}{(2\pi)^3}\frac{1}{2E}
\left[1-f(E-\mu)-f(E+\mu)\right]
```
即
```math
A(m,\mu,\beta,\Lambda)
= -4\int_0^\Lambda dp\,\frac{p^2}{E}
\left[1-f(E-\mu)-f(E+\mu)\right].
```

对普通 NJL 记号，上式里的 `f(E-\mu)` 与 `f(E+\mu)` 就是正能量支上的粒子/反粒子占据数。进入 PNJL 后，它们要替换成背景 Polyakov 环下的有效分布。

对当前项目的 PNJL 记号，应明确理解为：
```math
f(E-\mu)\;\widehat{=}\;n_q(E;\mu,\Phi,\bar\Phi),
\qquad
f(E+\mu)\;\widehat{=}\;n_{\bar q}(E;\mu,\Phi,\bar\Phi),
```
而不是 `n_q(-E;\mu,\Phi,\bar\Phi)` 或其他负能量支写法。负能量支恒等式只在某些 Matsubara 求和重排或 B0 实现改写中作为技术手段出现，不应反向替代这里 `2.2` 的正能量占据数定义。

## 3. 与当前 PNJL 代码实现的精确对应

### 3.1 当前代码中的两类 PNJL 分布

项目当前实现见 [src/QuarkDistribution.jl](/C:/Users/Wmzx/.codex/worktrees/346c/Julia_RelaxTime/src/QuarkDistribution.jl:50)：

```math
n_q(E;\mu,\Phi,\bar\Phi)=
\frac{\Phi e^{-\beta(E-\mu)}+2\bar\Phi e^{-2\beta(E-\mu)}+e^{-3\beta(E-\mu)}}
{1+3\Phi e^{-\beta(E-\mu)}+3\bar\Phi e^{-2\beta(E-\mu)}+e^{-3\beta(E-\mu)}},
```

```math
n_{\bar q}(E;\mu,\Phi,\bar\Phi)=
\frac{\bar\Phi e^{-\beta(E+\mu)}+2\Phi e^{-2\beta(E+\mu)}+e^{-3\beta(E+\mu)}}
{1+3\bar\Phi e^{-\beta(E+\mu)}+3\Phi e^{-2\beta(E+\mu)}+e^{-3\beta(E+\mu)}}.
```

对应到代码入口分别是：
- `n_q` -> [quark_distribution(E, μ, T, Φ, Φbar)](/C:/Users/Wmzx/.codex/worktrees/346c/Julia_RelaxTime/src/QuarkDistribution.jl:50)
- `n_{\bar q}` -> [antiquark_distribution(E, μ, T, Φ, Φbar)](/C:/Users/Wmzx/.codex/worktrees/346c/Julia_RelaxTime/src/QuarkDistribution.jl:56)

### 3.2 PNJL 下的完整展开

在 PNJL 中，可先把背景时间规范场写成颜色依赖的虚化学势移位；对每个颜色本征值 `\varphi_c`，有
```math
f_c^+(E)=\frac{1}{e^{\beta(E-\mu-i\varphi_c)}+1},
\qquad
f_c^-(E)=\frac{1}{e^{\beta(E+\mu+i\varphi_c)}+1}.
```

对颜色求迹并用 `\Phi,\bar\Phi` 表示后，就得到当前项目实现的两类 PNJL 有效分布：
```math
n_q(E;\mu,\Phi,\bar\Phi)=
\frac{\Phi e^{-\beta(E-\mu)}+2\bar\Phi e^{-2\beta(E-\mu)}+e^{-3\beta(E-\mu)}}
{1+3\Phi e^{-\beta(E-\mu)}+3\bar\Phi e^{-2\beta(E-\mu)}+e^{-3\beta(E-\mu)}},
```
```math
n_{\bar q}(E;\mu,\Phi,\bar\Phi)=
\frac{\bar\Phi e^{-\beta(E+\mu)}+2\Phi e^{-2\beta(E+\mu)}+e^{-3\beta(E+\mu)}}
{1+3\bar\Phi e^{-\beta(E+\mu)}+3\Phi e^{-2\beta(E+\mu)}+e^{-3\beta(E+\mu)}}.
```

因此，`2.2` 在本项目中的正确完整展开式是
```math
A_f
= -16\pi^2 \int \frac{d^3p}{(2\pi)^3}\frac{1}{2E_f}
\left[
1-n_q(E_f;\mu_f,\Phi,\bar\Phi)-n_{\bar q}(E_f;\mu_f,\Phi,\bar\Phi)
\right].
```

把角积分做掉后，
```math
A_f
= -4\int_0^\Lambda dp\,\frac{p^2}{E_f}
\left[
1-n_q(E_f;\mu_f,\Phi,\bar\Phi)-n_{\bar q}(E_f;\mu_f,\Phi,\bar\Phi)
\right].
```

### 3.3 物理展开形式与代码实现形式

当前代码实际采用的是完全等价的重排形式
```math
A_f
= 4\int_0^\Lambda dp\,\frac{p^2}{E_f}
\left[n_q(E_f;\mu_f,\Phi,\bar\Phi)+n_{\bar q}(E_f;\mu_f,\Phi,\bar\Phi)-1\right].
```

因为括号满足
```math
n_q+n_{\bar q}-1
=-\left(1-n_q-n_{\bar q}\right).
```

因此要特别注意：
- `- [1-n_q-n_{\bar q}]` 与 `n_q+n_{\bar q}-1` 才是完全等价的；
- 当前代码 [A(...)](/C:/Users/Wmzx/.codex/worktrees/346c/Julia_RelaxTime/src/relaxtime/OneLoopIntegrals.jl:536) 的实现路径是“先减常数项积分，再加热分布项积分”，因此与上式严格一致；
- 之前若把 `A` 直接写成 `+ [1-n_q-n_{\bar q}]`，那就是符号写反了。

代码级写法可直接从实现读出：
```julia
integral = -const_integral_term_A(m)
dist_quark = distribution_value(:pnjl, :plus, E, μ, T, Φ, Φbar)
dist_antiquark = distribution_value(:pnjl, :minus, E, μ, T, Φ, Φbar)
integral += weight_p * p^2 / E * (dist_quark + dist_antiquark)
return 4.0 * integral
```

这里：
- `distribution_value(:pnjl, :plus, ...)` 返回 `n_q`
- `distribution_value(:pnjl, :minus, ...)` 返回 `n_{\bar q}`

## 4. 与夸克凝聚 `\phi_f=\langle\bar q_f q_f\rangle` 的关系

若采用项目主线物理口径
```math
\phi_f \equiv \langle \bar q_f q_f \rangle,
```
则标准有限温度凝聚写法为
```math
\phi_f
= -\frac{N_c}{\pi^2}m_f\int_0^\Lambda dp\,\frac{p^2}{E_f}
\left[1-n_q(E_f;\mu_f,\Phi,\bar\Phi)-n_{\bar q}(E_f;\mu_f,\Phi,\bar\Phi)\right].
```

与当前代码实现的 `A_f` 对照，可写成
```math
A_f
= -4\int_0^\Lambda dp\,\frac{p^2}{E_f}
\left[1-n_q(E_f;\mu_f,\Phi,\bar\Phi)-n_{\bar q}(E_f;\mu_f,\Phi,\bar\Phi)\right],
```
因此两者满足
```math
\phi_f=\frac{N_c}{4\pi^2}m_fA_f.
```

这说明：
- 在当前实现口径下，`A_f` 与物理凝聚 `\phi_f` 直接同号对应；
- 真空下两者都应为负；
- 后续若某处代码引入 `-\phi_f`，那是历史 helper 变量口径，不应反向污染这里的物理定义。

## 5. 常数项积分的解析形式

对代码里的常数项，需要用到
```math
\int \frac{x^2}{\sqrt{x^2+m^2}}dx
=\frac{x}{2}\sqrt{x^2+m^2}
-\frac{m^2}{2}\ln\left|x+\sqrt{x^2+m^2}\right|+C.
```

因此
```math
\int_0^\Lambda dp\,\frac{p^2}{\sqrt{p^2+m^2}}
```
可以解析计算，这也是当前实现把常数项单独拿出来处理的原因。
