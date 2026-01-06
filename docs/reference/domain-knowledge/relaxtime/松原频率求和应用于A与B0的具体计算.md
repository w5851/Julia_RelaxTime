松原频率求和是有限温度场论中处理离散虚数频率的常用技术。表达式  
\[
A(m, \mu, \beta, \Lambda) = 16\pi^2 \int \frac{d^3p}{(2\pi)^3} \left[ f(E - \mu) \frac{1}{2E} - f(-E - \mu) \frac{1}{2E} \right]
\]  
来源于对一个包含费米子松原频率的求和进行围道积分，具体步骤如下。

### 1. 松原频率求和的一般形式
考虑费米子的松原频率 \(\omega_n = (2n+1)\pi T\)，其中 \(T = 1/\beta\)。常见的一类求和形式为  
\[
T \sum_{n=-\infty}^{\infty} \frac{1}{(i\omega_n + \mu)^2 - E^2},
\]  
其中 \(E = \sqrt{\mathbf{p}^2 + m^2}\)。这种形式出现在自由费米子传播子或单圈图的动量积分中。

### 2. 部分分式分解
对求和项进行分解：  
\[
\frac{1}{(i\omega_n + \mu)^2 - E^2} = \frac{1}{2E} \left( \frac{1}{i\omega_n + \mu - E} - \frac{1}{i\omega_n + \mu + E} \right).
\]  
于是，  
\[
T \sum_{n} \frac{1}{(i\omega_n + \mu)^2 - E^2} = \frac{1}{2E} \left[ T \sum_{n} \frac{1}{i\omega_n + \mu - E} - T \sum_{n} \frac{1}{i\omega_n + \mu + E} \right].
\]

### 3. 应用标准求和公式
对于费米子，标准公式为  
\[
T \sum_{n} \frac{1}{i\omega_n - \xi} = f(\xi) = \frac{1}{e^{\beta \xi} + 1}.
\]  
将公式应用于两项：
- 第一项：令 \(\xi = E - \mu\)，则  
  \[
  T \sum_{n} \frac{1}{i\omega_n + \mu - E} = T \sum_{n} \frac{1}{i\omega_n - (E - \mu)} = f(E - \mu).
  \]
- 第二项：令 \(\xi = -E - \mu\)，则  
  \[
  T \sum_{n} \frac{1}{i\omega_n + \mu + E} = T \sum_{n} \frac{1}{i\omega_n - (-E - \mu)} = f(-E - \mu).
  \]

因此，  
\[
T \sum_{n} \frac{1}{(i\omega_n + \mu)^2 - E^2} = \frac{1}{2E} \left[ f(E - \mu) - f(-E - \mu) \right].
\]

### 4. 物理意义
- \(f(E - \mu)\) 是夸克（粒子）的费米-狄拉克分布，化学势为 \(+\mu\)。
- \(f(-E - \mu)\) 与反夸克分布相关：利用 \(f(-x) = 1 - f(x)\)，有  
  \[
  f(-E - \mu) = 1 - f(E + \mu),
  \]  
  其中 \(f(E + \mu)\) 是反夸克的分布（化学势为 \(-\mu\)）。  
  因此，  
  \[
  f(E - \mu) - f(-E - \mu) = f(E - \mu) + f(E + \mu) - 1,
  \]  
  这一组合常出现在热力学量的计算中，如净粒子数密度以外的其他量（如自能、有效势等）。

### 5. 与原文表达式的联系
原文中 \(A\) 的定义可能对应于某个物理量（如自由能密度的一部分），其动量积分包含上述松原求和的结果，因此直接导出  
\[
A \propto \int \frac{d^3p}{(2\pi)^3} \left[ f(E - \mu) \frac{1}{2E} - f(-E - \mu) \frac{1}{2E} \right].
\]

这种形式是松原求和的标准结果，反映了粒子与反粒子贡献的对称性。在更复杂的模型（如PNJL）中，分布函数会修正，但自由费米气体下的推导仍为理解其起源提供了基础。
---
对于给定的表达式：

\[
B_0(k, m_1, \mu_1, m_2, \mu_2, i\nu_m, \beta, \Lambda) = \frac{16\pi^2}{\beta} \sum_n \exp(i\omega_n \eta) \int \frac{d^3p}{(2\pi)^3} \frac{1}{((i\omega_n + \mu_1)^2 - E_1^2)} \frac{1}{((i\omega_n - i\nu_m + \mu_2)^2 - E_2^2)},
\]

其中 \( E_1 = \sqrt{p^2 + m_1^2} \), \( E_2 = \sqrt{(\mathbf{p}-\mathbf{k})^2 + m_2^2} \)，\(\omega_n = (2n+1)\pi T\) 为费米子松原频率，\(\nu_m = 2m\pi T\) 为玻色子松原频率，\(\eta \to 0^+\) 为收敛因子，处理步骤如下：

### 1. 部分分式分解
将每个传播子分母分解为单极点形式：
\[
\frac{1}{(i\omega_n + \mu_1)^2 - E_1^2} = \frac{1}{2E_1} \left( \frac{1}{i\omega_n + \mu_1 - E_1} - \frac{1}{i\omega_n + \mu_1 + E_1} \right),
\]
\[
\frac{1}{(i\omega_n - i\nu_m + \mu_2)^2 - E_2^2} = \frac{1}{2E_2} \left( \frac{1}{i\omega_n - i\nu_m + \mu_2 - E_2} - \frac{1}{i\omega_n - i\nu_m + \mu_2 + E_2} \right).
\]

乘积展开为四个项的求和：
\[
\frac{1}{((i\omega_n + \mu_1)^2 - E_1^2)((i\omega_n - i\nu_m + \mu_2)^2 - E_2^2)} = \frac{1}{4E_1 E_2} \sum_{s_1=\pm} \sum_{s_2=\pm} \sigma(s_1, s_2) \frac{1}{(i\omega_n + \mu_1 - s_1 E_1)(i\omega_n - i\nu_m + \mu_2 - s_2 E_2)},
\]
其中符号因子 \(\sigma(s_1, s_2) = s_1 s_2\)（即 \(s_1, s_2\) 分别取 ±1，乘积决定正负）。

### 2. 对每个项进行松原频率求和
对每个 \((s_1, s_2)\) 组合，将分式拆分为单极点之差：
\[
\frac{1}{(i\omega_n + \mu_1 - s_1 E_1)(i\omega_n - i\nu_m + \mu_2 - s_2 E_2)} = \frac{1}{D(s_1, s_2)} \left( \frac{1}{i\omega_n - i\nu_m + \mu_2 - s_2 E_2} - \frac{1}{i\omega_n + \mu_1 - s_1 E_1} \right),
\]
其中
\[
D(s_1, s_2) = (s_2 E_2 - s_1 E_1) + (\mu_2 - \mu_1 - i\nu_m).
\]

利用费米子松原频率求和公式：
\[
T \sum_n \frac{e^{i\omega_n \eta}}{i\omega_n - \xi} = -f(\xi),
\]
其中 \(f(\xi) = 1/(e^{\beta \xi} + 1)\) 为费米分布函数。可得：
\[
T \sum_n \frac{e^{i\omega_n \eta}}{i\omega_n + \mu_1 - s_1 E_1} = -f(s_1 E_1 - \mu_1),
\]
\[
T \sum_n \frac{e^{i\omega_n \eta}}{i\omega_n - i\nu_m + \mu_2 - s_2 E_2} = -f(i\nu_m - \mu_2 + s_2 E_2).
\]

因此，每个 \((s_1, s_2)\) 项的求和结果为：
\[
T \sum_n \frac{e^{i\omega_n \eta}}{(i\omega_n + \mu_1 - s_1 E_1)(i\omega_n - i\nu_m + \mu_2 - s_2 E_2)} = \frac{1}{D(s_1, s_2)} \left( -f(i\nu_m - \mu_2 + s_2 E_2) + f(s_1 E_1 - \mu_1) \right).
\]

### 3. 组合并积分
将所有项组合，并代入原表达式：
\[
B_0 = 16\pi^2 \int \frac{d^3p}{(2\pi)^3} \frac{1}{4E_1 E_2} \sum_{s_1, s_2} s_1 s_2 \frac{1}{D(s_1, s_2)} \left( f(s_1 E_1 - \mu_1) - f(i\nu_m - \mu_2 + s_2 E_2) \right),
\]
其中已取 \(\eta \to 0^+\) 忽略指数因子。

### 4. 说明
- 该结果为松原频率 \(i\nu_m\) 的函数，实际计算中常通过解析延拓 \(i\nu_m \to \omega + i\epsilon\) 得到物理关联函数。
- 分布函数中的变量可能为复数，但通常通过关系 \(f(z) = 1/(e^{\beta z} + 1)\) 定义解析延拓。
- 此表达式可用于计算两极子图的热力学量或响应函数。