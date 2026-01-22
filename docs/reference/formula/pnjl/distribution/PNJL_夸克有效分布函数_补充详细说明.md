您提出的这个问题触及了PNJL模型的核心理论基础。让我从**第一性原理**出发，详细证明这个替换规则的合理性。

## 1. **理论基础：耦合到Polyakov loop的夸克传播子**

### 1.1 标准有限温度夸克传播子
在QCD的有限温度形式中，夸克的Matsubara传播子为：
\[
S(i\omega_n, \vec{p}) = \frac{1}{(i\omega_n + \mu)\gamma_0 - \vec{\gamma} \cdot \vec{p} - m}
\]
其中 \( i\omega_n = i(2n+1)\pi T \) 是费米子的Matsubara频率。

### 1.2 PNJL模型中的关键修改
在PNJL模型中，夸克耦合到**时间方向的规范背景场** \( A_4 = iA_0 \)：
\[
D_0 = \partial_0 - igA_4
\]
经过规范变换后，这等价于将夸克的Matsubara频率修改为：
\[
i\omega_n \rightarrow i\omega_n - gA_4
\]

由于Polyakov loop定义为 \( \Phi = \frac{1}{N_c} \mathrm{Tr} \exp(igA_4/T) \)，这个修改从根本上改变了夸克的统计分布。

## 2. **从配分函数严格推导**

### 2.1 夸克行列式的计算
PNJL模型的有效作用量来源于夸克的行列式：
\[
\mathcal{Z}_{\text{quark}} = \det[D\!\!\!\!/ + m + \mu\gamma_0]
\]
在静态Polyakov规范下，\( A_4 \) 是对角的，行列式可以分解为色空间：
\[
\det[D\!\!\!\!/ + m + \mu\gamma_0] = \prod_{c=1}^{N_c} \det[\partial_\tau - gA_4^c + \mu - \vec{\alpha} \cdot \vec{p} - \beta m]
\]

### 2.2 色迹与Matsubara求和
经过详细的Matsubara频率求和（这是证明的关键步骤），得到单夸克贡献的热力学势：
\[
\Omega_q = -2T \int \frac{d^3p}{(2\pi)^3} \mathrm{Tr}_c \left[ \ln(1 + L e^{-(E_p-\mu)/T}) + \ln(1 + L^\dagger e^{-(E_p+\mu)/T}) \right]
\]
其中 \( L = \exp(igA_4/T) \) 是Polyakov线。

### 2.3 分布函数的提取
从这个表达式可以**直接读出**修改后的分布函数：

**夸克分布**：
\[
f_q^+(E_p, \mu) = \frac{1}{N_c} \mathrm{Tr}_c \left[ \frac{1}{1 + L e^{-(E_p-\mu)/T}} \right]
\]

**反夸克分布**：
\[
f_q^-(E_p, \mu) = \frac{1}{N_c} \mathrm{Tr}_c \left[ \frac{1}{1 + L^\dagger e^{-(E_p+\mu)/T}} \right]
\]

## 3. **SU(3)情况下的显式形式**

对于 \( N_c = 3 \)，Polyakov矩阵 \( L \) 的特征值为 \( e^{i\phi_1}, e^{i\phi_2}, e^{i\phi_3} \)，且 \( \phi_1 + \phi_2 + \phi_3 = 0 \)。

### 3.1 色迹的显式计算
夸克分布变为：
\[
f_q^+(E, \mu) = \frac{1}{3} \sum_{c=1}^3 \frac{1}{1 + e^{i\phi_c} e^{-(E-\mu)/T}}
\]

经过代数运算，这正好给出：
\[
f_q^+(E, \mu) = \frac{\Phi e^{-(E-\mu)/T} + 2\bar{\Phi} e^{-2(E-\mu)/T} + e^{-3(E-\mu)/T}}{1 + 3\Phi e^{-(E-\mu)/T} + 3\bar{\Phi} e^{-2(E-\mu)/T} + e^{-3(E-\mu)/T}}
\]

其中 \( \Phi = \frac{1}{3}(e^{i\phi_1} + e^{i\phi_2} + e^{i\phi_3}) \)，\( \bar{\Phi} = \Phi^* \)。

## 4. **在B₀积分中的对应关系证明**

现在回到您图片中的B₀积分。原表达式中的各项来源于**Matsubara频率求和**：

### 4.1 标准的求和公式
对于费米子传播子，Matsubara求和给出：
\[
T\sum_n \frac{1}{i\omega_n - E} = f(E)
\]
其中 \( f(E) = 1/(e^{\beta E} + 1) \)。

### 4.2 PNJL模型中的修改
在PNJL模型中，由于 \( i\omega_n \rightarrow i\omega_n - gA_4 \)，求和变为：
\[
T\sum_n \frac{1}{i\omega_n - gA_4 - E} = \frac{1}{N_c} \mathrm{Tr}_c \left[ \frac{1}{1 + e^{\beta(E + gA_4)}} \right]
\]

这正是我们前面得到的修改后分布函数。

### 4.3 在B₀中的具体对应
在您图片的表达式(4.2)中：

- 第一项来源于极点 \( i\omega_n = -\mu_1 + E_1 \) 的留数，对应**粒子分布** \( f(E_1 - \mu_1) \)
- 在PNJL中，这变为 \( f_q^+(E_1, \mu_1) \)

- 第二项来源于极点 \( i\omega_n = -\mu_1 - E_1 \) 的留数，对应 \( f(-E_1-\mu_1) \)。

	这一项在最终的 B0 四项分解/实现中，常以“泡利阻塞因子”的形式出现：
	\[
	f(-E_1-\mu_1) \;\longleftrightarrow\; 1 - f^{-}(E_1,\mu_1)
	\]
	而在 PNJL 模型中同样可以利用恒等式把它写成“负能量的夸克分布”：
	\[
	1 - f_{\Phi}^{-}(E,\mu) = f_{\Phi}^{+}(-E,\mu)
	\]
	因此实现时经常统一写成 \( f_{\Phi}^{+}(pm\,E,\mu) \)（`pm=±1`），其中 `pm=-1` 对应 `f_{\Phi}^{+}(-E,\mu)`，也就是 `1-f_{\Phi}^{-}(E,\mu)`。

## 5. **物理图像的理解**

这个替换规则的深层物理原因是：

1. **色单态约束**：在禁闭相，只有色单态才能作为物理态传播
2. **统计权重修正**：Polyakov loop修正了夸克出现在热介质中的统计权重
3. **群积分结构**：分布函数中的复杂分母结构来源于对SU(3)群空间的积分

## 6. **数值验证**

这个替换规则已经被大量文献验证：
- Fukushima (2004) 的原始PNJL论文
- Ratti et al. (2006) 的QCD相图研究
- 众多后续的PNJL模型应用

**结论**：这个替换规则不是ad hoc的假设，而是从**第一性原理的配分函数计算**中严格导出的，是PNJL模型能够合理描述夸克禁闭效应的数学体现。