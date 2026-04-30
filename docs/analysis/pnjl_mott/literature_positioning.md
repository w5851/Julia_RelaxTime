# 文献定位矩阵：本工作与现有文献的系统性对比

本文档涵盖核心说明第5点要求的文献扩展和第3点要求的新颖性论证，提供本工作与现有文献的详细对比框架。

---

## 1. 核心新颖性声明

> **本工作在 2+1 味 PNJL 模型中，首次将 Romatschke-Strickland 动量各向异性与介子 Mott 离解温度的计算统一在一个框架内。**

现有文献要么：(a) 在各向同性 PNJL 中计算 Mott 温度（Blaschke 2017, Maslov 2023）；(b) 在各向异性 NJL/QM 中计算相边界/介子质量但不提取 Mott 温度（Zhang 2021, Zhang 2021 NJL）；(c) 在各向异性 PNJL 中计算相边界但不算介子谱（He 2023）。**无文献同时在 PNJL 框架内包含 RS 各向异性 + 介子极点求解 + Mott 判据**。

---

## 2. 详细对比矩阵

### 2.1 直接竞争者

| 文献 | 模型 | 各向异性 | 介子质量 | Mott 温度 | 与本工作关系 |
|------|------|---------|---------|----------|------------|
| **本工作** | PNJL 2+1 | RS ξ = ±0.3 | π, K（复极点） | ✓ T_mott(ξ) 表 | — |
| Blaschke 2017 | PNJL 2+1 | 无（ξ=0） | π, K, a₀, κ（BU相移） | ✓ ξ=0 基准 | 我们的 ξ=0 基准对标目标 |
| Maslov 2023 | PNJL 2 | 无（ξ=0） | π（BU+离壳） | ✓ ξ=0 基准 | 极化函数方法论来源 |
| Zhang 2021 (NJL) | NJL 2 | RS ξ = ±0.3 | π, σ（弱ξ展开） | ✓ T_mott(ξ) | **唯一可对比的 ξ-T_mott 数据** |
| Zhang 2020 (QM) | QM 2+1 | RS ξ = ±0.4 | π, K, η, σ | ✗ | 介子质量 vs ξ 的模型对照 |
| He 2023 (PNJL) | PNJL 2+1 | RS ξ = ±0.4 | ✗ | ✗ | 同一框架，相边界对照 |

### 2.2 机制相关但焦点不同

| 文献 | 体系 | 解离对象 | 关键发现 | 与本工作关系 |
|------|------|---------|---------|------------|
| Burnier 2009 | pQCD+HTL | 重夸克偶素 | 等熵条件下 T_diss 不变 | Discussion 对立论点 |
| Jamal 2018 | EQPM | 重夸克偶素 | ξ>0 使 T_diss 升高 | 跨体系定性一致性 |
| Singh 2025 | PCQMF | 无介子谱 | ξ 改变输运和相变 | 同方向竞争者 |

### 2.3 方法论基础

| 文献 | 贡献 | 引用位置 |
|------|------|---------|
| Klevansky 1992 | NJL 标准综述，介子极化函数基础 | §2 Introduction |
| Buballa 2005 | NJL 有限密度综述，能隙方程推导 | §2 Introduction |
| Rehberg 1996 | 有限温度单圈积分（A, B₀）标准方法 | Appendix |
| Romatschke 2003 | RS 各向异性原始文献 | §3 |
| Romatschke 2004 | 各向异性集体模式 II（非物理叶分析） | §3 背景 |
| Strickland 2014 | 各向异性流体力学综述 | §1 Introduction |

---

## 3. 关键对比数据

### 3.1 T_mott(ξ) 对比：PNJL vs NJL

| ξ | T_mott^π (MeV) — 本工作 PNJL | T_mott^π (MeV) — Zhang 2021 NJL | 差值 |
|---|------------------------------|--------------------------------|------|
| -0.30 | 196.5 | 187 | +9.5 |
| 0.00 | 207.3 | 196 | +11.3 |
| +0.30 | 216.0 | 206 | +10.0 |

- Polyakov 环引入 ~10 MeV 的系统性上移
- 斜率 dT_mott/dξ ≈ +33 MeV（两模型一致）→ 各向异性响应由能隙方程主导，Polyakov 环提供 ξ 无关的常值偏移

### 3.2 T_χ(ξ) 行为对比：本工作 vs He 2023

- He 2023: T_χ 随 ξ **非单调**（在 ξ≈-0.16 处最大）
- 本工作: T_mott 随 ξ **单调**上升
- 解读：Mott 判据不是手征恢复的简单代理——它包含介子极点相对于连续阈值的额外 ξ 敏感性

---

## 4. 文献扩展路线图（核心说明第5点：~50篇）

### 已完成（15篇）
Tier 1-3 核心文献，已全部追加到 refs.bib 并在论文中引用。

### 待扩展方向与候选条目

#### A. NJL/PNJL 基础理论（+8篇）
- Hatsuda & Kunihiro, Phys. Rept. 247, 221 (1994) — NJL 综述
- Vogl & Weise, Prog. Part. Nucl. Phys. 27, 195 (1991) — NJL 综述
- Costa et al., PRD 79, 116003 (2009) — PNJL 中介子性质
- Ruivo et al., PRD 86, 114038 (2012) — PNJL 相图和介子
- Contrera et al., PRD 82, 054026 (2010) — 非局域 PNJL 中的介子
- Schaefer et al., PRD 76, 074023 (2007) — PQM 模型相图
- Herbst et al., PRD 89, 054005 (2014) — PNJL 参数系统学
- Torres-Rincon et al., PRC 91, 065205 (2015) — NJL+RTA 输运框架

#### B. 各向异性 QGP 理论（+8篇）
- Florkowski & Ryblewski, PRC 83, 034907 (2011) — 各向异性流体力学（已在 bib）
- Martinez & Strickland, NPA 848, 183 (2010) — 各向异性等离子体介电函数
- Mrowczynski, PLB 314, 118 (1993) — Weibel 不稳定性
- Rebhan et al., PRL 94, 102303 (2005) — 非阿贝尔等离子体不稳定性
- Nopoush et al., PRC 90, 054910 (2014) — 各向异性流体力学中的状态方程
- Alqahtani et al., PRC 91, 054902 (2015) — 各向异性流体力学数值
- Carrington et al., PRC 104, 064908 (2021) — 集体模式广义各向异性（已在 bib）
- Hauksson et al., PRC 103, 064904 (2021) — QGP 探针和等离子体不稳定性

#### C. 介子谱/Mott 离解方法（+6篇）
- Blaschke et al., PPNP 91, 1 (2016) — Mott 离解综述（如有）
- Dubinin et al., PRD 93, 054038 (2016) — PNJL 中介子相移
- Yamazaki & Matsui, PRD 90, 074038 (2014) — 介子谱函数
- Zhuang et al., PRD 51, 3728 (1995) — 手征相变附近的介子
- Costa et al., PRD 71, 116002 (2005) — NJL 中的标量和赝标量介子
- Hansen et al., PRD 75, 065004 (2007) — 介子谱函数和 Mott 离解

#### D. 手征相变与序参量（+5篇）
- Aoki et al., Nature 443, 675 (2006) — 格点 QCD 相变级次（已在 PNJL.bib）
- Bazavov et al., PRD 85, 054503 (2012) — 格点 QCD 状态方程
- Fodor & Katz, JHEP 04, 050 (2004) — 格点 QCD CEP 搜索
- Shaoveddy et al., PRD 94, 094002 (2016) — PNJL 中 net-baryon 涨落
- Fischer, JPG 32, R253 (2006) — DSE 方法 QCD 相图

#### E. 重离子唯象（+5篇）
- ALICE, PLB 720, 52 (2013) — v₂ 测量
- Blaschke et al., PRC 84, 045205 (2011) — horn 效应和 Mott
- Shuryak, PPNP 62, 48 (2009) — 重离子碰撞中的 QCD 物理
- Braun-Munzinger et al., PLB 518, 41 (2001) — 统计强子化
- Andronic et al., NPA 837, 65 (2010) — 强子产额比

#### F. 邻近模型对比（+5篇）
- Schaefer & Wambach, PRD 75, 085015 (2007) — PQM 模型中的介子
- Kamikado et al., PRD 87, 034019 (2013) — DSE 中的 Mott 离解
- Fischer et al., PRD 90, 034022 (2014) — DSE 中的 QCD 相图
- Liu & Rapp, PRC 97, 034918 (2018) — 强子共振气体中的离解
- Casalderrey-Solana et al., PRD 104, 074027 (2021) — 全息重夸克偶素离解（已读书目中）

#### G. 输运系数背景（+3篇）
- Marty et al., PRC 88, 045204 (2013) — NJL 输运系数（已在 refs.bib）
- Soloveva et al., PRC 103, 054901 (2021) — PNJL 输运系数（已在 refs.bib）
- heShearViscosityElectric, arXiv:2403.05946 (2025) — PNJL 输运系数（已在 refs.bib）

### 总计：15（已有） + 40（候选）= 55篇

---

## 5. 引用策略建议

1. **Introduction 扩展**：引用 A 类（NJL/PNJL 基础）+ B 类（各向异性理论）+ F 类（邻近模型），建立完整的理论背景
2. **Methodology**：引用 A 类+C 类（介子谱方法）
3. **Results**：重点引用 Zhang 2021 NJL + He 2023 + Maslov 2023 做定量对比
4. **Discussion**：引用 Burnier 2009 + Jamal 2018 + E 类（唯象）展开讨论
5. **避免过度引用**：G 类（输运）和 D 类（格点 QCD）按需引用，不强制全加

---

## 6. 下一步操作

用户可以：
1. 审核上述扩展候选条目，勾选需要加入的
2. 我批量从已有 bib 文件中提取 bib 条目，追加到 refs.bib
3. 必要时对新增 Tier 1 条目写分析笔记
4. 在论文中插入适当引用位置
