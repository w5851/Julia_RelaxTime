# ADR-0007: 公式路线闭合门禁与 charged-RPA/BU 试点

## 状态

已接受

## 背景

项目在 PNJL、介子传播子、BU 相移和输运等领域同时使用项目内公式与外部文献
约定。不同来源可能对凝聚符号、KMT 有效耦合、极化归一化、传播子分母、相位
分支和单位采用不同记号。若只把某一段公式或代码接入，而没有从微观模型追到
最终输出的闭合链，数值上“能运行”并不能说明路线自洽。

本轮 charged-RPA/BU 审查具体暴露出：`K_{45}`/`K_{67}` 的味道映射、旧
`H_f=-phi_f` helper、scalar `1-4KPi` 与 neutral matrix `1-2KPi` 的差异、
`num_s_quark=1` 对三动量截断的兼容处方，以及 `x_min_cut` 的诊断边界，需要
放在同一个可逐项质询的规范中。

## 决策

1. 在 `docs/guides/sop/formula_route_closure.md` 建立跨领域的公式路线闭合 SOP。
2. 在 `config/governance/formula_route_closure.toml` 登记每条路线的文档、状态、
   外部来源、代码/测试映射和未决项，并用
   `scripts/dev/check_formula_route_closure.jl` 做结构化检查。
3. 新领域或新生产候选默认遵循
   `draft -> candidate -> production_authorized`；检查器只阻止缺失/矛盾，不自动
   进行人工升格。
4. 以 `ChargedRPA_BU_ProductionRoute.md` 作为首个试点：固定 BQS quark-only
   背景，明确 `pi^\pm -> K12`、`K^\pm -> K45`，并用显式 charged ladder trace
   与 Goldstone 条件闭合 scalar/matrix 分母归一化；严格 ordered retarded bubble、
   相位/Levinson、单电荷 BU 测度实现、Bose 支撑和数值收敛仍保留为 candidate
   未决项，second-sheet 极点只作为 strict-BW 路线的 oracle。
5. 本决策只增加治理和公式审查层，不修改 `PNJLCore`、稳定 `MesonDensity` 默认
   语义、transport 路线或正式数值 baseline。

## 理由

- 将公式、来源、单位、代码和测试放在同一可审阅包中，可以在实现前发现符号与
  归一化断裂。
- 机读登记让未来新领域可重复执行同一最低门槛，同时保持不同专题 SOP 的数值
  收敛与产物治理独立。
- charged-RPA 试点足以验证流程，但不把尚未完成的严格数值层包装成生产结果。

## 备选方案

### 方案A：现在创建独立 Codex skill

**不选择原因**：目前只有一个跨领域试点；先用 SOP、registry 和 checker 验证流程，
待至少三个新领域重复使用后再评估是否自动化为 skill。

### 方案B：把全部公式写入一个总 SOP

**不选择原因**：总 SOP 会混合模型、输运、介子和绘图的职责；公式规范应留在
`docs/reference/formula/`，专题执行仍由各自 SOP 负责。

### 方案C：只依赖代码 review 和现有文献链接

**不选择原因**：无法保证每个外部公式的方程号、单位转换、替代约定和未决项都
被记录，也无法由 CI/本地检查快速发现文档与登记脱节。

## 后果

**正面影响**：

- 新路线必须明确模型起点、近似层次、来源、代码、测试和 production 边界；
- 可在不运行昂贵求解器的情况下先发现路径/公式登记错误；
- 便于独立审阅者按固定快照复核，保留失败和 diagnostic 证据。

**负面影响**：

- 新领域初期需要编写额外的公式和 provenance 文档；
- checker 只能验证结构和文字标记，不能替代物理推导或数值验证。

**风险与缓解**：

- 标记完整但公式仍有误：要求独立人工审阅、维度/极限测试和外部固定点；
- 路线长期停在 candidate：保留未决项清单，禁止更新 production baseline；
- 来源版本变化：在快照中记录 DOI/arXiv、日期和基线 SHA，必要时创建新 route id。

## 相关决策

- [ADR-0005](0005-scientific-sop-and-document-authority.md) - 科学计算 SOP 与文档权威治理
- [公式路线闭合 SOP](../guides/sop/formula_route_closure.md)
- [charged-RPA/BU 公式路线](../reference/formula/relaxtime/ChargedRPA_BU_ProductionRoute.md)

## 参考资料

- Rehberg, Klevansky, Hüfner, *Nucl. Phys. A* 608 (1996), DOI [10.1016/0375-9474(96)00247-3](https://doi.org/10.1016/0375-9474(96)00247-3)
- Rehberg & Klevansky, *Ann. Phys.* 268 (1998), DOI [10.1006/aphy.1996.0140](https://doi.org/10.1006/aphy.1996.0140)
- Tian et al., *Phys. Rev. D* 114 (2026), DOI [10.1103/d7nm-y2vp](https://doi.org/10.1103/d7nm-y2vp)
- Blaschke et al., *Phys. Rev. D* 96 (2017), DOI [10.1103/PhysRevD.96.094008](https://doi.org/10.1103/PhysRevD.96.094008)
- Blaschke et al., *Particles* 3 (2020), DOI [10.3390/particles3010014](https://doi.org/10.3390/particles3010014)

---

**日期**：2026-08-30
**作者**：项目维护者
**审查者**：待独立审阅（PR290 后续对话）
