---
title: Dependency Analysis & Recommendations
archived: true
original: docs/dev/任务2.md
archived_date: 2026-01-19
---

概要：

- 本文档总结了为本项目引入模块依赖可视化与审计的价值、实践建议与分层方案（L1/L2/L3），并给出基于仓库现状的落地建议与最小可行方案（MVP）。
- 已在仓库中实现：`scripts/dev/gen_deps.jl`（生成 Mermaid）、`scripts/dev/analyze_deps.jl`（审计/违规检测）、`docs/architecture/dependencies.manual.md`（L1/L3 手工内容）、CI 工作流 `.github/workflows/dependency-audit.yml`（严格审计）。

建议后续行动（优先级）：

- 确认并固定 `docs/architecture/dependency_rules.md` 中的允许依赖矩阵（若需调整，更新规则文档）。
- 在 CI 中根据团队共识决定是否开启严格模式（`DEPS_STRICT=1`）。
- 对审计报告中的违规点优先做小范围重构或将共享逻辑上移到 `src/utils/` 或 `src/core/`。

以下为原始内容（保留，以便审阅与历史参考）：

***

根据以下分析内容，结合项目本身分析适合本项目的依赖分析工具：
**是的，这是一个非常好的实践，但具体实现方式需要根据项目情况灵活选择。** 记录模块依赖关系对于项目的可维护性、新人上手和架构理解至关重要，关键在于如何高效、可持续地维护这份文档。

### 为什么值得做？（价值）

1.  **架构清晰化**：可视化依赖关系有助于理解系统分层、职责划分和核心数据流。
2.  **降低认知成本**：新成员能快速了解代码结构，避免在迷宫中摸索。
3.  **影响分析**：修改某个模块时，能快速评估影响范围，避免 unintended consequences。
4.  **发现坏味道**：能识别出循环依赖、过度耦合、上帝模块等问题，驱动重构。
5.  **文档化设计决策**：依赖图本身就是架构设计的一种体现。

### 如何实现？（实践建议）

**优先考虑自动化生成，辅以必要的手动文档。**

#### 1. **自动化工具（推荐首选）**
   - **代码分析工具**：使用工具直接从源码生成依赖图，保证准确性且可持续更新。
     - **JavaScript/TypeScript**：`madge`、`dependency-cruiser`、`codeviz`
     - **Java**：`Structure101`、`JDepend`、`ArchUnit`（用于校验）
     - **Python**：`pydeps`、`snakefood`
     - **通用/多语言**：`CodeSee`、`Sourcegraph`、`Lattix`
   - **构建工具集成**：一些构建系统（如 Gradle/Maven）可生成模块依赖报告。
   - **IDE 插件**：许多 IDE 支持生成局部的依赖图。

#### 2. **半自动文档**
   - **使用文本/图表代码**：用代码定义依赖图，可版本化管理。
     - **Mermaid**（推荐）：在 Markdown 中绘制，易维护。
       ```mermaid
       graph TD
           A[前端UI模块] --> B[API服务层]
           B --> C[数据访问层]
           C --> D[(数据库)]
       ```
     - **Graphviz/DOT**：更强大但学习曲线略高。
   - **结合架构即代码工具**：如 `Structurizr`，可将架构定义为代码并生成多视图文档。

#### 3. **手动文档（谨慎使用）**
   - 适用于小型项目或高层抽象模块图。
   - **必须与代码同步更新**：否则很快过时，比没有文档更有害。
   - 建议作为自动生成图的**补充说明**，解释设计意图、关键接口等。

### 最佳实践组合方案

1.  **分层记录**：
    - **高层架构图**（L1）：手绘/工具绘制的子系统/微服务级依赖，说明业务边界。
    - **模块依赖图**（L2）：**自动化生成**的 `src/` 内模块/包级关系。
    - **关键类/接口依赖**（L3）：对核心复杂部分，可生成类图或补充说明。

2.  **集成到开发流程**：
    - 将自动生成的依赖图嵌入到项目 `README` 或文档站点。
    - 在 CI 中运行依赖规则检查（如禁止循环依赖）。
    - **变更时更新**：将依赖图更新作为架构重大改动的一部分。

3.  **内容建议**：
    - 标注**依赖方向**（谁依赖谁）。
    - 区分**编译时依赖**与**运行时依赖**。
    - 注明**关键接口/协议**。
    - 对于循环依赖、技术债处添加注释。

### 何时可以不做？

- **极小型项目**（几个文件）依赖关系一目了然。
- **原型阶段**：结构变化过于频繁，维护成本过高。
- **已有更优替代**：如果项目有完善的接口文档、服务契约（如 OpenAPI）且团队已习惯使用。

### 推荐方案

**对于绝大多数生产级项目：**
1.  设置一个自动化脚本/任务，一键生成最新的模块依赖图（如使用 `madge` 或 `pydeps`）。
2.  将生成的 SVG/PNG 或 Mermaid 代码放入 `/docs/architecture/dependencies.md`。
3.  在 PR 模板中提醒：若模块关系发生变化，请重新生成依赖图。
4.  在 README 中提供指向该文档的链接。

这样你既获得了依赖可视化的好处，又避免了手动维护的负担，使得文档能够**随代码共同演进**。

**总结：不是“应不应该”，而是“如何有效地做”**。自动化生成的可视化依赖图是现代软件工程中一项高性价比的实践，值得投入。

---

## 结合本项目的实践落地（建议执行顺序）

> 目标：在不引入重型依赖的前提下，先覆盖 `src/` 的 Julia 模块依赖，再逐步补齐 `web/` 与脚本层，最终形成“可自动更新 + 可读”的依赖文档。

### 1) 依赖图分层与覆盖范围（对应 L1/L2/L3）

- **L1 高层架构图（手动 Mermaid）**
  - 覆盖：`src/`、`scripts/server/`、`web/`、`tests/`、`docs/`、`data/`、`config/`
  - 建议放置：`docs/architecture/dependencies.md`
  - 目的：新人快速理解“主链路 + 辅助资产”的关系。

- **L2 模块依赖圖（自动生成，Julia 为主）**
  - 重点：`src/` 内模块（`pnjl/`、`relaxtime/`、`integration/`、`simulation/`、`utils/`）
  - 依赖来源：
    - `include("...")` 链路（文件级）
    - `using .Module` / `import .Module`（模块级）
  - 输出形式：Mermaid `graph TD`，按目录分组（subgraph）

- **L3 关键链路补充（半自动/手动）**
  - 重点：
    - 计算链路：`ScatteringAmplitude → DifferentialCrossSection → TotalCrossSection → AverageScatteringRate → RelaxationTime`
    - PNJL：`SeedStrategies → Solver → Scan`
  - 可直接复用 README 中“计算链路概览”的描述作为注释。

### 2) 自动化生成策略（不引入新依赖）

**推荐：新增一个轻量 Julia 脚本用于解析依赖并输出 Mermaid。**

- 位置建议：`scripts/dev/gen_deps.jl`
- 解析规则（足够覆盖核心场景）：
  - `include("...")` 视为文件依赖
  - `using .X` / `import .X` 视为模块依赖
  - 忽略第三方包（仅关注仓库内部）
- 输出：`docs/architecture/dependencies.md`（包含 L2 Mermaid 图）

> 优点：零额外包依赖、可在 Windows 直接运行、适配当前 Julia 项目形态。

### 3) 与仓库规范对齐

- **文档放置**：符合 [docs/dev/项目结构约定.md](docs/dev/%E9%A1%B9%E7%9B%AE%E7%BB%93%E6%9E%84%E7%BA%A6%E5%AE%9A.md) 的 `docs/` 规则。
- **更新时机**：
  - 结构性重构（如 `src/` 模块拆分）后生成
  - PR 模板加入提示：若模块关系变更，需更新依赖图

### 4) 可选增强（后续再做）

- **Web 层依赖**：后续可补充 `web/js/` 的 ES Module `import` 解析，生成单独图。
- **外部依赖清单**：用 Julia 内置 `Pkg.dependencies()` 输出第三方包清单，放到 `docs/architecture/external_deps.md`。
- **循环依赖检查**：在自动生成脚本中检测强连通分量，并在图上标注。

---

## 对本项目的最小可行版本（MVP）

1. 新增 `docs/architecture/dependencies.md`
2. 手写 L1 高层图（1 页）
3. 使用 `scripts/dev/gen_deps.jl` 自动生成 L2 Mermaid 图
4. 在 README 增加链接指向依赖文档

这样即可满足“可持续维护 + 可读可用”的实践要求。