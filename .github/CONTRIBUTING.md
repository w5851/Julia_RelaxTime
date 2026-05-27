# 贡献指南（Contributing）

欢迎贡献。这个仓库包含科研/数值计算相关代码（Julia 为主），我们非常重视：可复现、单位一致性、以及对现有 API 的最小破坏。

## 开始之前

- **单位约定**：项目内部统一使用自然单位制；面向外部的 MeV 字段和内部 `fm^-1` 字段需按 `AGENTS.md` 与 `docs/dev/README.md` 的命名约定显式标注。
- **范围意识**：`src/` 里的模块用于可复用逻辑；探索性/一次性脚本请放到 `scripts/analysis/`，性能探针放到 `scripts/perf/` 或 `benchmark/`，非测试脚本不要放入 `tests/`。
- **入口优先级**：新增稳定调用方优先使用 `Models` 与 `src/models/entrypoints.jl`；历史兼容路径只作为迁移说明或对照验证使用。

## 开发环境

Julia 项目使用 `Project.toml` / `Manifest.toml` 管理依赖。

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

## 如何运行测试

- 本地改动优先运行对应分层 smoke 或单文件选择器，例如：

```powershell
julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'
julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'
```

- 运行单个测试脚本（示例）：

```powershell
julia --project=. -e 'include("tests/unit/pnjl/test_solver_seed_strategies.jl")'
```

- 文档或入口治理相关改动建议同时运行：

```powershell
julia --project=. scripts/dev/check_docs_consistency.jl
julia --project=. scripts/dev/check_script_entrypoints.jl
julia --project=. scripts/dev/check_models_entry_contract.jl
```

- 全量 package test 仍可作为更宽的补充验证：

```powershell
julia --project=. -e 'using Pkg; Pkg.test()'
```

请确保你的改动至少覆盖相关单元测试，或在 PR 中说明为什么无法测试。

## 提交与分支

- 建议使用清晰的 topic 分支；自动化协作分支可使用 `codex/<topic>`，人工 feature/fix 分支可使用 `feature/<topic>` 或 `fix/<topic>`。
- 提交信息建议包含清晰动词与影响范围，例如：
  - `pnjl: fix seed_state propagation in workflow`
  - `docs: update PNJL api docs`

## Pull Request（PR）期望

PR 请尽量包含：
- 问题背景/动机（为什么要改）
- 改动摘要（改了什么）
- 验证方式（跑了哪些测试/脚本、输出在哪里）
- 若涉及物理/数值假设：给出推导或参考文档链接（放到 `docs/`）

### 变更边界（重要）

- **尽量不要**在同一个 PR 中混合“功能改动 + 大范围重排/格式化”。
- 稳定入口、脚本契约、`Models` 公共 API 或前后端数据契约变更需要：
  - 更新对应 `docs/api/` 文档
  - 必要时同步更新 `docs/guides/`、`docs/dev/`、`docs/notes/` 或 `CHANGELOG.md`

## 文档贡献

- API 文档位于 `docs/api/`。
- 我们不强制要求 `docs/api/` 与 `src/` 逐文件 1:1 镜像；更推荐按“稳定模块/公共入口”组织。
- 若你新增模块或重构目录结构，请同步更新：
  - `docs/dev/README.md`
  - 受影响的 `docs/api/*` 页面

## 数据与结果

- 原始数据与运行输出一般不建议直接提交到 Git（尤其是大文件）。
- 若必须提交：请说明来源、生成方式与尺寸，并在 PR 描述中标注。

## 行为准则
参与贡献即表示你同意遵守本仓库的行为准则：`.github/CODE_OF_CONDUCT.md`。
