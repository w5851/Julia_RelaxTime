# Mott 验证与 xi 扫描 Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在不重写核心求解器的前提下，完成 isotropic Mott/介子质量验证测试与 `mu=0` 的 anisotropic `xi` 扫描脚本，形成可复现输出。

**Architecture:** 复用 `src/models/workflows/MesonMassWorkflow.jl` 作为单一计算入口；新增参考数据映射层负责数据清洗与单位统一；测试负责口径锁定；分析脚本负责扫描与结果落盘。通过状态标记和 metadata 提升结果可追溯性。

**Tech Stack:** Julia 1.10+, Test stdlib, 现有 Models/RelaxTime 模块, CSV 输出（按项目现有依赖与风格）

---

## Chunk 1: 参考数据映射与口径函数

### Task 1: 新增参考数据映射工具

**Files:**
- Create: `scripts/analysis/mott_reference_mapping.jl`
- Modify: `docs/dev/active/2026-03-20_Mott温度与介子质量验证及各向异性扫描任务单.md`
- Spec: `docs/superpowers/specs/2026-03-20-mott-validation-xi-scan-design.md`

- [ ] **Step 1: 写最小失败测试（映射层）**

在 `tests/validation/` 新建对应测试文件骨架，先断言字段缺失会抛错。

- [ ] **Step 2: 运行单测确认失败**

Run: `julia --project=. -e 'include("tests/validation/<new_mapping_test>.jl")'`
Expected: FAIL（缺少实现或字段验证）。

- [ ] **Step 3: 实现最小映射函数**

实现函数建议：
- `load_reference_table(path)::Vector{NamedTuple}`
- `normalize_reference_units(rows; mass_unit=:MeV, temp_unit=:MeV)`
- `validate_reference_schema(rows, required_cols)`

并补充：
- 验证通道按“老程序理论可计算能力”锁定，不依赖当前导出文件是否已存在。
- 首批目标通道：`pi`, `K`, `eta`, `eta_prime`, `sigma_pi`, `sigma_K`, `sigma`, `sigma_prime`。

### Task 1b: 补齐老程序介子量导出规范

**Files:**
- Create: `docs/dev/active/2026-03-20_legacy_meson_export_spec.md`（或并入现有任务单）
- External: `D:/Desktop/Fortran/NJL-ReTime/...`（记录路径，不在本仓库直接改）

- [ ] **Step 1: 定义导出 schema（字段、单位、分隔符、文件命名）**
- [ ] **Step 2: 定义最小验证样例点（T, muB, xi）**
- [ ] **Step 3: 明确 Fortran 导出入口函数位置**
- [ ] **Step 4: 生成导出规范文档并回链到主任务单**

当前状态：

- [x] 导出规范文档已落地：`docs/dev/active/2026-03-20_legacy_meson_export_spec.md`
- [ ] 待你在 Fortran 侧按该规范补导出后，再回填到 Julia 仓库 `tests/validation/data/` 对应目录。

- [ ] **Step 4: 运行测试确认通过**

Run: `julia --project=. -e 'include("tests/validation/<new_mapping_test>.jl")'`
Expected: PASS。

- [ ] **Step 5: 提交本任务变更**

```bash
git add scripts/analysis/mott_reference_mapping.jl tests/validation/<new_mapping_test>.jl docs/dev/active/2026-03-20_Mott温度与介子质量验证及各向异性扫描任务单.md
git commit -m "feat: add reference mapping helpers for mott validation"
```

### Task 2: 固化 Mott 温度判据函数

**Files:**
- Create: `scripts/analysis/mott_criterion_utils.jl`
- Test: `tests/validation/<new_mott_criterion_test>.jl`

- [ ] **Step 1: 写失败测试（过零/最近点逻辑）**
- [ ] **Step 2: 运行确认失败**
- [ ] **Step 3: 实现最小判据函数**

建议函数：
- `estimate_mott_temperature(T_grid, gap_grid; method=:zero_crossing_or_minabs)`

建议扩展：
- `refine_mott_temperature_bisection(gap_func, T_lo, T_hi; tol=1e-4, max_iter=50)`（仅在括号区间可用时调用）
- 统一返回 `method_used`（`linear`, `bisection`, `minabs_approx`）以便结果审计。

- [ ] **Step 4: 运行测试确认通过**
- [ ] **Step 5: 提交本任务变更**

```bash
git add scripts/analysis/mott_criterion_utils.jl tests/validation/<new_mott_criterion_test>.jl
git commit -m "feat: add deterministic mott temperature criterion utilities"
```

---

## Chunk 2: 各向同性验证测试落地

### Task 3: 新增 isotropic 对照验证测试

**Files:**
- Create: `tests/validation/relaxtime/test_mott_meson_validation_isotropic.jl`
- Modify: `tests/validation/runtests.jl`
- Read/Reuse: `src/models/workflows/MesonMassWorkflow.jl`

- [ ] **Step 1: 写失败测试（读取参考 + 计算 + 误差断言）**
- [ ] **Step 2: 运行单文件测试确认失败**

Run: `julia --project=. -e 'include("tests/validation/relaxtime/test_mott_meson_validation_isotropic.jl")'`
Expected: FAIL（接口未接好或阈值未配置）。

- [ ] **Step 3: 接入 workflow 与映射层完成最小实现**

测试至少包含：
- `xi=0, mu=0`
- 通道：`pi`, `K`（可扩展 `eta`, `eta_prime`）
- `atol/rtol` 双阈值断言
- 失败时输出：温度点、通道、参考值、计算值、误差

- [ ] **Step 4: 运行测试确认通过**

Run: `julia --project=. tests/validation/runtests.jl`
Expected: 新增测试通过，不破坏现有用例。

- [ ] **Step 5: 提交本任务变更**

```bash
git add tests/validation/relaxtime/test_mott_meson_validation_isotropic.jl tests/validation/runtests.jl
git commit -m "test: validate isotropic mott and meson masses against references"
```

---

## Chunk 3: 各向异性 xi 扫描脚本

### Task 4: 新增 xi 扫描脚本（mu=0）

**Files:**
- Create: `scripts/analysis/scan_mott_meson_vs_xi_mu0.jl`
- Optionally Create: `config/analysis/mott_xi_scan_mu0.toml`（若仓库已有同类配置习惯）

- [ ] **Step 1: 写最小失败测试或 smoke 驱动检查**

若不便写完整单测，至少写一个可执行 smoke 检查：
- 小网格 `xi=[-0.1,0.0,0.1]`
- 断言输出 CSV 存在且字段齐全

- [ ] **Step 2: 运行确认失败**

Run: `julia --project=. scripts/analysis/scan_mott_meson_vs_xi_mu0.jl --smoke`
Expected: FAIL（脚本未实现或字段缺失）。

- [ ] **Step 3: 实现脚本最小可用版本**

功能要求：
- 调用 `solve_gap_and_meson_point`
- 固定 `mu=0`
- 可配置 `xi` 网格与温度区间（默认 `xi_min=-0.4`, `xi_max=0.4`, `xi_step=0.05`）
- 结果写出 CSV
- 输出状态字段：`success/fallback/fail`

- [ ] **Step 4: 运行 smoke 与常规扫描**

Run:
- `julia --project=. scripts/analysis/scan_mott_meson_vs_xi_mu0.jl --smoke`
- `julia --project=. scripts/analysis/scan_mott_meson_vs_xi_mu0.jl`

Expected: 两者均完成，smoke 快速通过，常规扫描可产出结果。

- [ ] **Step 5: 提交本任务变更**

```bash
git add scripts/analysis/scan_mott_meson_vs_xi_mu0.jl config/analysis/mott_xi_scan_mu0.toml
git commit -m "feat: add mu0 anisotropy xi scan for mott and meson masses"
```

---

## Chunk 4: 文档、治理与回归闭环

### Task 5: 更新开发文档与执行记录

**Files:**
- Modify: `docs/dev/active/2026-03-20_Mott温度与介子质量验证及各向异性扫描任务单.md`
- Modify: `docs/superpowers/specs/2026-03-20-mott-validation-xi-scan-design.md`（若实现后有偏差）

- [ ] **Step 1: 回填已完成项与命令记录**
- [ ] **Step 2: 记录风险与已知限制（失败点、fallback 比例）**
- [ ] **Step 3: 运行文档治理检查**

Run: `julia --project=. scripts/dev/check_active_docs_governance.jl`
Expected: OK。

- [ ] **Step 4: 提交文档变更**

```bash
git add docs/dev/active/2026-03-20_Mott温度与介子质量验证及各向异性扫描任务单.md docs/superpowers/specs/2026-03-20-mott-validation-xi-scan-design.md
git commit -m "docs: align mott validation task doc with implementation status"
```

### Task 6: 最终验证与汇总

**Files:**
- No new files required; verify all modified paths

- [ ] **Step 1: 运行最小验证矩阵**

Run:
- `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
- `julia --project=. tests/validation/runtests.jl`

- [ ] **Step 2: 核验脚本产物路径与字段**
- [ ] **Step 3: 输出最终变更摘要与后续建议**

- [ ] **Step 4: 最终提交（如存在未提交变更）**

```bash
git add <remaining-files>
git commit -m "refactor: finalize mott validation and xi scan workflow integration"
```

---

## 执行注意事项

- 仅在用户明确要求时执行 `git push`。
- 若参考数据包含缺失/异常值，测试应 fail-fast 并给出具体字段。
- 若扫描出现不收敛点，不应中断全局任务，必须落盘 `status` 与错误摘要。
- 所有新增脚本放在 `scripts/analysis/`，避免放到 `tests/` 目录。
