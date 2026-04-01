# Models 固定维度耦合点基线清单（Gate A）

## 1. 目的与范围

本清单作为 Gate A 开工前基线附件，用于冻结当前固定维度耦合点，避免迁移过程中范围漂移。

覆盖范围：

- `src/models/state.jl`
- `src/models/solver/Conditions.jl`
- `src/models/solver/ImplicitSolver.jl`
- `src/models/constraint_solver.jl`
- `src/models/scans/TmuScan.jl`
- `src/models/scans/TrhoScan.jl`

---

## 2. 基线快照（v0）

### 2.1 状态与化学势硬编码契约

- `src/models/state.jl:40`：`x_state length must be 3 or >= 5`（状态维度硬约束）。
- `src/models/state.jl:65`：`state_vector(st::MeanFieldState) -> SVector{5}`（固定 5 维输出）。
- `src/models/state.jl:91`：`mu_vec must have length 3`（固定 3 味约束）。

### 2.2 Solver 条件组装硬切片

- `src/models/solver/Conditions.jl:134`：`gap_conditions(x_state::SVector{5,...})`。
- `src/models/solver/Conditions.jl:191`：`x_state = SVector{5}(x[1],...,x[5])`。
- `src/models/solver/Conditions.jl:192`：`mu_vec = SVector{3}(x[6],x[7],x[8])`。
- `src/models/solver/Conditions.jl:332`：`F[1:5] = gap_conditions(...)`（残差前 5 维固定语义）。

### 2.3 ImplicitSolver 结果与中间量固定维度

- `src/models/solver/ImplicitSolver.jl:223`：`x_state::SVector{5, Float64}`。
- `src/models/solver/ImplicitSolver.jl:224`：`mu_vec::SVector{3, Float64}`。
- `src/models/solver/ImplicitSolver.jl:230`：`masses::SVector{3, Float64}`。
- `src/models/solver/ImplicitSolver.jl:419`：`x_state = SVector{5}(Tuple(x_sol[1:5]))`。
- `src/models/solver/ImplicitSolver.jl:420`：`mu_vec = SVector{3}(x_sol[6],x_sol[7],x_sol[8])`。

### 2.4 Constraint solver 固定切片与固定向量构造

- `src/models/constraint_solver.jl:16`：`seed[6:8]`。
- `src/models/constraint_solver.jl:147`：`x_state = SVector{5}(Tuple(state_vector(st)))`。
- `src/models/constraint_solver.jl:158`：`calculate_mass_vec(..., SVector{3}(x_state[1],x_state[2],x_state[3]))`。
- `src/models/constraint_solver.jl:319`：`x_state=SVector{5}(fill(NaN, 5))`。
- `src/models/constraint_solver.jl:320`：`mu_vec=SVector{3}(fill(NaN, 3))`。

### 2.5 Scan 结果回写固定维度强制转换

- `src/models/scans/TmuScan.jl:494`：`SVector{5}(result.x_state)`。
- `src/models/scans/TmuScan.jl:495`：`SVector{3}(result.mu_vec)`。
- `src/models/scans/TrhoScan.jl:659`：`SVector{5}(result.x_state)`。
- `src/models/scans/TrhoScan.jl:660`：`SVector{3}(result.mu_vec)`。

---

## 3. 与 7 模型覆盖关系

当前固定维度耦合主要沿 PNJL 语义展开（5+3）；
对 `NJL2/Rotation/GasLiquid` 等非同构状态语义模型，属于高风险迁移点。

Gate A 要求：

- 所有上述耦合点需迁移到 schema-driven 契约；
- 迁移后不再新增同类硬编码点；
- 7 模型至少各有一条主链路通过验证。

---

## 4. 迁移后判定规则（用于对照）

- 不再出现新增 `x[1:5]` / `x[6:8]` / `SVector{5}` / `SVector{3}` 作为主流程强约束。
- 状态维度由 `ModelStateSchema` 派生，而非写死常量。
- `SolverResult`/scan 回写不再强制固定维度转换。

---

## 5. 基线生成命令（复核用）

```powershell
julia --project=. -e 'include("src/constants/Constants_PNJL.jl"); include("src/models/Models.jl"); println(Models.registered_model_kinds())'
```

```powershell
rg -n "SVector\{5|SVector\{3|\[1:5\]|\[6:8\]|fill\(NaN, 5\)|fill\(NaN, 3\)|x_state length must be 3 or >= 5|mu_vec must have length 3" src/models
```

---

## 6. 执行日志

### 2026-04-01

- [x] 建立 Gate A 固定维度耦合点基线清单（v0）。
