# 配置 TOML 权威来源审计

**日期**: 2026-03-09  
**范围**: `config/` 下全部 17 个 TOML  
**目标**: 明确哪些文件属于运行时物理/模型权威配置，哪些只是测试覆写或 CI 元配置，避免将测试参数误当成物理真值。

---

## 一、审计结论

当前仓库的 TOML 可以分为三层：

1. **运行时权威配置层**  
   仅包括 `config/physics/default.toml` 与 `config/models/*/default.toml`，以及 `config/models/pnjl/magnetic_default.toml` 这类明确承载物理模型默认参数的专用 profile。

2. **测试覆写层**  
   所有 `unittest*.toml` 与 `config/physics/unittest.toml` 都是测试/回归 profile，用来验证继承、覆写与调用路径，不是物理真值来源。

3. **CI 元配置层**  
   `config/ci/*.toml` 仅用于仓库治理与门禁，不属于模型参数来源。

因此，运行时参数权威来源应解释为：

- 共享物理常数与共享三味基线参数：`config/physics/default.toml`
- 各模型默认参数：`config/models/<model>/default.toml`
- 专项物理场景参数：如 `config/models/pnjl/magnetic_default.toml`
- 测试与 CI 文件：只对验证流程负责，不参与“物理权威来源”定义

---

## 二、权威来源判定规则

### 1. 公式文档优先用于物理模型参数

下列键应优先以 `docs/reference/formula/models/...` 为权威来源：

- NJL / PNJL / rPNJL 的 `Lambda_MeV`、`G_over_Lambda2`、`K_over_Lambda5`
- 流夸克质量 `m_ud0_MeV`、`m_s0_MeV`
- Polyakov 势参数 `T0_MeV`、`a0`、`a1`、`a2`、`b3`、`b4`
- rPNJL 扩展参数 `g1_MeV_inv8`、`g2_MeV_inv8`、`kappa`
- 磁场 PNJL 的 IMC 参数 `a`、`b`、`c`、`d`、`Lambda_QCD_MeV`

### 2. 仓库实现默认值用于非公式型运行参数

下列键不是公式文档中的核心模型参数，应视为仓库运行默认值，而不是公式真值：

- `transport_workflow.prefer_energy_aniso`
- `transport_workflow.a_builder.*`
- `label`、`version`
- `rho0_fm3`
- 测试 profile 中故意设置的 `N_color = 7`、`N_color = 9` 等覆写值

### 3. CI 配置单独治理

下列 TOML 仅由仓库规则与 CI 门禁负责：

- `config/ci/unit_skip_policy.toml`
- `config/ci/models_invokelatest_allowlist.toml`

---

## 三、运行时权威配置层

### 3.1 `config/physics/default.toml`

**角色**: 共享物理常数 + 共享三味基线参数 + transport workflow 默认值。

**参数来源判定**:

- `[physical].hbarc`, `[physical].alpha_em`  
  当前属于共享物理常数层的仓库权威值。现有 `docs/reference/formula/` 尚未单列基础常数文档，因此它们目前以该文件为运行时单一事实来源。

- `[model_shared].N_color`, `Lambda_MeV`, `G_over_Lambda2`, `K_over_Lambda5`, `m_ud0_MeV`, `m_s0_MeV`  
  其中三味 NJL / PNJL 基线参数可由以下公式文档支撑：
  - `docs/reference/formula/models/njl/NJL_core.md`
  - `docs/reference/formula/models/pnjl/Omega_各向同性.md`

- `[transport_workflow]` 与 `[transport_workflow.a_builder]`  
  属于数值工作流实现默认值，不属于公式文档中的物理参数；权威来源是工作流实现与测试，而不是模型公式文档。

### 3.2 `config/models/njl/default.toml`

**角色**: 三味 NJL 默认 profile。

**键分类**:

- `N_flavor = 3` 由三味 NJL 模型定义支撑，来源可追溯到 `docs/reference/formula/models/njl/NJL_core.md`
- `rho0_fm3 = 0.16` 是经验性的核物质参考密度，当前未在公式文档中作为核心模型参数表给出，应视为仓库约定默认值
- `label`, `version` 为元数据，不属于物理参数

### 3.3 `config/models/njl2/default.toml`

**角色**: 两味 NJL 默认 profile。

**键分类**:

- `N_color = 3`, `N_flavor = 2`, `Lambda_MeV = 631.0`, `G_over_Lambda2 = 2.14`, `m_ud0_MeV = 5.5`  
  可由 `docs/reference/formula/models/njl/NJL_core.md` 中“两味NJL参数表（Set-2f-A）”支撑
- `rho0_fm3 = 0.16` 仍属于仓库经验默认值
- `label`, `version` 为元数据

结论：`njl2/default.toml` 不是“无文档来源参数”，当前已有可接受的公式侧支撑。

### 3.4 `config/models/pnjl/default.toml`

**角色**: 各向同性 PNJL 默认 profile。

**键分类**:

- `[model].N_flavor`, `rho0_fm3` 的角色与 NJL 默认 profile 相同
- `[polyakov].T0_MeV`, `a0`, `a1`, `a2`, `b3`, `b4` 的权威来源为：
  - `docs/reference/formula/models/pnjl/Omega_各向同性.md`
  - `docs/reference/formula/models/pnjl/Polyakov有效势_多参数化方案.md`

### 3.5 `config/models/pnjl/magnetic_default.toml`

**角色**: 磁场 PNJL 专用 profile。

**键分类**:

- `[magnetic].eB_fm2`, `n_max`, `p_num`, `pz_max`, `cutoff_N`  
  属于磁场计算场景的运行参数；其中 `eB_fm2` 对应 `docs/reference/formula/models/pnjl_magnetic/PNJL_magnetic_core.md` 的磁场变量定义，积分截断与数值离散参数属于实现级默认值

- `[magnetic.imc].a`, `b`, `c`, `d`, `Lambda_QCD_MeV`  
  权威来源为 `docs/reference/formula/models/pnjl_magnetic/PNJL_magnetic_core.md` 中 IMC 经验公式

### 3.6 `config/models/rpnjl/default.toml`

**角色**: rPNJL 默认 profile。

**键分类**:

- `[rpnjl].Lambda_MeV`, `G_over_Lambda2`, `K_over_Lambda5`, `m_ud0_MeV`, `m_s0_MeV`, `g1_MeV_inv8`, `g2_MeV_inv8`, `kappa`, `T0_MeV`, `a0`, `a1`, `a2`, `b3`, `b4`  
  权威来源为 `docs/reference/formula/models/rpnjl/rPNJL_core.md`

结论：`rpnjl/default.toml` 当前已经与公式文档注释保持一致，适合作为运行时唯一默认 profile。

---

## 四、测试覆写层

以下文件都不是物理真值来源，只用于测试继承、覆写与路径选择：

- `config/physics/unittest.toml`
- `config/models/njl/unittest.toml`
- `config/models/njl/unittest_oldonly.toml`
- `config/models/njl2/unittest.toml`
- `config/models/pnjl/unittest.toml`
- `config/models/pnjl/unittest_gk_polyakov.toml`
- `config/models/pnjl/unittest_lambda.toml`
- `config/models/pnjl/unittest_oldonly.toml`
- `config/models/rpnjl/unittest.toml`

说明：

- 这些文件中的值只需要满足“可验证覆写逻辑”和“可稳定触发测试路径”。
- 它们可以故意偏离公式文档参数，例如更改 `N_color`、`Lambda_MeV` 或 Polyakov 参数。
- `config/models/rpnjl/unittest.toml` 虽然为空覆写，但语义仍属于测试 profile，而不是默认权威配置。

---

## 五、CI 元配置层

以下文件不应纳入物理配置审计争议：

- `config/ci/unit_skip_policy.toml`
- `config/ci/models_invokelatest_allowlist.toml`

它们的权威来源分别是：

- 测试治理策略
- 迁移期代码治理策略

因此，T5-04 的“权威来源”不应错误地要求它们服从公式文档。

---

## 六、全量文件分类清单

| 文件 | 分类 | 权威来源 |
|------|------|----------|
| `config/physics/default.toml` | 运行时权威配置 | 共享物理常数文件本身 + `docs/reference/formula/models/njl/NJL_core.md` + `docs/reference/formula/models/pnjl/Omega_各向同性.md` |
| `config/physics/unittest.toml` | 测试覆写 | 测试需求 |
| `config/models/njl/default.toml` | 运行时权威配置 | `docs/reference/formula/models/njl/NJL_core.md` + 仓库元数据约定 |
| `config/models/njl/unittest.toml` | 测试覆写 | 测试需求 |
| `config/models/njl/unittest_oldonly.toml` | 测试覆写 | 测试需求 |
| `config/models/njl2/default.toml` | 运行时权威配置 | `docs/reference/formula/models/njl/NJL_core.md` 两味补充参数表 |
| `config/models/njl2/unittest.toml` | 测试覆写 | 测试需求 |
| `config/models/pnjl/default.toml` | 运行时权威配置 | `docs/reference/formula/models/pnjl/Omega_各向同性.md` + `docs/reference/formula/models/pnjl/Polyakov有效势_多参数化方案.md` |
| `config/models/pnjl/magnetic_default.toml` | 运行时权威配置 | `docs/reference/formula/models/pnjl_magnetic/PNJL_magnetic_core.md` + 数值实现默认值 |
| `config/models/pnjl/unittest.toml` | 测试覆写 | 测试需求 |
| `config/models/pnjl/unittest_gk_polyakov.toml` | 测试覆写 | 测试需求 |
| `config/models/pnjl/unittest_lambda.toml` | 测试覆写 | 测试需求 |
| `config/models/pnjl/unittest_oldonly.toml` | 测试覆写 | 测试需求 |
| `config/models/rpnjl/default.toml` | 运行时权威配置 | `docs/reference/formula/models/rpnjl/rPNJL_core.md` |
| `config/models/rpnjl/unittest.toml` | 测试覆写 | 测试需求 |
| `config/ci/unit_skip_policy.toml` | CI 元配置 | 测试治理策略 |
| `config/ci/models_invokelatest_allowlist.toml` | CI 元配置 | 迁移治理策略 |

---

## 七、后续建议

1. 若希望 `config/physics/default.toml` 的 `[physical]` 常数也具备文档侧权威来源，建议后续补一份 `docs/reference/formula/constants/` 或 `docs/reference/domain-knowledge/constants.md`。
2. `rho0_fm3` 当前属于仓库约定默认值；若后续它被用于更严格的物理解释，建议补一份参考说明文档。
3. 审计后不建议把测试 profile 的数值“校正回公式值”，否则会削弱继承/覆写测试的有效性。