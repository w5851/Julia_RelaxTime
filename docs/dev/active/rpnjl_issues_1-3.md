## rPNJL Issue 清单（1-3）

简要：以下为可直接分配的三个 issue，覆盖参数注册、八夸克质量项与 Polyakov 势（含 Vandermonde）。每项包含修改要点、测试要点与估时。

---

### Issue 1 — 添加 rPNJL 参数常量
- 摘要：在配置与常量文件中注册 rPNJL 所需参数，并实现 MeV->fm 的加载转换。  
- 修改文件：
  - `config/pnjl/default.toml`（新增 `g1_MeV`, `g2_MeV`, `kappa`, `T0_MeV`, `a0,a1,a2,b3,b4`）
  - `src/Constants_PNJL.jl`（新增常量、加载与转换函数 `load_rpnjl_params!(cfg)`）
- 主要改动说明：
  - 在配置解析处把以 `MeV` 标记的参数转换为内部 `fm^-1`（或与项目内部单位一致）。  
  - 增加注释与单元，明确默认值（参见文档表3.1/3.2）。
- 测试：加载不同 `default.toml` 配置，断言常量已注册且数值单位转换正确（用已知换算因子检查）。
- 估时：0.5 天
- 验收条件：`Constants_PNJL` 能以配置覆盖默认值并在核心计算前提供 `g1,g2,kappa,...`（单位为内部单位）。

---

### Issue 2 — 实现八夸克质量项（式3.31）
- 摘要：在质量计算中加入 `g1,g2` 的贡献，使自洽质量方程匹配文档。  
- 修改文件（建议位置）：
  - `src/pnjl/core/Core.jl` 或 `src/pnjl/core/MassEq.jl`（若存在）
  - 相关测试目录 `tests/unit/` 增加 `mass_eq_rpnjl_test.jl`
- 主要改动说明：
  - 在 `calculate_mass_vec(σ::SVector{3})` 内实现
    ```text
    M_f = m_f - 2*gS*σ_f + (gD/4)*σ_{f+1}σ_{f+2} - 2*g1*σ_f*(Σσ_{f'}^2) - 4*g2*σ_f^3
    ```
  - 保持函数签名不变，确保 `gap_conditions` 与求解器无需改动。
- 测试：比较开启/关闭 `g1,g2` 时的质量解差异；在 T=0, μ=0 下检验无 NaN 并收敛。  
- 估时：1–2 天
- 验收条件：在典型参数下（config 中）运行 `solve`，在 T=0 时返回物理解且质量满足新公式（通过单元测试验证表达式）。

---

### Issue 3 — 实现 rPNJL Polyakov 势（式3.27–3.29）与 `safe_vandermonde`
- 摘要：实现带温度依赖 `b2(T)` 与 Vandermonde 项的 Polyakov 势，并在数值异常时发出警告（不自动修正）。
- 修改文件：
  - `src/pnjl/core/Thermodynamics.jl`（或 `Polyakov.jl`）——新增 `calculate_U_rpnjl(T,Φ,Φ̄)` 或在现有 `calculate_U` 中按 `model` 分支。  
  - 日志/工具模块（若无）新增 `safe_vandermonde(value; ε=1e-12)`，在 `value<=0` 时记录 `@warn` 并返回原值（不夹紧）。
- 主要改动说明：
  - 实现：
    - `b2(T) = a0 + a1*(T0/T)*exp(-a2*T/T0)`
    - `J(Φ,Φ̄) = 27/(24π^2) * [1 - 6ΦΦ̄ + 4(Φ^3+Φ̄^3) - 3(ΦΦ̄)^2]`
    - `U/T^4 = -b2(T)/2 ΦΦ̄ - b3/6(Φ^3+Φ̄^3) + b4/4 (ΦΦ̄)^2 - κ*log(J)`
  - 在 `log(J)` 处使用 `@warn` 当 `J<=0`（按用户决策，不修正）。
- 测试：验证 `calculate_U_rpnjl` 在一系列 (T,Φ,Φ̄) 下返回有限值；在 `J<=0` 时能发出警告日志；与旧 `calculate_U` 在 κ=0、g1=g2=0 时一致性检查。  
- 估时：1–2 天
- 验收条件：`calculate_U` 支持 rPNJL 形式并由 `model` 开关选用；日志在 `J<=0` 时记录警告且程序不崩溃。

---

## 后续
完成上述三个 issue 后，下一步将拆分并生成 Issue 4–7（接口重构、AD 验证、测试脚本与基准）。如需我现在将 Issue 1–3 转为仓库中的 issue 模板或直接生成 patch 草案，请回复“生成 issue”或“生成 patch”。
