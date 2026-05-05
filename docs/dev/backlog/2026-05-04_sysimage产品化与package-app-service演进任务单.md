# sysimage 产品化与 package / app / service 演进任务单

更新日期：2026-05-05

当前状态：Phase A 已完成首轮收口；Phase B 已完成首轮收口并具备 metadata / bootstrap / release workflow 基础；Phase C 已完成 package 化主路径收敛与首轮性能精修验证；Phase D 已完成 D1 launcher 设计、D2 `create_app(...)` 适配性评估与 D3 最小 launcher app PoC 主收口；Phase E 已完成 E1 服务化候选能力盘点、E2 服务预热策略首批收口与 E3 CLI/服务职责边界首批固化。本轮 sysimage 产品化主线已完成阶段性收尾，后续若继续做 `transport-point` endpoint 或 phase job/artifact 扩展，应另起新的 active 任务单承接。

> 目的：把当前已验证有效的 sysimage / precompile 冷启动优化，从“仓库内工程能力”推进到“用户默认入口、跨重启稳定复用、跨机器可分发、长期架构演进”。

---

## 1. 背景与目标

- 当前仓库已经完成：
  - 稳定 CLI 冷启动基线
  - `scan_pipeline_cli` precompile capability
  - 本地 sysimage 构建
  - `run_with_sysimage.ps1` wrapper 原型
- 已验证结论：
  - 冷启动时间大头来自启动/JIT
  - sysimage 对稳定 CLI 有显著改善，但还没有逼近热态
  - 若不统一入口，用户与零知识 agent 无法稳定吃到优化

下一阶段目标：

- [ ] 把 sysimage / wrapper 从“可选优化”提升为“稳定入口默认路径”
- [ ] 让 fresh clone / 新机器 / 零知识 agent 也有可操作的 sysimage 获取方式
- [ ] 为后续 package 化 / app 化 / service 化提供明确分阶段路线

---

## 2. 范围与非范围

### 2.1 本期范围

- [ ] wrapper 默认入口产品化
- [ ] 预构建 sysimage 资产分发方案
- [ ] launcher / bootstrap 脚本方案
- [ ] package / app / service 三阶段演进设计

### 2.2 非范围

- [ ] 不直接重构当前活跃物理主线
- [ ] 不在本期强推所有稳定脚本改造成 app
- [ ] 不在没有统一契约的前提下直接 service 化全部能力

---

## 3. 阶段拆分

### Phase A：wrapper 默认化

- [x] A1 盘点所有稳定 CLI，定义哪些必须默认经 `run_with_sysimage` 启动
  - 验收：形成白名单映射表
- [x] A2 为 Windows 用户入口补齐统一命令示例
  - 验收：README / guides / scripts 文档不再以裸 `julia --project=.` 作为稳定入口主示例
- [x] A3 评估是否补 `run_with_sysimage.sh`
  - 验收：明确 Linux/macOS 是否进入同一产品化路径
  - 2026-05-04：Phase A 已拆到 `docs/dev/active/2026-05-04_sysimage产品化_PhaseA_wrapper默认化任务单.md` 并完成首轮收口；若后续需要跨平台实机验证，可作为补充批次继续推进。

### Phase B：sysimage 获取与分发

- [x] B1 定义 sysimage metadata 契约
  - 验收：至少包含 `julia_version / git_commit / platform / generated_at / artifact_name`
- [x] B2 设计 GitHub Release 资产命名与平台矩阵
  - 验收：明确 win/linux/macos 与 x64/arm64 支持策略
  - 2026-05-04：Phase B 已拆到 `docs/dev/active/2026-05-04_sysimage产品化_PhaseB_metadata与release资产任务单.md`；metadata 契约与 release 资产命名 / 平台矩阵在该文档中持续维护。
- [x] B3 设计 bootstrap 下载脚本
  - 验收：fresh clone 可通过一条命令获取匹配 sysimage
- [x] B4 设计“版本不匹配时”的回退策略
  - 验收：明确 strict / fallback / rebuild 行为
  - 2026-05-04：已新增 `scripts/dev/bootstrap_sysimage.ps1/.sh`，并为 `run_with_sysimage.ps1/.sh` 固定 mismatch policy 契约；后续只需补 release workflow 即可落真实分发。
- [x] B5 实现 release 资产打包脚本
  - 验收：本地可基于 metadata 产出可上传的压缩资产与校验文件
- [x] B6 实现 GitHub Actions release workflow
  - 验收：平台矩阵构建结果可作为 workflow artifact / release asset 发布
  - 2026-05-04：已新增 `scripts/dev/package_sysimage_release.jl` 与 `.github/workflows/sysimage-release.yml`；本地已验证 Windows 资产打包输出。

### Phase C：package 化

- [x] C1 盘点当前稳定脚本中仍过厚的顶层逻辑
  - 验收：形成 thin-wrapper 候选列表
  - 2026-05-04：已拉起 `docs/dev/active/2026-05-04_sysimage产品化_PhaseC_package化_C1薄封装盘点任务单.md` 并完成首轮盘点；`run_gap_transport_scan.jl` 被识别为首要收敛对象。
- [x] C2 把稳定 CLI 执行路径进一步收敛到 `src/models` / `src/simulation` 暴露的稳定 API
  - 验收：脚本更多承担参数解析与入口职责
  - 2026-05-04：已拉起 `docs/dev/active/2026-05-04_sysimage产品化_PhaseC_package化_C2首批边界收敛任务单.md` 并完成首批收敛；`run_gap_transport_scan.jl` 的 CLI 与脚本 IO 层已抽离为 helper，但核心求解执行体仍保留在脚本中。
  - 2026-05-04：已继续拉起 `docs/dev/active/2026-05-04_sysimage产品化_PhaseC_package化_C2第二批provenance壳层收敛任务单.md` 并完成第二批收敛；`effective_config / summary / artifacts / sidecar` 汇总壳层已抽离为 helper。
  - 2026-05-04：已继续拉起 `docs/dev/active/2026-05-04_sysimage产品化_PhaseC_package化_C2第三批phase-equilibrium边界收敛任务单.md` 并完成第三批收敛；`phase tracker + equilibrium dispatch` 已形成 script-private helper 边界，但暂不升格为 `src/` 稳定 API。
  - 2026-05-04：已继续拉起 `docs/dev/active/2026-05-04_sysimage产品化_PhaseC_package化_C2第四批orchestration壳层收敛任务单.md` 并完成第四批收敛；`run_scan` 重复的单点执行链已抽离为 orchestration helper，且 transport 后处理已对齐 `Models.solve_transport_from_equilibrium`。
  - 2026-05-04：已继续拉起 `docs/dev/active/2026-05-04_sysimage产品化_PhaseC_package化_C2第五批scan-plan收口任务单.md` 并完成第五批收口；扫描计划与循环控制已独立成 helper，`run_gap_transport_scan.jl` 已基本退回到装配职责。
- [x] C3 扩展 precompile capability 与 workload，使其绑定稳定 API 而非脚本偶然路径
  - 验收：precompile registry 与 stable API 对齐
  - 2026-05-04：已拉起 `docs/dev/active/2026-05-04_sysimage产品化_PhaseC_package化_C3首批precompile-workload对齐任务单.md` 并完成首批对齐；`scan/core/full` profile 已新增 `transport_point_api` capability，开始命中稳定 transport point API。
  - 2026-05-04：已继续拉起 `docs/dev/active/2026-05-04_sysimage产品化_PhaseC_package化_C3第二批transport-point收益测量任务单.md` 并完成第二批测量；`transport_point_api` 对 `no-sys -> with-sys` 的 focused workload 有明确收益，但相对当前主 sysimage 未显示额外的 residual-trace 压降。
  - 2026-05-04：已在第二批测量中继续完成 focused residual trace diff；`main_with_sys` 与 `probe_with_sys` 的 focused residual 完全相同，说明仅新增 `transport_point_api` 未继续压缩 residual。
  - 2026-05-04：已进一步确认 `Constants_PNJL.ConfigLoader` / `Models.TransportWorkflow.ConfigLoader` 双实例是 `config_loader` 桶的主要根因；`TransportWorkflow` 复用 `Main.Constants_PNJL.load_config` 后，`with-sys` focused residual common 从 `97` 降到 `85`，`config_loader` 桶从 `30` 降到 `18`，下一步优先级应转向 `solver_ad`。
  - 2026-05-04：已继续拉起 `docs/dev/active/2026-05-04_sysimage产品化_PhaseC_package化_C3第三批solver-ad残余收敛任务单.md`；本批目标是把 `Models.solve` 内部 `ForwardDiff.Tag / DifferentiationInterface.prepare_jacobian / NLSolversBase.OnceDifferentiable` 相关 residual 从混合 transport focus 中单独拆出，并用更窄 probe 判断哪些签名还能继续被 sysimage 吃掉。
  - 2026-05-04：C3 第三批已完成首轮 solver-only probe；`build/solver_ad_probe/JuliaRelaxTimeSolverAdProbe.dll` 对 `Models.solve(FixedMu)` workload 可把 solver focus 从 `79` 压到 `42`，当前共同 residual 已收缩为 `ForwardDiff.Tag = 38`、`prepare_jacobian = 2`、`solve_kwcall = 2`，下一步应继续围绕 `ForwardDiff.Tag` 主导签名做更细 probe，而不是先扩 capability 集合。
  - 2026-05-04：C3 第三批已继续完成 residual-only probe；`Conditions.build_residual! + ForwardDiff.jacobian` workload 可把更窄 focus 从 `65` 压到 `33`，新增吃掉的 `32` 条主要是 `ForwardDiff misc / partials`，共同 residual 则进一步集中到 `Models.Conditions.var"#gap_conditions..."` 的 nested-tag 主干。基于该证据，`src/models/precompile/registry.jl` 已新增 `solver_residual_ad` capability，并接入 `scan/core/full/all` profile，下一步应验证其对主 sysimage focused residual 的实际收益。
  - 2026-05-04：主 sysimage 复测已完成；`solver_residual_ad` 接入 profile 后，`solver_ad` focused residual 与 `transport_point` focused residual 均未继续下降，且新旧 `with-sys` trace 逐行一致。这说明当前通过“追加 workload / capability”还能吃掉的稳定签名已基本见顶；后续若继续推进 C3，应转向 `residual!#140 / _gap_conditions_with_model` 的 closure/tag 形状治理，而不是继续叠加 probe。
  - 2026-05-04：已继续拉起 `docs/dev/active/2026-05-04_sysimage产品化_PhaseC_package化_C3第四批closure-tag形状治理任务单.md` 并完成首轮收口；`GapResidualInplace / GapOmegaFn / GapPressureFn` 已替换匿名 solver/conditions closure family。solver-focused probe 的 with-sys focus `42 -> 41`，且主 sysimage 真实路径也同步下降：`solver_ad main with-sys 42 -> 41`、`transport_point main with-sys residual common 87 -> 86`。结论是 closure/tag 形状治理有效，但已进入“单批边际压降”阶段；C3 后续若继续推进，应围绕剩余 `config_loader` / `solver_ad` / transport helper 桶精修，而不是继续扩 capability 集。
  - 2026-05-05：已继续拉起 `docs/dev/active/2026-05-05_sysimage产品化_PhaseC_package化_C3第五批config-defaults与冷启动指标任务单.md`。本批将 `TransportWorkflow` 默认 TOML 配置收敛为单次 `load_config` + 双字段复用，并新增 `scripts/perf/relaxtime/compare_sysimage_summary_walltime.jl`，正式把 `wall_ms` 升为 C3 的主展示指标。当前证据显示：`transport_point` 主 workload 的 with-sys focused residual 完全未变，因此本批代码整理尚未命中主路径 residual family；后续 C3 更应直接针对 `transport_api` 桶里的 kwargs/generator family，而不是继续处理未命中的 config-defaults 分支。
  - 2026-05-05：已继续拉起 `docs/dev/active/2026-05-05_sysimage产品化_PhaseC_package化_C3第六批transport-kwargs形状收敛任务单.md`。本批将 `TransportWorkflow.solve_transport_from_equilibrium(...)` 内 `transport_kwargs` 的匿名 `Filter/Generator/merge` 链收敛为显式 materialize helper；主 sysimage 的 `transport_point` with-sys focused residual `86 -> 84`，且 transport-point focus lens `20 -> 18`，说明本批首次直接命中了 `transport_api` 桶的真实残余家族。wall time 单次结果仍有宿主波动，因此继续采用“`wall_ms` 主展示、with-sys trace diff 归因”的口径更稳妥。
  - 2026-05-05：阶段判断更新：C3 主目标已达成，即“稳定 API 已命中主 workload，后续优化进入边际残余治理”；后续若继续推进 C3，应视为可选 perf refinement，而非当前主线。

### Phase D：app 化（当前主线）

- [x] D1 设计 launcher 层
  - 验收：用户可通过固定 launcher 启动稳定能力
  - 2026-05-05：已拉起 `docs/dev/active/2026-05-05_sysimage产品化_PhaseD_D1_launcher层设计任务单.md` 并完成首轮设计收口；launcher 被明确限定为“统一稳定入口层”，首批白名单以 `phase / unified-scan / transport-scan / relaxtime-orchestrator / suscept / server` 为主，`transport-point` 暂不提升为 launcher 顶层命令。
- [x] D2 评估 `PackageCompiler.create_app(...)` 是否适合当前仓库
  - 验收：给出适合/不适合与约束
  - 2026-05-05：已拉起 `docs/dev/active/2026-05-05_sysimage产品化_PhaseD_D2_create_app适配性评估任务单.md` 并完成首轮评估；结论为“当前根仓库不适合直接 `create_app(...)`，但薄 launcher package 壳适合作为 app 承载形态”，因此 D3 应以最小 launcher package app PoC 为目标，而不是直接对根仓库做 app 化。
- [x] D3 形成最小 app 分发 PoC
  - 验收：至少一个稳定入口可通过 app/launcher 形式运行
  - 2026-05-05：已拉起 `docs/dev/active/2026-05-05_sysimage产品化_PhaseD_D3最小launcher-app-PoC任务单.md`；本批严格收敛为“薄 launcher package 壳 + `phase` 子命令 + `create_app(...)` 本地可运行验证”，并显式排除 `transport-scan` / `server` / installer / auto-update。
  - 2026-05-05：D3 首批已完成 launcher package 壳、`phase_cli_support` 复用层与本地解释器链路验证；`create_app(...)` 构建脚本已落地，且 `build/launcher-app/` bundle 目录能够开始产出，但自定义 launcher 可执行文件尚未在当前验证窗口内落盘，因此 D3 仍保持进行中，下一批应专门收口 app bundle 最终产物与构建耗时问题。
  - 2026-05-05：D3 第二批已拉起 `docs/dev/active/2026-05-05_sysimage产品化_PhaseD_D3第二批app产物收口任务单.md` 并完成首轮诊断；已确认前序 bundle 步骤正常，最初的结构性失败来自 launcher package 缺少 `Pkg` 依赖与同步 manifest，现已修复。当前真实停留点是 `create_sysimage(...)` 的长耗时编译阶段，尚未走到 `create_executable_from_sysimg(...)`，因此后续应把关注点放在本地 PoC 构建时长/参数，而不是怀疑 `executables` 命名契约失效。
  - 2026-05-05：D3 第三批已把 `scripts/dev/build_launcher_app.jl` 从 `create_app(...)` 黑盒调用改为“手工 bundle + repo sysimage 作为 base sysimage 的极薄 launcher sysimage”路线；当前该路线已能推进到 `build/launcher-app/lib/julia/` 目标路径创建，但在现有验证窗口内仍停留在 sysimage 长耗时编译阶段，因此下一批应把重点从‘结构是否正确’完全切到‘完成条件 / wall time / 是否需要更长窗口或更轻编译参数’。
  - 2026-05-05：D3 第四批已拉起 `docs/dev/active/2026-05-05_sysimage产品化_PhaseD_D3第四批launcher-sysimage完成验证任务单.md` 并完成首轮轻量参数验证；即使切到 `base_sysimage + include_transitive_dependencies=false + -O0 --compile=min --strip-metadata --strip-ir`，10 分钟窗口内仍停留在 `create sysimage start` 之后，未拿到 `sys.dll` / `jrt-launcher(.exe)` 证据。当前可认为 D3 的剩余阻塞已明确是“本机短窗口内 launcher sysimage 编译 wall time 不可接受”，而不是结构或命名契约问题。
  - 2026-05-05：阶段决策更新：接受 D3 以“仓库内 launcher PoC 已成立、最终可执行产物长时验证另列”为收口条件。即，D3 主目标已完成，后续若仍需拿 `sys.dll` / `jrt-launcher(.exe)` 证据，应单独开长时离线构建批次，而不再阻塞 Phase E 主线。

### Phase E：service 化

- [x] E1 盘点哪些能力最适合长驻服务化
  - 验收：至少区分 scan / point / phase / orchestrator
  - 2026-05-05：已拉起 `docs/dev/active/2026-05-05_sysimage产品化_PhaseE_E1服务化候选能力盘点任务单.md`，用于承接 D3 收口后的主线切换；本批先做服务化候选白名单 / 暂缓名单与边界判断，不直接实现服务进程。
  - 2026-05-05：E1 已完成首轮盘点收口；结论为“现有 `server_full` / `FullServerApp` 直接作为 service 壳继续演进，point-family 为第一批同步服务化白名单，phase/scan/orchestrator 统一按异步 job/worker 语义推进”。因此 E2 应优先设计 `solve_pnjl_point`、`solve_transport_from_equilibrium`、`run_phase_pipeline` 的预热路径，而不是先扩更多顶层 endpoint。
- [x] E2 设计服务预热策略
  - 验收：明确启动时预热哪些稳定路径
  - 2026-05-05：已拉起 `docs/dev/active/2026-05-05_sysimage产品化_PhaseE_E2服务预热策略任务单.md`，并完成首批收口；当前已形成 `none / point / service_core` warmup profile 契约，`localhost / staging / remote` 的默认映射，以及接入 `ServerLauncher.run_full_server(...)` 的最小 warmup hook。后续 E3 可直接在该契约之上继续落 endpoint 与职责边界。
- [x] E3 设计“CLI 入口”和“服务入口”的职责边界
  - 验收：避免双套契约长期漂移
  - 2026-05-05：已拉起 `docs/dev/active/2026-05-05_sysimage产品化_PhaseE_E3_CLI与服务职责边界任务单.md`，并完成首批收口；当前已把 `point` 与 `job` 的 service 表面边界固化到 `/api/modules` discovery 契约，明确 `point-family` 可同时保留 CLI/API 但 service 为默认在线入口，而 `scan / phase / orchestrator` 继续以 CLI 为权威执行入口、service 仅提供 job 壳。

---

## 4. 推荐顺序

- [x] Step 1：先做 Phase A wrapper 默认化
- [x] Step 2：再做 Phase B sysimage 分发
- [x] Step 3：然后推进 Phase C package 化
- [x] Step 4：在入口收敛后做 Phase D app 化
- [ ] Step 5：最后做 Phase E service 化

### 4.1 离开 C3 的门槛

- [x] 主 workload 已命中稳定 API，而不是脚本偶然路径
- [x] 继续追加 capability / workload 已基本见顶
- [x] 后续 residual 收敛进入边际区间，需要 trace diff 才能证明命中效果
- [x] 单次冷启动 wall time 已可作为展示指标，但批间强结论需要重复测量/中位数口径

结论：

- [x] 满足以上条件后，C3 不再作为主线阶段
- [x] 后续未关闭的 C3 工作应转入 perf refinement 支线，仅在有明确 KPI 或新热点出现时再继续

---

## 5. 验收标准

- [ ] 用户和零知识 agent 能明确知道默认稳定入口应走 wrapper
- [ ] 同机跨重启可稳定复用 sysimage
- [ ] fresh clone 有清晰的 sysimage 获取路径
- [ ] stable CLI 不再过度依赖脚本顶层逻辑
- [ ] package / app / service 三阶段边界清晰，不混成一次性大重构

---

## 6. 风险与约束

- [ ] 风险 R1：提前进入 service 化，导致契约与入口层混乱
  - 约束：必须先做 wrapper / package 化
- [ ] 风险 R2：sysimage 资产与 Julia 版本、平台不匹配
  - 约束：metadata + strict/fallback 策略必须先定义
- [ ] 风险 R3：把性能优化路径做成平台私有技巧
  - 约束：至少保留跨平台兼容方案说明
- [ ] 风险 R4：文档口径与真实默认入口不一致
  - 约束：每次入口策略变更都同步 README / QUICKSTART / scripts guide

---

## 7. DoD

- [ ] 已形成下一阶段可执行任务分解
- [ ] sysimage 产品化、package 化、app 化、service 化的先后顺序明确
- [ ] 可在当前任务单收尾后直接拉起执行

### 7.1 当前阶段决策

- [x] Phase C 主线收口，保留 C3 为可选性能精修支线
- [x] Phase D 主线收口，接受 D3 为“仓库内 launcher PoC 已成立，长时最终产物验证另列”
- [x] Phase E 已完成 E1：服务化候选白名单 / 暂缓名单与同步-异步边界盘点
- [x] Phase E 已完成 E2：服务预热 profile / 默认映射 / 最小启动 hook 首批收口
- [x] Phase E 已完成 E3：CLI 入口与服务入口职责边界首批固化
- [x] 当前总任务单已完成本轮主线收尾；后续扩展改由新任务单承接
 - [x] 后续扩展已起新单：`docs/dev/active/2026-05-05_sysimage产品化_PhaseE_E4_transport-point服务入口任务单.md`
