# sysimage 产品化与 package / app / service 演进任务单

更新日期：2026-05-04

当前状态：Phase A 已拉起到 `docs/dev/active/2026-05-04_sysimage产品化_PhaseA_wrapper默认化任务单.md` 并完成首轮收口；Phase B 已拆成两批活动文档并完成首轮收口：`docs/dev/active/2026-05-04_sysimage产品化_PhaseB_metadata与release资产任务单.md`（B1-B4）与 `docs/dev/active/2026-05-04_sysimage产品化_PhaseB_release-workflow与资产打包任务单.md`（B5-B6）。

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
- [ ] C2 把稳定 CLI 执行路径进一步收敛到 `src/models` / `src/simulation` 暴露的稳定 API
  - 验收：脚本更多承担参数解析与入口职责
- [ ] C3 扩展 precompile capability 与 workload，使其绑定稳定 API 而非脚本偶然路径
  - 验收：precompile registry 与 stable API 对齐

### Phase D：app 化

- [ ] D1 设计 launcher 层
  - 验收：用户可通过固定 launcher 启动稳定能力
- [ ] D2 评估 `PackageCompiler.create_app(...)` 是否适合当前仓库
  - 验收：给出适合/不适合与约束
- [ ] D3 形成最小 app 分发 PoC
  - 验收：至少一个稳定入口可通过 app/launcher 形式运行

### Phase E：service 化

- [ ] E1 盘点哪些能力最适合长驻服务化
  - 验收：至少区分 scan / point / phase / orchestrator
- [ ] E2 设计服务预热策略
  - 验收：明确启动时预热哪些稳定路径
- [ ] E3 设计“CLI 入口”和“服务入口”的职责边界
  - 验收：避免双套契约长期漂移

---

## 4. 推荐顺序

- [x] Step 1：先做 Phase A wrapper 默认化
- [x] Step 2：再做 Phase B sysimage 分发
- [ ] Step 3：然后推进 Phase C package 化
- [ ] Step 4：在入口收敛后做 Phase D app 化
- [ ] Step 5：最后做 Phase E service 化

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
