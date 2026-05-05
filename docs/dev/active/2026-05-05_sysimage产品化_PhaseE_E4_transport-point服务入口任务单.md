# sysimage 产品化 Phase E：E4 transport-point 服务入口任务单

更新日期：2026-05-05

当前状态：已启动；承接 E3 收口后的下一批 service 扩展，实现 `transport-point` 作为 point-family 服务入口的最小落地。

> 目的：在不破坏 E3 已固化的 CLI/service 边界前提下，把 `transport-point` 从“稳定 API + warmup capability”推进到“现有 server 可调用的同步 point endpoint + discovery contract”。

---

## 1. 前置结论

- [x] E1 已确认 `transport-point` family 是第一批最值得补进 server 的同步 point 候选
- [x] E2 已把 `transport_point` 纳入 `point / service_core` warmup profile
- [x] E3 已固化：
  - `point-family` 可同时保留 CLI/API，但 service 为默认在线入口
  - `/api/modules` discovery 是 service 边界的机器可读真源

---

## 2. 本批目标

- [x] E4-1 为 `transport-point` 设计最小 HTTP 契约
- [x] E4-2 复用现有 `Models.solve_transport_from_equilibrium` 稳定 API，不再引入脚本耦合
- [x] E4-3 接入 `FullServerApp` 路由与 module discovery
- [x] E4-4 补 service contract / integration test

---

## 3. 范围与非范围

### 3.1 本批范围

- [x] `transport-point` request/response contract
- [x] `FullServerApp` handler / route / module registry 接入
- [x] `docs/api/data_contracts.md` 更新
- [x] integration contract / smoke test

### 3.2 非范围

- [ ] 本批不改造 `phase` 为 job endpoint
- [ ] 本批不改造 `transport-scan` 为同步接口
- [ ] 本批不重写 transport 物理解算逻辑
- [ ] 本批不处理 release / installer / app 分发

---

## 4. 初步边界

- [x] 入口类型：`sync + point`
- [x] discovery：纳入 `/api/modules`
- [x] 默认在线入口：service
- [ ] CLI 仍保留现有批处理入口，不以本批新增 transport-point CLI 为目标

---

## 5. 验收标准

- [x] `/api/modules` 可发现 `transport-point`
- [x] 存在稳定 `POST` point endpoint
- [x] 请求/响应字段、错误语义、单位约定已写入数据契约
- [x] integration 测试覆盖最小 happy path 与参数校验

---

## 6. DoD

- [x] `transport-point` 已成为 E3 边界内的真实 service module
- [x] 不与现有 CLI/job 路径冲突
- [ ] 后续若继续扩 `phase` job / artifact endpoint，可在此基础上独立推进
