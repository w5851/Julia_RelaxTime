# 项目审查：修复与待审清单（2026-11-17）

## 范围与结论

审查起点：`main` / `cbb62cb3`（已包含 #349、#350、#351）。

对起点的 8,566 个受控文件做了目录、入口、语法、依赖和引用盘点；按调用关系检查了 Models 工作流、扫描任务服务、前端存储、测试入口、CI 和当前文档。重点阅读与执行验证集中在实际异常及其调用方。

| 范围 | 覆盖情况 |
| --- | --- |
| 核心与工具 | `src` 218、`scripts` 635、`tests` 697、`web` 16、`app` 4、`benchmark` 41、`config` 80 个受控文件纳入盘点 |
| Julia | 1,200 个受控 `.jl` 全量静态解析；修复一个现行 CLI，剩余一个错误位于保留的历史产物副本，见 R8 |
| Python / TOML | 485 个 `.py`、270 个 `.toml` 静态解析通过；运行 Python unit 全套及修复后的定向测试 |
| 文档 | 检查当前导航、API、guides、active/backlog 中 168 份 Markdown 的普通链接；另查反引号中的路径，修复确认为过时的源码/命令引用 |
| 数据与历史 | `data` 4,993、`docs` 1,783 个文件纳入盘点；通过现有产物测试检查相关来源与验收语义，没有逐件重新验收科研数据 |

这次覆盖不等于逐行人工读完全部历史文件，也不证明所有缺陷均已排除。没有修改物理公式、基线、数值容差或正式图件。本地按受影响范围验证，完整 Nightly 以远端结果为准。缓存对照仅代表下述测试负载。

开始时已有两份 SOP 的本地修改；次轮复核沿用其分析建议，并补充适用范围，见 R7。外部 worktree 保持原样。

## 已直接修复

| 编号 | 问题与修复 | 主要文件 / 验证 |
| --- | --- | --- |
| F1 | crossover meson-density CLI 的 `@__FILE__ && ...` 写法无法解析；改为正常 `if` 入口 | [run_crossover_meson_density_scan.jl](../../scripts/relaxtime/run_crossover_meson_density_scan.jl)；实际 `--help` 退出 0 |
| F2 | Julia 脚本入口检查只核对路径，能放过语法错误；加入只解析、不执行的语法检查，并覆盖坏 macro guard、未完成语法和 include 引用 | [check_script_entrypoints.jl](../../scripts/dev/check_script_entrypoints.jl)、[fixture 测试](../../tests/unit/config/test_script_entrypoints.jl)；370 个现行脚本通过，CI 接入 fixture |
| F3 | unit full 只遍历一层并特殊处理 solver，漏掉 derivatives/workflow 下 7 个测试；改为逐层递归，移除 solver 重复入口 | [unit runner](../../tests/unit/runtests.jl)；7 个恢复的测试文件及入口检查 fixture 共 189/189；临时嵌套 fixture 证明顺序、排除规则和无重复 |
| F4 | `PipelineSpec.stages` 中重复 id 被 Set 校验掩盖，同一阶段可能执行两次；执行前拒绝重复 id | [PipelineRunner.jl](../../src/models/workflow_engine/PipelineRunner.jl)、[行为测试](../../tests/unit/models/workflow/test_pipeline_runner_behavior.jl)；验证零阶段副作用、无 manifest |
| F5 | 禁用 localStorage 时 getter 在 try 外抛错；坏 history/template 记录可能使后续代码崩溃；默认模板对象可被调用方污染 | [task_store.js](../../web/js/task_store.js)；getter 纳入异常处理、过滤坏记录、默认模板深拷贝、历史读写均限 50 条，前端测试通过 |
| F6 | scan 被取消后失败仍可重试；超时完成可能被标成成功；timeout 缺少结束事件 | [pnjl_scan_jobs.jl](../../src/simulation/fullserver/pnjl_scan_jobs.jl)；增加尝试边界状态/期限检查及 TIMEOUT 事件，worker 行为测试通过。扫描停止边界见 R1 |
| F7 | 输出路径校验在 Linux 上错误忽略大小写，未解析已有 symlink，且非法路径先登记再异步失败 | 同上；按平台比较路径、解析已有祖先、拒绝目录；显式非法路径登记前返回 `400 INVALID_REQUEST`。Windows 常规边界与 HTTP 响应通过，POSIX 分支待 Linux CI |
| F8 | 普通 Node 测试自动访问服务；联调默认最多提交 60 × 2,304 个扫描点且不清理；错误路径断言仍期待异步失败 | [api.scan.integration.test.mjs](../../web/js/api.scan.integration.test.mjs)；联调和队列压力分别显式启用，指定测试服务，在 finally 取消本测试创建的任务，兼容结果已快速完成的情况；移除失效的终态轮询 helper |
| F9 | 冻结 phase-map 来源 JSON 在 Windows checkout 时被转换成 CRLF，导致既有 SHA256 校验失败 | [.gitattributes](../../.gitattributes) 对精确文件加 `-text`；工作文件恢复为 Git 原始 blob，验收哈希与图件未改；来源校验复验通过 |
| F10 | 6 个分析脚本直接依赖 pandas，但安装清单没有它；Python 语法 gate 只覆盖零散旧清单 | [requirements.txt](../../requirements.txt) 补 pandas；[Python syntax gate](../../scripts/dev/check_python_script_syntax.py) 覆盖整个 scripts，213 个脚本通过；CI 触发路径包含所有脚本 Python 文件 |
| F11 | 一个 unit 测试硬编码个人 `D:\Temp` 历史目录，当前环境直接 skip | [raw-curve archive tests](../../tests/unit/python/test_pnjl_raw_curve_archive.py) 改为两个临时 fixture，验证有/无 production 目录但缺少来源 manifest 时均不可恢复；无须本机历史目录 |
| F12 | solver/workflow 源码路径、已删除的 TCS 测试路径、相图命令和“latest”性能数字过时 | 同步相关 API/指南；历史 timing 明确标为历史运行。新相图命令实际生成 3 张临时 PNG，并说明 accepted 层仍可包含未认证行 |

F9 涉及的文件是 `data/outputs/results/pnjl/phase_diagram/figure4_phase_diagram_prod_v1/phase_reference_source_manifest.json`。
恢复后的 working tree、index、HEAD 三者字节一致，SHA256 均为 `58f8b59d6e63deaa4d87fed30b0427caa930d966f9e16dc66b68d389cc47476f`；没有用更新预期哈希的方式消除失败。

## 审阅后的处理（2026-11-18）

### R1 · 运行中的 PNJL 扫描在下一处边界停止

取消、超时共用独立的停止回调。当前单点求解可以完成，扫描在下一点、延拓步或写入前响应。已刷新的 CSV 行按原样保留，可按同一网格续扫；未写入的计算结果可以丢弃。超时从取得 worker 起用单调时钟计时，并阻止后续重试。状态与取消响应提供 `worker_active`；它为 `true` 时仍占用名额。

[真实小网格测试](../../tests/integration/pnjl/test_scan_stop_boundaries.jl)验证 Tμ/Tρ 各写一行后停止，续扫只补缺行，PALC 与自然延拓均能响应停止。[worker 测试](../../tests/unit/simulation/test_pnjl_scan_job_execution.jl)覆盖取消、超时、无重试和释放名额。

### R2 · 取消脚本任务并等待进程退出

[启动门](../../src/simulation/fullserver/script_task_process.jl)在运行原有 wrapper 前纳管进程：Windows 使用 [Job Object](https://learn.microsoft.com/en-us/windows/win32/procthread/job-objects)，POSIX 经 [libuv detached 启动](https://docs.libuv.org/en/v1.x/process.html)建立独立进程组。启动和登记与取消使用同一把锁。取消先阻止尚未启动的命令，或终止已纳管的任务进程；实际退出前继续占用名额，结果返回 `409 JOB_NOT_FINISHED`。退出后可查看日志和已留文件；取消留下的文件可能不完整。脚本任务没有自动超时策略。

[Windows 服务测试](../../tests/unit/simulation/test_fullserver_script_task_jobs.jl)通过真实 wrapper 启动三层 Julia，再逐个验证进程退出、结果开放和参数中的空格得到保留。Ubuntu WSL 的 10 项探针验证正常退出码、日志、组内后代退出、SIGKILL 退出码和组外进程不受影响。POSIX 进程组覆盖继承该组的子进程。

### R3 · 比较后采用 Julia Task 独立缓存

基于五个既有截面固定点测试串行、四 Task 同背景和四 Task 不同背景。五次热运行取中位数，单位为秒；现有无锁共享字典只跑串行。短锁原型仅锁字典/统计，数值求值放在锁外。

| 热缓存负载 | 原无锁字典 | 共享短锁原型 | Task 独立原型 |
| --- | ---: | ---: | ---: |
| 单 Task，2,000 个截面 | 0.125 | 0.304 | 0.207 |
| 四 Task，同背景，合计 8,000 个截面 | — | 1.315 | 0.433 |
| 四 Task，不同背景，合计 8,000 个截面 | — | 1.173 | 0.360 |

同背景热缓存下，两个原型均命中超过 99.8%；独立缓存累计 640 项，共享缓存 160 项。同背景和不同背景的冷运行都是毫秒量级，五次波动区间重叠。各负载的截面校验和一致。Task 独立缓存避免共享锁竞争，也避免一个计算的 `reset_cache!` 清掉另一任务的条目。正式实现比对照原型更简单，保留量化精度和公开签名；测试负载的分配量与原串行版相同。正式实现的墙钟时间受同机并行验证影响，本轮不据此承诺加速比。

[缓存](../../src/relaxtime/PolarizationCache.jl)的统计和重置现在只作用于调用 Task。两处使用旧字典注入候选极化值的诊断脚本使用临时独立缓存，退出后恢复调用方缓存和统计。该结论仅覆盖此缓存；其他输运共享状态不在此项并行验证范围内。

### R4 · 保留成功判据与标签

[现有列](../api/models/scans/Algorithms.md)足以区分已落盘行上留存的求解器结果、残差、物理有效性和扫描接受；`acceptance_reason` 只显示优先级最高的主理由。例如物理无效与残差超限可以同时发生，摘要只显示物理无效，组合各列仍可读出两种证据。未知历史信息和未选候选并未被记录，标签不能重建它们。

“收紧”指改变扫描成功条件，例如要求有限残差和物理有效性都通过，因而可能改变续用行和相图下游输入。减少或合并标签只影响报告粒度。保留 #351 的口径，包括 Tμ 已有物理检查及 Tρ 各路线自己的筛选。

### R5 · 清理重复历史重跑，Nightly 只验证关键产物合同

修复后全 Python unit 原有 666 项、80 子断言通过，595.44 秒。最慢三项共 326.81 秒，重复生成已完成迁移的历史 v2/v3 相图审计，其中两项执行同一个 builder。删除这三项，保留生成脚本、冻结审计表和现行迁移/来源验证。删减后全套 **663 项、80 子断言通过，250.91 秒**。这是同机两次观测，不是 GitHub runner 计时。

[Nightly Python job](../../.github/workflows/nightly-full-regression.yml)只执行原始曲线归档和相图交付两份测试，安装 Poppler 以实际渲染并验证 PDF。覆盖归档恢复、篡改拒绝、冻结产物字节和交付负例。本地以固定输入复验 11 项通过，72.13 秒。一次复验期间两个阶段复制到的图件 SOP 哈希不同，交付校验按预期拒绝了输出。全量 Python suite 继续供本地选择。

### R6 · 移除只匹配源码字符串的介子扫描测试

删去 `test_meson_scan_workflow_entry.jl` 及其 selector，共 33 条字符串断言。保留现行脚本语法检查、目录规则和工作流行为测试。不用等价的字符串断言替换原有测试。

### R7 · 保留可解释的物理对照，明确适用路线

[图件 SOP](../guides/sop/workflows/figure_production.md)的无质量、固定背景 RS 参考按研究问题选用。核对有效温度展开后保留原建议，字形例外回到版式条款。[介子密度 SOP](../guides/sop/workflows/meson_density.md)保留同背景的 bound/pair/Landau 分项和原始失败状态；明确通用 phase-shift 路线的 `dω/(2π)` 历史默认仍是比值口径，独立 charged GBU 后端的默认测度是 `dω/π`。生产资格由具体路线验收确定。

### R8 · 给历史来源副本标明现行入口

历史产物 `freezeout_envelope_20261018_selected_design` 中的脚本副本有旧 macro guard 语法错误，仍按历史字节保留。工作树、index 和 HEAD 的副本逐字节一致，SHA256 为 `80a1045177404ca18d1d7681cd3e3e50890bb3307e391896f8b99cadf4182d98`。[指南](../guides/scripts/meson_freezeout_comparison.md)已链接到[现行分析脚本](../../scripts/analysis/relaxtime/run_freezeout_envelope_refinement.jl)，并说明精确回放历史批次需恢复当时的代码、目录布局与输入。

## 验证结果与边界

| 检查 | 结果 |
| --- | --- |
| Julia unit smoke | 200/200，4m46s |
| 恢复的 7 个嵌套单元测试 + script-entrypoint fixture | 189/189，5m54s |
| PNJL jobs、worker、HTTP utils 与脚本任务定向单测 | 最终组合 251/251，11m08s |
| 下一处扫描边界停止与续扫 | Tμ/Tρ、两种延拓和 HTTP 的定向 integration 50/50；最终排队取消和超时 HTTP fixture 单独复验 11/11；确定性运行中取消/超时及名额检查 19/19 |
| 脚本任务进程隔离 | Windows 服务定向测试 94/94；Ubuntu WSL 进程组探针 10/10 |
| Task 缓存、诊断脚本及输运采样定向单测 | 213/213，Julia 4 线程 |
| 截面和输运固定点回归 | 132/132；未调整基线或容差 |
| 输运模块加载顺序与工作流冒烟 | 三种加载顺序的子进程分别通过 28/28、41/41、41/41；integration 主进程 90/90，22m53s |
| unit 递归发现探针 | 通过；实际调用 runner，验证本层/子层顺序、performance/helper 排除和不重复 |
| Python unit 全套，修复前 | 663 passed、1 failed、1 skipped，80 subtests passed，691 秒；失败为 F9，跳过为 F11 |
| Python 初次修复后定向复验 | raw-curve archive 整个测试文件 + 原失败的来源校验，共 7 passed |
| Python 全套复验 | 清理前 666 项通过（595.44 秒）；清理后 663 项通过（250.91 秒）；均有 80 子断言通过 |
| Nightly Python 选择器 | Windows、Ubuntu WSL 和 GitHub Actions 均 11/11，含真实 PNG/PDF 和拒绝用例 |
| Node 默认全套 | 4 passed、1 skipped（显式保留的 live integration），无需运行服务 |
| 真 HTTP socket + Node 显式联调 | 健康、静态页面、400 路径拒绝、结果已完成/未完成、队列 429、测试任务清理通过；Julia 探针 10/10。实际 job/routing 实现配合确定性数值 stub，未进行真实大扫描压力测试；临时服务已关闭 |
| 现行 Julia script-entrypoint gate | 370 个脚本通过 |
| Python scripts syntax gate | 213 个脚本通过 |
| 全量受控 Julia 静态解析 | 1,200 个文件，仅剩 R8 的历史副本失败；新增测试另经实际执行 |
| 相图文档示例 | accepted 层、ξ=0 生成 3 张临时 PNG；日志保留 `uncertified_rows` 和 `runtime_consumption=false` |
| 治理 | task ledger、docs consistency、active docs、script entrypoints、models entry、solver leakage、relaxtime script governance、PNJL migration、dependency policy、agent instructions、model profile matrix、data output guard、precompile profile coverage、skill governance 共 14 项通过 |
| Git | `git diff --check` 通过；两份冻结来源文件与原始 blob 一致 |

本机临时日志和探针在 `tmp/project_audit_20261117/`，不作为新的长期科研证据包。关键修复有受控行为测试；本报告保留结论与必要的复验入口。

### 主要复验命令

```powershell
julia --project=. scripts/dev/check_script_entrypoints.jl
python -X utf8 scripts/dev/check_python_script_syntax.py
node --test web/js/*.test.mjs
python -X utf8 -m pytest -q tests/unit/python/test_pnjl_raw_curve_archive.py tests/unit/python/test_phase_diagram_prd.py
$env:UNIT_FILES = "simulation/test_pnjl_scan_jobs.jl,simulation/test_pnjl_scan_job_execution.jl,simulation/test_fullserver_script_task_jobs.jl,relaxtime/test_polarization_cache.jl,relaxtime/test_cross_section_loop_probe.jl"
julia --threads=4 --project=. tests/unit/runtests.jl
Remove-Item Env:UNIT_FILES
$env:INTEGRATION_FILES = "pnjl/test_scan_stop_boundaries.jl;relaxtime/test_pnjl_scan_cancel_timeout.jl"
julia --project=. tests/integration/runtests.jl
Remove-Item Env:INTEGRATION_FILES
$env:REGRESSION_FILES = "relaxtime/test_total_cross_section_fixedpoint_regression.jl;relaxtime/test_transport_fixedpoint_regression.jl"
julia --project=. tests/regression/runtests.jl
Remove-Item Env:REGRESSION_FILES
```

完整 live integration 的显式启用方式见 [command reference](agent_command_reference.md#frontend-tests)。定向验证覆盖本次修改触及的边界。完整 CI 状态以修复提交对应 GitHub Actions run 为准，未做全参数域数值重算。

## 提交与复验期间发现的 CI 失败（2026-11-18）

- [旧版 Nightly](https://github.com/w5851/Julia_RelaxTime/actions/runs/44059054844) 的 unit 和 Models/PNJL integration 都报同名 `SolverConvergence` 无法派发。相图 CLI 在进程内 `include` 时重新定义了 `Models`，两份已退役的 `run_tmu_scan.jl` 测试也有同样的问题，先加载的平衡态 facade 仍持有旧模块。现行相图 CLI 复用已加载模块并校验模块对象一致；删除两份旧测试，保留已有的退役与统一 CLI 行为覆盖。本地复现原错误，修复后同一模块中的平衡态计算正常，相关 unit 60/60 通过。
- [旧版预编译检查](https://github.com/w5851/Julia_RelaxTime/actions/runs/44069204694) 成功构建 sysimage，同一组测试在两种模式下均为 59/59。无 sysimage 1047 行、有 sysimage 401 行，只有覆盖差值 656 超出原上限 650。差值可能在覆盖改善时增大。改为约束两种模式的实际编译记录，分别最多 1100 和 450 行；保留总差值、重点差值和重点明细作为诊断。用 [11 月 12 日成功 run](https://github.com/w5851/Julia_RelaxTime/actions/runs/43196387013)、[11 月 13 日失败 run](https://github.com/w5851/Julia_RelaxTime/actions/runs/43404316698) 和上述最新失败数据复验，连同改善覆盖及超限负例共 7 项通过。[修复后远端复验](https://github.com/w5851/Julia_RelaxTime/actions/runs/44245187683) 构建成功，两种模式仍各通过 59/59，1047 和 401 行均在预算内。
- [首次 Nightly Python 检查](https://github.com/w5851/Julia_RelaxTime/actions/runs/44245187455/job/134602319216) 在 Ubuntu 上无法匹配冻结清单中的 Windows 相对路径。规范化读取后的路径，再按完整的相对路径核对记录、文件大小和 SHA256。真实渲染还复现了 Matplotlib 3.11 增加数学字形字段的问题；验收现只读取字体、字号和字符。Ubuntu WSL 配合 Matplotlib 3.11.2 和 Poppler 复验 11/11（26.12 秒）；Windows 绘图相关复验 69/69。[修复后的远端 Python job](https://github.com/w5851/Julia_RelaxTime/actions/runs/44247186578/job/134608866249) 11/11（测试耗时 4.80 秒），含真实 PNG/PDF 和拒绝用例。冻结来源与已验收图件的字节保持原样。
- [随后 Nightly 的 Julia unit](https://github.com/w5851/Julia_RelaxTime/actions/runs/44245187455/job/134602319682) 报告台账仍引用已删除的历史重跑测试，已指向仍验证 accepted 默认入口的测试。mixed-meson 诊断脚本也在测试进程中重新定义了 `Models`，导致后续介子热力学无法派发；本地复现后改为复用统一模块，并让诊断脚本的测试直接验证加载前后的模块对象一致。在同一进程中按顺序复验台账、平衡态入口、诊断脚本和介子热力学，共 284/284；独立运行诊断脚本 `--help` 和台账校验器均通过。
- [`b463b0d6` 的 Nightly unit](https://github.com/w5851/Julia_RelaxTime/actions/runs/44248993393/job/134614781045) 已通过台账校验，但仍有一处 `SolverConvergence` 派发错误。全 RS 自检加载联合验证及背景尾部验证脚本，最终由 [run_transport_background_validation.jl](../../scripts/analysis/relaxtime/run_transport_background_validation.jl) 再次定义 `Models`。本地按此顺序复现：全 RS 自检 10 项通过，但加载前后模块对象不同，紧接的平衡态求解失败。背景验证脚本改为复用 `Main.Models`，现有全 RS 测试补查这一行为。同一进程复验平衡态入口、全 RS 自检、mixed-meson 诊断及介子热力学，共 85/85；现行脚本入口检查 370 个文件通过。数值计算与容差沿用原实现。
