# 使用指南（Web + API）

本指南提供当前仓库结构下的标准使用方式，适用于：

- 本地启动 API 与前端联调
- 基础功能验证
- 常见错误排查

## 1. 启动方式

### 方式 A（推荐）

```powershell
julia --project=. scripts/server/server_full.jl
```

### 方式 B（Windows 一键）

```powershell
.\scripts\server\start.bat
```

### 方式 C（仅 API）

```powershell
julia --project=. scripts/server/server.jl
```

## 2. 页面入口

- 主页面：`web/index.html`
- 简化测试：`web/simple_test.html`
- API 交互测试：`web/test_api.html`

## 3. 验证流程

1. 健康检查：访问 `http://localhost:8080/health`，应返回 `OK`
2. 打开 `web/index.html`，确认页面状态指示器显示在线
3. 使用默认参数发起一次 `/compute` 请求并观察结果
4. 切到“扫描任务中心”，可触发：
   - `pnjl-scan` 异步任务：`tmu/trho` + `scan/point`
   - `pnjl-gap` 同步单点：`POST /api/modules/pnjl-gap/run`
   - `transport-point` 同步单点：`POST /api/modules/transport-point/run`
   - `script-tasks` 异步脚本长任务：前端先展示任务用途、关键参数、输出与 preset，再通过 `POST /api/modules/script-tasks/jobs` 提交

## 4. 测试建议

```powershell
# 单元测试入口（默认 smoke）
julia --project=. tests/unit/runtests.jl

# 显式 smoke
$env:UNIT_PROFILE='smoke'
julia --project=. tests/unit/runtests.jl
```

## 5. 常见问题

- **端口占用**：`julia --project=. scripts/server/server_full.jl 8081`
- **前端显示离线**：先检查 `http://localhost:8080/health`
- **Three.js 网络加载失败**：参考 `web/THREEJS_LOCAL.md`
- **依赖未安装**：执行 `julia --project=. -e "using Pkg; Pkg.instantiate()"`

## 6. 状态口径说明

- Web/API 演示链路可用
- PNJL 求解与扫描链路可用
- 截面/弛豫时间链路已验证可用（建议结合对比报告持续回归，以 `README.md` 为准）

## 7. 输出目录

- 默认：`data/outputs/`
- 前端脚本长任务默认写入：`data/outputs/frontend_jobs/{job_id}`
- 前端脚本长任务日志默认写入：`data/outputs/logs/frontend_jobs/{job_id}.out.log` 与 `.err.log`
- `outputs/` 仅作历史兼容目录，不作为新流程默认落盘路径

## 8. 脚本长任务安全口径

- `smoke` 是默认 preset，适合先理解任务作用与产物形状。
- `canonical` 和 `custom` 属于重任务，前端必须勾选 `confirm_heavy` 才会提交。
- 前端脚本任务使用 `run_with_sysimage` 的 fallback policy：有兼容 sysimage 时复用，没有时回退到普通 Julia，不在点击时自动重建 sysimage。
- `custom_args` 按“每行一个 argv 参数”填写，不按 shell 字符串拆分。
- `scripts/dev/`、`scripts/analysis/`、`scripts/debug/`、`scripts/perf/` 默认只作为诊断/治理信息，不作为普通前端可点选执行入口。
- 若任务需要已有输入文件（例如 Mott 派生 CSV、plot modes、离线补点），前端只提供 `custom` 入口，用户需显式填写 `--in/--input` 与输出参数。
