# 快速启动指南

本指南对应当前仓库结构（`scripts/server/`、`tests/unit/`、`docs/reference/`）。

说明：根项目运行兼容范围仍以仓库根 `Project.toml` 为准（当前 `julia = "1.10"`）；但 sysimage / launcher 的预构建与发布流程统一固定在 CI 基线 Julia `1.12.5`，以保证 `Manifest.toml`、预编译产物和 release 资产口径一致。

## 1) 初始化环境

```powershell
julia --project=. -e "using Pkg; Pkg.instantiate(); Pkg.precompile()"
```

如需在本机稳定复用冷启动优化，建议后续稳定 CLI 统一通过：

```powershell
powershell -ExecutionPolicy Bypass -File scripts/dev/run_with_sysimage.ps1 <script> [args...]
```

或在 Linux / macOS 环境下使用：

```bash
sh scripts/dev/run_with_sysimage.sh <script> [args...]
```

如需在 fresh clone / 新机器上直接获取匹配的预构建 sysimage，可先执行：

```powershell
powershell -ExecutionPolicy Bypass -File scripts/dev/bootstrap_sysimage.ps1
```

或：

```bash
sh scripts/dev/bootstrap_sysimage.sh
```

wrapper 默认采用 `fallback` 策略；如需强制要求兼容 sysimage，可改用 `-MismatchPolicy strict` 或 `--mismatch-policy=strict`。

## 2) 启动服务（API + 前端静态资源）

```powershell
powershell -ExecutionPolicy Bypass -File scripts/dev/run_with_sysimage.ps1 scripts/server/server_full.jl
```

或在 Windows 下使用：

```powershell
.\scripts\server\start.bat
```

默认监听 `http://localhost:8080`。

## 3) 打开前端

- 主界面：`web/index.html`
- 最小测试页：`web/simple_test.html`
- API 测试页：`web/test_api.html`

## 4) 基础校验

```powershell
julia --project=. tests/unit/runtests.jl
```

CI 默认使用 smoke 配置：

```powershell
$env:UNIT_PROFILE='smoke'
julia --project=. tests/unit/runtests.jl
```

## 5) 常见问题

- **端口占用**：`julia --project=. scripts/server/server_full.jl 8081`
- **前端离线**：访问 `http://localhost:8080/health` 检查服务是否存活
- **Three.js 加载失败**：查看 `web/THREEJS_LOCAL.md` 使用本地资源方案

## 6) 输出目录口径

- 默认结果目录：`data/outputs/`
- 根目录 `outputs/` 仅保留历史兼容，不作为默认落盘目录

## 7) 参考文档

- 安装与复现：`INSTALL.md`
- 使用说明：`docs/guides/USER_GUIDE.md`
- 状态说明：`docs/guides/STATUS.md`
- 脚本入口清单：`docs/guides/scripts/README.md`
- 守恒荷广义磁化率脚本：`docs/guides/scripts/pnjl_conserved_charge_susceptibilities.md`
- 参考资料：`docs/reference/`
