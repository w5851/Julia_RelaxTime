# 快速启动指南

本指南对应当前仓库结构（`scripts/server/`、`tests/unit/`、`docs/reference/`）。

## 1) 初始化环境

```powershell
julia --project=. -e "using Pkg; Pkg.instantiate(); Pkg.precompile()"
```

## 2) 启动服务（API + 前端静态资源）

```powershell
julia --project=. scripts/server/server_full.jl
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
- 参考资料：`docs/reference/`
