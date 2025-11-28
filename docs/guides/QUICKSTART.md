# 快速启动指南

## ✅ 当前状态

**Julia服务器**: 已成功启动并运行在 `http://localhost:8080`

## 🚀 使用系统

### 方法1: 直接打开前端（推荐）

1. **双击打开**: `web\index.html`
2. 在表单中输入参数（默认值已填充）
3. 点击 **"🚀 计算散射"** 按钮
4. 查看右侧3D可视化和左侧物理量

### 方法2: 使用一键启动脚本

```powershell
.\start.bat
```

这将自动：
- 检查并安装Julia依赖
- 启动HTTP服务器
- 打开前端页面

## 📦 依赖说明

### Julia依赖（已安装✅）
- `HTTP.jl` - HTTP服务器
- `JSON3.jl` - JSON序列化
- `LinearAlgebra` - 标准库

安装命令（已完成）:
```julia
using Pkg
Pkg.add(["HTTP", "JSON3"])
```

### JavaScript依赖（无需安装✅）
- **Three.js** - 从CDN加载（https://cdn.jsdelivr.net/npm/three@0.160.0/）
- **OrbitControls** - 从CDN加载

**优点**:
- ✅ 无需本地安装，打开即用
- ✅ 自动缓存，二次加载快
- ✅ 始终使用最新稳定版

**如需离线使用**: 参见 `web/THREEJS_LOCAL.md`

## 🧪 测试系统

运行单元测试：
```powershell
julia test_unit/test_frame_transformations.jl
julia test_unit/test_momentum_mapping.jl
```

## 📖 默认示例参数

- **p₁** = (0.5, 0.0, 1.8) fm⁻¹ - 入射粒子1
- **p₂** = (-0.5, 0.0, -1.8) fm⁻¹ - 入射粒子2
- **m** = 1.52 fm⁻¹ - 所有粒子质量
- **θ*** = π/4 = 0.785 rad - 质心系极角
- **φ*** = π/6 = 0.524 rad - 质心系方位角

预期结果：
- **√s** ≈ 4.81 fm⁻¹
- **p*** ≈ 1.23 fm⁻¹
- **β** ≈ 0 (质心系静止)
- **γ** ≈ 1.0

## 🎯 工作流程

```
输入动量参数
    ↓
[前端] 表单验证
    ↓
[API] POST /compute → Julia服务器
    ↓
[Julia] 计算质心系参数
    ↓
[Julia] Lorentz变换
    ↓
[Julia] 椭球参数计算
    ↓
[Julia] 物理验证
    ↓
[API] JSON响应 ← Julia服务器
    ↓
[Three.js] 渲染3D场景
    ↓
显示结果
```

## 🛠️ 故障排除

### 问题1: 服务器无法启动
**错误**: `Package HTTP not found`
**解决**: 
```powershell
julia --project=. -e "using Pkg; Pkg.add([\"HTTP\", \"JSON3\"])"
```

### 问题2: 端口被占用
**错误**: `Address already in use`
**解决**: 使用其他端口
```powershell
julia server.jl 8081
```
并在 `web/js/api.js` 中修改:
```javascript
const API_BASE_URL = 'http://localhost:8081';
```

### 问题3: 前端显示"服务器离线"
**检查**:
1. Julia服务器是否在运行？
2. 浏览器控制台是否有CORS错误？
3. 访问 http://localhost:8080/health 是否返回 OK？

### 问题4: 3D视图不显示
**原因**: Three.js未加载（网络问题）
**解决**: 
1. 检查网络连接
2. 查看浏览器控制台错误
3. 尝试其他浏览器（推荐Chrome/Edge）

## 📚 完整文档

- **使用指南**: `examples/web_demo.md`
- **API文档**: `examples/web_demo.md` → API文档章节
- **物理原理**: `doc/domain-knowledge/从质心系到实验室系动量转换.md`
- **实现计划**: `.github/prompts/plan-twoToTwoScatteringMomentumVisualization.prompt.md`

## 🎉 开始使用

现在可以直接打开 `web\index.html` 开始使用系统！

服务器已在后台运行，无需额外操作。
