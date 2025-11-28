# 2→2散射Web可视化模块

本模块为Julia_RelaxTime项目新增了2→2散射过程的动量计算和椭球可视化功能。

## 快速开始

```powershell
# 1. 启动Julia HTTP服务器
julia server.jl

# 2. 打开浏览器访问
web/index.html
```

## 项目结构

```
├── server.jl                          # HTTP服务器启动脚本
├── src/simulation/                    # Julia计算模块
│   ├── FrameTransformations.jl       # Lorentz变换
│   ├── EllipsoidCalculation.jl       # 椭球参数计算
│   ├── MomentumMapping.jl            # 动量映射主模块
│   └── HTTPServer.jl                  # REST API服务器
├── web/                               # 前端文件
│   ├── index.html                     # 主页面
│   ├── css/style.css                  # 样式
│   └── js/
│       ├── api.js                     # API通信
│       ├── visualization.js           # Three.js可视化
│       └── ui.js                      # UI控制器
├── test_unit/                         # 单元测试
│   ├── test_frame_transformations.jl
│   └── test_momentum_mapping.jl
└── examples/
    └── web_demo.md                    # 完整使用文档
```

## 功能特性

✅ **物理计算**
- Lorentz坐标系变换 (CMS ↔ Lab)
- 动量可达域椭球计算
- 能量/动量守恒验证
- Mandelstam变量计算

✅ **Web可视化**
- Three.js 3D场景渲染
- 动量箭头显示
- 椭球表面网格
- 交互式相机控制

✅ **HTTP API**
- RESTful接口 (POST /compute)
- JSON序列化
- CORS跨域支持
- 健康检查端点

## 使用示例

### 命令行使用

```julia
using LinearAlgebra
include("src/simulation/MomentumMapping.jl")
using .MomentumMapping

# 设置入射动量
p1 = [0.5, 0.0, 1.8]  # fm⁻¹
p2 = [-0.5, 0.0, -1.8]
m = 1.52  # fm⁻¹

# 计算散射
result = calculate_outgoing_momenta(p1, p2, m, m, m, m, π/4, π/6)

# 验证物理约束
validate_kinematics(result, m, m, m, m)

# 打印结果
print_kinematics_info(result)
```

### Web界面使用

1. 启动服务器: `julia server.jl`
2. 打开 `web/index.html`
3. 输入参数（默认值已填充）
4. 点击"计算散射"
5. 查看右侧3D可视化和左侧物理量

## API文档

### POST /compute

**请求体**:
```json
{
  "p1x": 0.5, "p1y": 0.0, "p1z": 1.8,
  "p2x": -0.5, "p2y": 0.0, "p2z": -1.8,
  "m1": 1.52, "m2": 1.52, "m3": 1.52, "m4": 1.52,
  "theta_star": 0.785,
  "phi_star": 0.524
}
```

**响应**:
```json
{
  "success": true,
  "data": {
    "ellipsoid": {...},
    "momenta": {...},
    "physics": {...},
    "validation": {...}
  }
}
```

详见 `examples/web_demo.md` 完整文档。

## 测试

```powershell
# 运行单元测试
julia test_unit/test_frame_transformations.jl
julia test_unit/test_momentum_mapping.jl
```

## 技术栈

- **后端**: Julia + HTTP.jl + JSON3.jl
- **前端**: Three.js + 原生JavaScript
- **架构**: REST API + 前后端分离

## 物理原理

使用仿射变换计算出射动量在实验室系的可达域：

```
p_lab = A · p_cms + b
椭球 = { A · p_cms + b | |p_cms| = p* }
```

通过特征分解 `S = A·Aᵀ` 得到椭球的主轴方向和半轴长。

详细推导见 `doc/domain-knowledge/从质心系到实验室系动量转换.md`。

## 注意事项

- 确保 √s > m₃ + m₄（高于散射阈值）
- 端口8080需要可用（或修改为其他端口）
- 浏览器需要支持ES6模块和Three.js

## 相关文档

- 完整使用指南: `examples/web_demo.md`
- 实现计划: `.github/prompts/plan-twoToTwoScatteringMomentumVisualization.prompt.md`
- 物理推导: `doc/domain-knowledge/从质心系到实验室系动量转换.md`

---

**状态**: ✅ 全部12项任务已完成
