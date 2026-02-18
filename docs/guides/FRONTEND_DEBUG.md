# 前端问题修复指南

## 当前状态

✅ **Julia服务器**: 正常运行在 http://localhost:8080
✅ **API测试**: 健康检查和计算端点均正常
❓ **前端界面**: 需要验证

## 已修复的问题

### 1. HTML表单验证
**问题**: 输入框的 `step="0.01"` 和严格的 `min/max` 导致"请输入有效值"错误

**解决方案**: 
- 所有输入框改为 `step="any"` 允许任意精度
- theta_star max改为 3.15 (略大于π)
- phi_star max改为 6.29 (略大于2π)

### 2. 添加调试日志
在以下位置添加了console.log:
- `ui.js` - checkServerStatus: 服务器检查状态
- `ui.js` - collectFormData: 表单数据收集
- `api.js` - checkHealth: API健康检查
- `api.js` - computeScattering: 计算请求

### 3. 创建简化测试页面
`web/simple_test.html` - 用于快速验证API功能

## 测试步骤

### 方法1: 使用简化测试页面（推荐）
1. 确保Julia服务器在运行
2. 打开浏览器访问: `web/simple_test.html`
3. 页面应显示"✓ 服务器在线"
4. 点击"快速计算"按钮
5. 应看到计算结果

### 方法2: 使用完整界面
1. 打开浏览器访问: `web/index.html`
2. 按F12打开开发者工具 → Console标签
3. 查看状态指示器（应显示"服务器在线"）
4. 检查Console日志:
   ```
   [UI] Checking server health...
   [API] Checking health at http://localhost:8080/health
   [API] Health check response: 200 OK
   [UI] Server status: ONLINE
   ```
5. 点击"计算散射"按钮
6. 检查Console日志查看表单数据是否正确收集

## 当前测试结果

### 后端API测试 ✅
```powershell
PS> Invoke-WebRequest http://localhost:8080/health
Status: 200
Content: OK

PS> # 计算测试
Status: 200
√s = 4.8168 fm⁻¹
p* = 1.8682 fm⁻¹
β = 0.0
γ = 1.0
```

## 调试步骤

如果前端仍然显示"检查服务器"或"请输入有效值":

### 1. 检查服务器状态
```powershell
# 测试服务器
Invoke-WebRequest -Uri http://localhost:8080/health -UseBasicParsing
```

### 2. 检查浏览器控制台
按F12打开开发者工具，查看:
- **Console标签**: 查找红色错误信息
- **Network标签**: 查看HTTP请求是否发送
  - 点击"计算散射"
  - 应看到 `/compute` 请求
  - 查看请求和响应内容

### 3. 检查CORS
如果看到CORS错误:
```
Access-Control-Allow-Origin
```
确认Julia服务器的CORS中间件正常工作。

### 4. 检查输入值
在Console中运行:
```javascript
// 检查theta_star输入
document.getElementById('theta_star').value
document.getElementById('theta_star').validity.valid

// 检查所有输入
['p1x','p1y','p1z','p2x','p2y','p2z','m1','m2','m3','m4','theta_star','phi_star']
  .forEach(id => {
    const el = document.getElementById(id);
    console.log(id, '=', el.value, 'valid:', el.validity.valid);
  });
```

## 预期行为

### 正常流程
1. 页面加载 → 检查服务器 → 显示"服务器在线" 🟢
2. 输入参数（或使用默认值）
3. 点击"计算散射"
4. 显示"计算中..."
5. 3-5秒后显示结果
6. 右侧3D视图显示椭球和箭头

### Console日志示例
```
[UI] Checking server health...
[API] Checking health at http://localhost:8080/health
[API] Health check response: 200 OK
[UI] Server status: ONLINE
[UI] Collected p1x = 0.5 (raw: "0.5")
[UI] Collected p1y = 0 (raw: "0.0")
...
[UI] Collected form data: {p1x: 0.5, p1y: 0, ...}
[API] Sending compute request...
[API] Compute response: 200 OK
```

## 快速修复命令

如果需要重启所有组件:
```powershell
# 停止Julia
Get-Process julia -ErrorAction SilentlyContinue | Stop-Process -Force

# 重启服务器
Start-Process cmd -ArgumentList "/k", "cd /d D:\Desktop\Julia_RelaxTime && julia --project=. scripts/server/server_full.jl"

# 等待10秒
Start-Sleep -Seconds 10

# 测试
Invoke-WebRequest -Uri http://localhost:8080/health -UseBasicParsing

# 打开简化测试页面
start web\simple_test.html
```

## 文件清单

**测试文件**:
- ✅ `web/simple_test.html` - 简化测试页面（新建）
- ✅ `web/test_api.html` - API测试页面
- ✅ `web/index.html` - 完整界面（已修复）

**已修复的文件**:
- ✅ `web/index.html` - 所有输入改为step="any"
- ✅ `web/js/api.js` - 添加调试日志
- ✅ `web/js/ui.js` - 添加调试日志
- ✅ `src/simulation/FrameTransformations.jl` - 导出boost_energy
- ✅ `src/simulation/HTTPServer.jl` - 添加using LinearAlgebra
- ✅ `start.bat` - 修复编码和退出问题

---

**下一步**: 打开 `web/simple_test.html` 验证功能是否正常
