# 🎯 使用指南 - 2→2散射可视化系统

## 📋 当前状态

✅ **后端服务器**: 运行正常 (http://localhost:8080)
✅ **API端点**: 测试通过
✅ **前端修复**: 表单验证问题已解决

## 🚀 快速启动

### 方法1: 一键启动（推荐）
```powershell
.\start.bat
```

### 方法2: 手动启动
```powershell
# 步骤1: 启动Julia服务器
julia server.jl

# 步骤2: 在浏览器中打开
start web\index.html
```

## 🧪 测试页面

### 简化测试页面（推荐新手使用）
```powershell
start web\simple_test.html
```
**特点**:
- 简洁界面
- 一键测试服务器
- 实时显示计算结果
- 适合调试和验证

### 完整3D界面
```powershell
start web\index.html
```
**特点**:
- Three.js 3D可视化
- 完整物理量显示
- 动量箭头和椭球渲染
- 交互式参数调整

## ✅ 验证清单

### 1. 服务器检查
打开浏览器后，应该看到：
- [ ] 左上角显示 "服务器在线" 🟢（绿色指示器）
- [ ] "计算散射"按钮可点击（未禁用）

**如果显示"检查服务器..."一直不变**:
1. 按F12打开开发者工具
2. 查看Console标签的错误信息
3. 检查Network标签是否有请求失败

### 2. 输入验证
默认参数应该全部有效：
- [ ] p₁ = (0.5, 0.0, 1.8)
- [ ] p₂ = (-0.5, 0.0, -1.8)  
- [ ] m = 1.52（所有粒子）
- [ ] θ* = 0.785 rad (= π/4)
- [ ] φ* = 0.524 rad (= π/6)

**如果显示"请输入一个有效的值"**:
- 已修复：所有输入框现在使用 `step="any"`
- 刷新页面（Ctrl+F5）清除缓存
- 或手动修改数值后再改回来

### 3. 计算测试
点击"计算散射"后：
- [ ] 按钮显示"计算中..."
- [ ] 3-5秒后显示结果面板
- [ ] 右侧3D视图出现椭球和箭头（完整界面）

**预期结果**:
```
√s ≈ 4.817 fm⁻¹
p* ≈ 1.868 fm⁻¹
β ≈ 0.0
γ ≈ 1.0
```

## 🐛 故障排除

### 问题1: "一直显示'检查服务器...'"

**原因**: 服务器未启动或CORS问题

**解决**:
```powershell
# 1. 检查服务器进程
Get-Process julia

# 2. 测试API
Invoke-WebRequest -Uri http://localhost:8080/health -UseBasicParsing

# 3. 如果失败，重启服务器
Get-Process julia | Stop-Process -Force
julia server.jl
```

### 问题2: "请输入一个有效的值"

**原因**: HTML5表单验证过于严格

**解决**: ✅ 已修复
- 刷新页面（Ctrl+F5）
- 确认所有输入框显示数值
- 检查浏览器Console是否有JavaScript错误

### 问题3: 点击计算后无响应

**诊断步骤**:
1. 打开开发者工具（F12）
2. 切换到Console标签
3. 点击"计算散射"
4. 查看日志输出

**正常日志**:
```
[UI] Collected p1x = 0.5 (raw: "0.5")
[UI] Collected p1y = 0 (raw: "0.0")
...
[UI] Collected form data: {p1x: 0.5, ...}
```

**如果看到错误**:
- `NaN detected` → 某个输入框为空或无效
- `CORS error` → 服务器CORS配置问题
- `Network error` → 服务器连接失败

### 问题4: 3D视图不显示

**原因**: Three.js加载失败（需要网络）

**解决**:
1. 检查网络连接
2. 查看Console是否有Three.js加载错误
3. 使用简化测试页面（不需要Three.js）
4. 或参考 `web/THREEJS_LOCAL.md` 安装本地版本

## 📊 使用技巧

### 修改参数探索
1. **对称散射**: p₁ = -p₂ → β = 0，椭球在原点
2. **增大动量**: |p| → 更大 → √s 增大
3. **改变方向**: θ*, φ* → 改变出射动量方向

### 查看详细信息
按F12打开开发者工具：
- **Console**: 查看计算日志
- **Network**: 查看API请求/响应
- **Elements**: 检查DOM元素

### 导出结果
在Console中运行：
```javascript
// 获取最后的计算结果（需要在计算后运行）
console.log(JSON.stringify(lastResult, null, 2));
```

## 📱 浏览器兼容性

**推荐浏览器**:
- ✅ Chrome 90+
- ✅ Edge 90+
- ✅ Firefox 88+
- ⚠️ Safari 14+（可能需要调整）

**要求**:
- ES6模块支持
- Fetch API
- WebGL（用于3D视图）

## 🔧 高级调试

### 启用详细日志
1. 打开 `web/js/api.js`
2. 所有console.log已添加
3. 查看浏览器Console输出

### 手动测试API
在浏览器Console中运行：
```javascript
// 测试健康检查
fetch('http://localhost:8080/health')
  .then(r => r.text())
  .then(console.log);

// 测试计算
fetch('http://localhost:8080/compute', {
  method: 'POST',
  headers: {'Content-Type': 'application/json'},
  body: JSON.stringify({
    p1x:0.5, p1y:0, p1z:1.8,
    p2x:-0.5, p2y:0, p2z:-1.8,
    m1:1.52, m2:1.52, m3:1.52, m4:1.52,
    theta_star:0.785, phi_star:0.524
  })
})
.then(r => r.json())
.then(console.log);
```

## 📞 获取帮助

如果问题仍然存在：

1. **查看日志文件**:
   - `FIXES.md` - 已修复问题列表
   - `FRONTEND_DEBUG.md` - 详细调试指南
   - `QUICKSTART.md` - 快速开始指南

2. **收集信息**:
   - Julia服务器输出
   - 浏览器Console日志
   - Network标签的请求详情

3. **简化测试**:
   - 使用 `web/simple_test.html`
   - 使用 `web/test_api.html`
   - 直接用PowerShell测试API

---

## ✨ 总结

**修复的问题**:
✅ bat脚本编码问题
✅ HTML表单验证过严
✅ Julia模块导出缺失
✅ JavaScript兼容性问题

**当前状态**:
✅ 服务器正常运行
✅ API测试通过
✅ 所有已知问题已修复

**下一步**:
1. 打开 `web/simple_test.html` 快速验证
2. 或打开 `web/index.html` 使用完整界面
3. 按F12查看Console确认一切正常

**祝使用愉快！** 🎉
