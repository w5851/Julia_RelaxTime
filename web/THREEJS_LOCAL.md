# Three.js 本地化说明

## 当前状态

前端使用CDN加载Three.js库：
```javascript
import * as THREE from 'https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js';
import { OrbitControls } from 'https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/controls/OrbitControls.js';
```

## 如需本地化（可选）

### 方案1：下载Three.js到本地

```powershell
# 在 web/ 目录下创建 lib 文件夹
mkdir web/lib

# 下载Three.js文件（需要手动从 https://threejs.org 下载）
# 将以下文件放入 web/lib/:
# - three.module.js
# - OrbitControls.js
```

然后修改 `web/js/visualization.js`:
```javascript
import * as THREE from '../lib/three.module.js';
import { OrbitControls } from '../lib/OrbitControls.js';
```

### 方案2：使用当前CDN方式（推荐）

**优点**:
- ✅ 无需下载，文件更小
- ✅ 自动缓存，加载快
- ✅ 使用最新稳定版本

**缺点**:
- ❌ 需要网络连接

### 方案3：使用npm + 打包工具

如需完全离线使用，可使用Vite或Webpack打包：

```powershell
# 初始化npm项目
cd web
npm init -y
npm install three

# 使用Vite打包
npm install -D vite
npx vite build
```

## 当前推荐

**保持CDN方式**，因为：
1. Three.js库较大（~1MB），CDN缓存更高效
2. 项目是本地开发工具，通常有网络连接
3. 代码简洁，无需构建步骤

如果您在离线环境使用，请选择方案1。
