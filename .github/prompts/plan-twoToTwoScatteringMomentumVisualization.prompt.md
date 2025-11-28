# Plan: 2→2散射动量计算与椭球可视化

根据入射粒子三动量和质心系方向参数,通过洛伦兹变换计算出射粒子三动量,并可视化实验室系中的动量可行域椭球。新建`src/simulation/`文件夹存放代码。

## Steps

### Phase 1: 后端核心计算 (Julia)

1. **创建 `src/simulation/FrameTransformations.jl`** - 实现洛伦兹变换核心功能:计算质心系速度 $\mathbf{v}_{\rm CM} = \frac{\mathbf{p}_1 + \mathbf{p}_2}{E_1 + E_2}$、洛伦兹因子 $\gamma$;实现boost公式 $\mathbf{p}_{\rm lab} = A\mathbf{p}^* + \mathbf{b}$ (基于`doc/domain-knowledge/从质心系到实验室系动量转换.md`);计算Källen函数和Mandelstam变量s

2. **创建 `src/simulation/EllipsoidCalculation.jl`** - 计算椭球参数:根据仿射变换矩阵 $A = I + \frac{\gamma-1}{\beta^2}\boldsymbol{\beta}\boldsymbol{\beta}^T$ 和偏移 $\mathbf{b} = \gamma E^* \boldsymbol{\beta}$ 计算椭球几何;通过特征分解 $S = AA^T = U\,\mathrm{diag}(\lambda_i)U^T$ 得到主轴方向和半轴长;支持在椭球表面均匀采样点(通过质心系角度 $\theta^*, \phi^*$)

3. **创建 `src/simulation/MomentumMapping.jl`** - 主控模块实现动量映射:输入实验室系入射三动量 $(\mathbf{p}_1, \mathbf{p}_2)$ 和粒子质量;根据质心系方向 $(\theta^*, \phi^*)$ 计算 $\mathbf{p}_3^* = p^*(\sin\theta^*\cos\phi^*, \sin\theta^*\sin\phi^*, \cos\theta^*)$;应用boost变换得到 $\mathbf{p}_3$ 和 $\mathbf{p}_4 = \mathbf{p}_1 + \mathbf{p}_2 - \mathbf{p}_3$;验证四动量守恒和质壳条件

4. **创建 `src/simulation/HTTPServer.jl`** - 实现轻量级HTTP服务器(使用`HTTP.jl`和`JSON3.jl`):提供REST API接口 `POST /compute`;接收JSON格式的入射动量和质量参数;调用`MomentumMapping`计算椭球参数;返回JSON响应(椭球中心、主轴方向、半轴长、物理量);支持CORS跨域请求;错误处理(阈值以下、无效输入)

5. **创建测试文件 `test_unit/test_simulation_*.jl`** - 遵循现有测试模式:验证洛伦兹变换的boost-inverse一致性;检查椭球参数计算(半轴长、主轴方向正交性);验证四动量守恒 $|E_1+E_2-E_3-E_4| < 10^{-10}$ 和质壳条件 $|E_i^2-\mathbf{p}_i^2-m_i^2| < 10^{-10}$;测试阈值边界情况 $s \to (m_3+m_4)^2$;性能基准(单次计算 <100μs)

### Phase 2: 前端可视化 (JavaScript + three.js)

6. **创建 `web/index.html`** - 主页面结构:输入表单(入射动量p1, p2,粒子质量);控制面板(显示/隐藏元素,调整透明度);3D画布容器;信息显示区域(物理量)

7. **创建 `web/js/api.js`** - 后端通信模块:封装fetch API调用Julia服务器;发送POST请求到`http://localhost:8080/compute`;处理JSON响应和错误;提供加载状态反馈

8. **创建 `web/js/visualization.js`** - three.js核心渲染:初始化场景、相机、渲染器;绘制椭球(通过主轴变换SphereGeometry);绘制动量矢量箭头(ArrowHelper);添加坐标轴和网格;实现鼠标控制(OrbitControls)

9. **创建 `web/js/ui.js`** - 用户交互逻辑:处理表单提交事件;解析输入并验证;调用api.js获取数据;更新visualization.js渲染;显示物理量信息

10. **创建 `web/css/style.css`** - 页面样式:布局(左侧控制面板+右侧3D视图);输入框和按钮样式;信息面板样式;响应式设计

### Phase 3: 集成与测试

11. **创建启动脚本 `server.jl`** - 一键启动服务器:加载所有模块;启动HTTP服务器(端口8080);打印访问地址;优雅关闭处理

12. **创建示例与文档** - `examples/web_demo.md`:说明如何启动服务器和打开浏览器;提供示例输入参数;截图展示效果;常见问题解答

## Further Considerations

1. **数值精度处理** - 在 $\beta^2 \to 0$ (近静止系)和 $s \to (m_3+m_4)^2$ (阈值附近)时的数值稳定性;是否需要特殊处理 $\gamma-1$ 项避免cancelation error?

2. **前端性能优化** - 椭球采样点数量(前端通过参数化曲面生成 vs 后端传输大量点);WebGL着色器优化;大量动量矢量的批量渲染(InstancedMesh)

3. **交互功能扩展** (按优先级排序):
   - **P0 (必需)**: 输入入射动量 → 计算 → 显示椭球和矢量
   - **P1 (重要)**: 在椭球上点击/拖动选择出射方向 → 更新p3, p4矢量
   - **P2 (有用)**: 实时调整入射动量(滑块) → 动态更新椭球
   - **P3 (可选)**: 动画演示洛伦兹boost过程(CMS → Lab渐变)
   - **P4 (未来)**: 并排显示质心系和实验室系视图

4. **部署与扩展** - 未来是否需要部署到云服务器(多用户支持)?是否需要保存/加载配置?是否需要导出图片/视频?

## 需求分析

### 核心功能需求

1. **入射态定义**
   - 输入:实验室系中两个入射粒子的三动量 $\mathbf{p}_1, \mathbf{p}_2$ (3-vector)
   - 粒子质量:m₁, m₂, m₃, m₄ (可从`Constants_PNJL.jl`获取或手动指定)
   - 能量计算: $E_i = \sqrt{\mathbf{p}_i^2 + m_i^2}$

2. **质心系参数**
   - 质心系速度(矢量): $\boldsymbol{\beta} = \mathbf{v}_{\rm CM} = \frac{\mathbf{p}_1 + \mathbf{p}_2}{E_1 + E_2}$
   - 洛伦兹因子: $\gamma = \frac{1}{\sqrt{1-\beta^2}}$, 其中 $\beta^2 = \boldsymbol{\beta} \cdot \boldsymbol{\beta}$
   - Mandelstam s: $s = (E_1 + E_2)^2 - |\mathbf{p}_1 + \mathbf{p}_2|^2$
   - 质心系动量模: $p^* = \frac{\sqrt{\lambda(s, m_3^2, m_4^2)}}{2\sqrt{s}}$
   - 质心系能量: $E_3^* = \sqrt{(p^*)^2 + m_3^2}$, $E_4^* = \sqrt{s} - E_3^*$

3. **椭球参数计算**
   - 仿射变换矩阵: $A = I + \frac{\gamma-1}{\beta^2}\boldsymbol{\beta}\boldsymbol{\beta}^T$ (当 $\beta^2 < 10^{-14}$ 时 $A = I$)
   - 偏移向量: $\mathbf{b} = \gamma E_3^* \boldsymbol{\beta}$
   - 形状矩阵: $S = AA^T$
   - 特征分解: $S = U\Lambda U^T$, 主轴方向为U的列向量,半轴长为 $p^* \sqrt{\lambda_i}$

4. **动量映射**
   - 质心系方向参数化: $\mathbf{p}_3^* = p^*(\sin\theta^*\cos\phi^*, \sin\theta^*\sin\phi^*, \cos\theta^*)$
   - Boost到实验室系: $\mathbf{p}_3 = A\mathbf{p}_3^* + \mathbf{b}$
   - 能量计算: $E_3 = \gamma(E_3^* + \mathbf{p}_3^* \cdot \boldsymbol{\beta})$
   - 四动量守恒: $\mathbf{p}_4 = \mathbf{p}_1 + \mathbf{p}_2 - \mathbf{p}_3$, $E_4 = E_1 + E_2 - E_3$

5. **可视化需求**
   - **椭球绘制**: 半透明3D椭球表面,显示出射粒子动量可行域
   - **动量矢量**: 绘制入射 $\mathbf{p}_1, \mathbf{p}_2$ (红色箭头)和出射 $\mathbf{p}_3, \mathbf{p}_4$ (蓝色箭头)
   - **采样点**: 在椭球表面随机采样点验证可行域
   - **投影视图**: xy, xz, yz平面投影(椭圆轮廓)
   - **物理量标注**: 显示s, $p^*$, $\gamma$, $\beta$等关键参数

### 技术架构设计

#### 模块职责划分

**`FrameTransformations.jl`** (洛伦兹变换工具)
```julia
using LinearAlgebra

# Källen函数
kallen_lambda(a, b, c) = a^2 + b^2 + c^2 - 2a*b - 2a*c - 2b*c

# 计算质心系参数
function calculate_cms_parameters(p1::Vector{Float64}, p2::Vector{Float64}, 
                                  m1::Float64, m2::Float64, 
                                  m3::Float64, m4::Float64)
    E1 = sqrt(dot(p1, p1) + m1^2)
    E2 = sqrt(dot(p2, p2) + m2^2)
    Ptot = p1 + p2
    Etot = E1 + E2
    
    s = Etot^2 - dot(Ptot, Ptot)
    @assert s > (m3 + m4)^2 "Below threshold"
    
    β = Ptot / Etot
    β2 = dot(β, β)
    γ = 1 / sqrt(1 - β2)
    
    sqrt_s = sqrt(s)
    p_star = sqrt(kallen_lambda(s, m3^2, m4^2)) / (2sqrt_s)
    E3_star = sqrt(p_star^2 + m3^2)
    E4_star = sqrt_s - E3_star
    
    return (s=s, p_star=p_star, E3_star=E3_star, E4_star=E4_star,
            beta=β, gamma=γ, beta2=β2)
end

# 构造仿射变换
function build_affine_transform(β::Vector{Float64}, γ::Float64, E_star::Float64)
    β2 = dot(β, β)
    if β2 < 1e-14
        A = Matrix{Float64}(I, 3, 3)
    else
        A = Matrix{Float64}(I, 3, 3) + ((γ - 1) / β2) * (β * β')
    end
    b = γ * E_star * β
    return A, b
end

# 质心系到实验室系boost
function boost_to_lab(p_star::Vector{Float64}, A::Matrix{Float64}, b::Vector{Float64})
    return A * p_star + b
end
```

**`EllipsoidCalculation.jl`** (椭球几何计算)
```julia
using LinearAlgebra

# 椭球参数结构
struct EllipsoidParams
    center::Vector{Float64}      # 中心坐标 b
    axes_directions::Matrix{Float64}  # 主轴方向 (列向量)
    half_lengths::Vector{Float64}     # 半轴长
end

# 计算椭球参数
function calculate_ellipsoid_parameters(A::Matrix{Float64}, p_star::Float64)
    S = A * A'
    eigen_result = eigen(Symmetric(S))
    evals = eigen_result.values
    evecs = eigen_result.vectors
    
    half_lengths = p_star * sqrt.(evals)
    
    return EllipsoidParams(zeros(3), evecs, half_lengths)
end

# 在椭球表面采样点(均匀角度采样)
function sample_ellipsoid_surface(N::Int, p_star::Float64, 
                                  A::Matrix{Float64}, b::Vector{Float64})
    points = Vector{Vector{Float64}}(undef, N)
    for i in 1:N
        cos_theta = 2rand() - 1  # uniform in [-1, 1]
        phi = 2π * rand()
        theta = acos(cos_theta)
        
        # 质心系单位向量
        p_cms = p_star * [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]
        # boost到实验室系
        points[i] = A * p_cms + b
    end
    return points
end
```

**`MomentumMapping.jl`** (主控模块)
```julia
include("FrameTransformations.jl")
include("EllipsoidCalculation.jl")

# 完整计算结果
struct ScatteringKinematics
    # 入射态
    p1_lab::Vector{Float64}
    p2_lab::Vector{Float64}
    E1::Float64
    E2::Float64
    
    # 质心系参数
    s::Float64
    p_star::Float64
    beta::Vector{Float64}
    gamma::Float64
    
    # 出射态 (实验室系)
    p3_lab::Vector{Float64}
    p4_lab::Vector{Float64}
    E3::Float64
    E4::Float64
    
    # 质心系角度
    theta_star::Float64
    phi_star::Float64
    
    # 椭球参数
    ellipsoid::EllipsoidParams
    affine_A::Matrix{Float64}
    affine_b::Vector{Float64}
end

# 主函数: 给定入射动量和质心系角度,计算出射动量
function calculate_outgoing_momenta(
    p1::Vector{Float64}, p2::Vector{Float64},
    m1::Float64, m2::Float64, m3::Float64, m4::Float64,
    theta_star::Float64, phi_star::Float64
)
    # 1. 计算质心系参数
    cms = calculate_cms_parameters(p1, p2, m1, m2, m3, m4)
    
    # 2. 构造仿射变换
    A, b = build_affine_transform(cms.beta, cms.gamma, cms.E3_star)
    
    # 3. 计算质心系出射动量
    p3_cms = cms.p_star * [sin(theta_star)*cos(phi_star), 
                            sin(theta_star)*sin(phi_star), 
                            cos(theta_star)]
    
    # 4. Boost到实验室系
    p3_lab = boost_to_lab(p3_cms, A, b)
    E3 = cms.gamma * (cms.E3_star + dot(p3_cms, cms.beta))
    
    # 5. 四动量守恒
    E1 = sqrt(dot(p1, p1) + m1^2)
    E2 = sqrt(dot(p2, p2) + m2^2)
    p4_lab = p1 + p2 - p3_lab
    E4 = E1 + E2 - E3
    
    # 6. 计算椭球参数
    ellipsoid = calculate_ellipsoid_parameters(A, cms.p_star)
    ellipsoid = EllipsoidParams(b, ellipsoid.axes_directions, ellipsoid.half_lengths)
    
    return ScatteringKinematics(
        p1, p2, E1, E2,
        cms.s, cms.p_star, cms.beta, cms.gamma,
        p3_lab, p4_lab, E3, E4,
        theta_star, phi_star,
        ellipsoid, A, b
    )
end

# 验证物理约束
function validate_kinematics(result::ScatteringKinematics, 
                             m1::Float64, m2::Float64, 
                             m3::Float64, m4::Float64)
    # 四动量守恒
    ΔE = abs(result.E1 + result.E2 - result.E3 - result.E4)
    Δp = norm(result.p1_lab + result.p2_lab - result.p3_lab - result.p4_lab)
    
    # 质壳条件
    shell3 = abs(result.E3^2 - dot(result.p3_lab, result.p3_lab) - m3^2)
    shell4 = abs(result.E4^2 - dot(result.p4_lab, result.p4_lab) - m4^2)
    
    @assert ΔE < 1e-10 "Energy conservation violated: ΔE = $ΔE"
    @assert Δp < 1e-10 "Momentum conservation violated: Δp = $Δp"
    @assert shell3 < 1e-10 "On-shell condition violated for p3: $shell3"
    @assert shell4 < 1e-10 "On-shell condition violated for p4: $shell4"
    
    return true
end
```

**`HTTPServer.jl`** (HTTP服务器模块)
```julia
using HTTP
using JSON3

include("MomentumMapping.jl")

# API响应结构
struct APIResponse
    success::Bool
    data::Union{Dict, Nothing}
    error::Union{String, Nothing}
end

# 处理计算请求
function handle_compute(req::HTTP.Request)
    try
        # 解析请求体
        body = JSON3.read(String(req.body))
        p1 = [body.p1x, body.p1y, body.p1z]
        p2 = [body.p2x, body.p2y, body.p2z]
        m1, m2, m3, m4 = body.m1, body.m2, body.m3, body.m4
        theta_star = get(body, :theta_star, π/4)
        phi_star = get(body, :phi_star, 0.0)
        
        # 计算椭球参数
        result = calculate_outgoing_momenta(p1, p2, m1, m2, m3, m4, 
                                           theta_star, phi_star)
        
        # 构造响应数据
        response_data = Dict(
            "ellipsoid" => Dict(
                "center" => result.ellipsoid.center,
                "axes_directions" => [result.ellipsoid.axes_directions[:, i] 
                                      for i in 1:3],
                "half_lengths" => result.ellipsoid.half_lengths
            ),
            "momenta" => Dict(
                "p1" => result.p1_lab,
                "p2" => result.p2_lab,
                "p3" => result.p3_lab,
                "p4" => result.p4_lab
            ),
            "physics" => Dict(
                "s" => result.s,
                "sqrt_s" => sqrt(result.s),
                "p_star" => result.p_star,
                "beta" => norm(result.beta),
                "gamma" => result.gamma
            )
        )
        
        resp = APIResponse(true, response_data, nothing)
        return HTTP.Response(200, JSON3.write(resp))
        
    catch e
        error_msg = sprint(showerror, e)
        resp = APIResponse(false, nothing, error_msg)
        return HTTP.Response(400, JSON3.write(resp))
    end
end

# 启动服务器
function start_server(port::Int=8080)
    # CORS中间件
    function cors_middleware(handler)
        return function(req::HTTP.Request)
            if req.method == "OPTIONS"
                return HTTP.Response(200, [
                    "Access-Control-Allow-Origin" => "*",
                    "Access-Control-Allow-Methods" => "POST, OPTIONS",
                    "Access-Control-Allow-Headers" => "Content-Type"
                ])
            end
            
            resp = handler(req)
            HTTP.setheader(resp, "Access-Control-Allow-Origin" => "*")
            return resp
        end
    end
    
    # 路由
    router = HTTP.Router()
    HTTP.register!(router, "POST", "/compute", handle_compute)
    
    # 启动
    server = HTTP.serve!(cors_middleware(router), "0.0.0.0", port)
    println("\n" * "="^60)
    println("服务器已启动: http://localhost:$port")
    println("API端点: POST http://localhost:$port/compute")
    println("按 Ctrl+C 停止服务器")
    println("="^60 * "\n")
    
    return server
end
```

### 前端架构设计 (JavaScript + three.js)

**`web/js/api.js`** (后端API通信)
```javascript
const API_BASE_URL = 'http://localhost:8080';

class ScatteringAPI {
    async compute(params) {
        const response = await fetch(`${API_BASE_URL}/compute`, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify(params)
        });
        
        if (!response.ok) {
            throw new Error(`HTTP error! status: ${response.status}`);
        }
        
        const result = await response.json();
        if (!result.success) {
            throw new Error(result.error || 'Computation failed');
        }
        
        return result.data;
    }
}

export default new ScatteringAPI();
```

**`web/js/visualization.js`** (three.js可视化)
```javascript
import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';

class ScatteringVisualizer {
    constructor(containerId) {
        this.container = document.getElementById(containerId);
        this.scene = new THREE.Scene();
        this.camera = new THREE.PerspectiveCamera(
            75, 
            this.container.clientWidth / this.container.clientHeight,
            0.1, 
            1000
        );
        
        this.renderer = new THREE.WebGLRenderer({ antialias: true });
        this.renderer.setSize(this.container.clientWidth, this.container.clientHeight);
        this.container.appendChild(this.renderer.domElement);
        
        // 控制器
        this.controls = new OrbitControls(this.camera, this.renderer.domElement);
        this.camera.position.set(5, 5, 5);
        this.controls.update();
        
        // 坐标轴和网格
        this.scene.add(new THREE.AxesHelper(5));
        this.scene.add(new THREE.GridHelper(10, 10));
        
        // 光源
        const ambientLight = new THREE.AmbientLight(0xffffff, 0.6);
        const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
        directionalLight.position.set(10, 10, 10);
        this.scene.add(ambientLight);
        this.scene.add(directionalLight);
        
        // 存储对象引用以便更新
        this.ellipsoidMesh = null;
        this.arrowHelpers = [];
        
        this.animate();
    }
    
    // 绘制椭球
    drawEllipsoid(center, axesDirections, halfLengths) {
        // 清除旧椭球
        if (this.ellipsoidMesh) {
            this.scene.remove(this.ellipsoidMesh);
        }
        
        // 创建单位球
        const geometry = new THREE.SphereGeometry(1, 32, 32);
        const material = new THREE.MeshPhongMaterial({
            color: 0x00aaff,
            transparent: true,
            opacity: 0.3,
            side: THREE.DoubleSide
        });
        
        this.ellipsoidMesh = new THREE.Mesh(geometry, material);
        
        // 缩放至半轴长
        this.ellipsoidMesh.scale.set(
            halfLengths[0],
            halfLengths[1],
            halfLengths[2]
        );
        
        // 旋转矩阵 (axesDirections的列向量为主轴)
        const rotationMatrix = new THREE.Matrix4();
        rotationMatrix.set(
            axesDirections[0][0], axesDirections[1][0], axesDirections[2][0], 0,
            axesDirections[0][1], axesDirections[1][1], axesDirections[2][1], 0,
            axesDirections[0][2], axesDirections[1][2], axesDirections[2][2], 0,
            0, 0, 0, 1
        );
        this.ellipsoidMesh.applyMatrix4(rotationMatrix);
        
        // 平移至中心
        this.ellipsoidMesh.position.set(center[0], center[1], center[2]);
        
        this.scene.add(this.ellipsoidMesh);
    }
    
    // 绘制动量矢量箭头
    drawMomentumArrows(momenta) {
        // 清除旧箭头
        this.arrowHelpers.forEach(arrow => this.scene.remove(arrow));
        this.arrowHelpers = [];
        
        const origin = new THREE.Vector3(0, 0, 0);
        const arrows = [
            { momentum: momenta.p1, color: 0xff0000, label: 'p1' },
            { momentum: momenta.p2, color: 0xff0000, label: 'p2' },
            { momentum: momenta.p3, color: 0x0000ff, label: 'p3' },
            { momentum: momenta.p4, color: 0x0000ff, label: 'p4' }
        ];
        
        arrows.forEach(({ momentum, color }) => {
            const dir = new THREE.Vector3(...momentum).normalize();
            const length = Math.sqrt(
                momentum[0]**2 + momentum[1]**2 + momentum[2]**2
            );
            const arrowHelper = new THREE.ArrowHelper(
                dir, origin, length, color, 0.2, 0.1
            );
            this.scene.add(arrowHelper);
            this.arrowHelpers.push(arrowHelper);
        });
    }
    
    // 更新可视化
    update(data) {
        this.drawEllipsoid(
            data.ellipsoid.center,
            data.ellipsoid.axes_directions,
            data.ellipsoid.half_lengths
        );
        this.drawMomentumArrows(data.momenta);
    }
    
    // 动画循环
    animate() {
        requestAnimationFrame(() => this.animate());
        this.controls.update();
        this.renderer.render(this.scene, this.camera);
    }
}

export default ScatteringVisualizer;
```

**`web/js/ui.js`** (用户界面逻辑)
```javascript
import ScatteringAPI from './api.js';
import ScatteringVisualizer from './visualization.js';

class UIController {
    constructor() {
        this.visualizer = new ScatteringVisualizer('canvas-container');
        this.setupEventListeners();
    }
    
    setupEventListeners() {
        const form = document.getElementById('input-form');
        form.addEventListener('submit', async (e) => {
            e.preventDefault();
            await this.handleCompute();
        });
    }
    
    async handleCompute() {
        const statusDiv = document.getElementById('status');
        const physicsDiv = document.getElementById('physics-info');
        
        try {
            statusDiv.textContent = '计算中...';
            statusDiv.className = 'status loading';
            
            // 读取表单输入
            const params = {
                p1x: parseFloat(document.getElementById('p1x').value),
                p1y: parseFloat(document.getElementById('p1y').value),
                p1z: parseFloat(document.getElementById('p1z').value),
                p2x: parseFloat(document.getElementById('p2x').value),
                p2y: parseFloat(document.getElementById('p2y').value),
                p2z: parseFloat(document.getElementById('p2z').value),
                m1: parseFloat(document.getElementById('m1').value),
                m2: parseFloat(document.getElementById('m2').value),
                m3: parseFloat(document.getElementById('m3').value),
                m4: parseFloat(document.getElementById('m4').value)
            };
            
            // 调用API
            const data = await ScatteringAPI.compute(params);
            
            // 更新可视化
            this.visualizer.update(data);
            
            // 显示物理量信息
            physicsDiv.innerHTML = `
                <h3>物理量</h3>
                <p>s = ${data.physics.s.toFixed(4)} fm⁻²</p>
                <p>√s = ${data.physics.sqrt_s.toFixed(4)} fm⁻¹ 
                   (${(data.physics.sqrt_s * 197.327).toFixed(2)} MeV)</p>
                <p>p* = ${data.physics.p_star.toFixed(4)} fm⁻¹</p>
                <p>β = ${data.physics.beta.toFixed(4)}</p>
                <p>γ = ${data.physics.gamma.toFixed(4)}</p>
            `;
            
            statusDiv.textContent = '计算成功！';
            statusDiv.className = 'status success';
            
        } catch (error) {
            statusDiv.textContent = `错误: ${error.message}`;
            statusDiv.className = 'status error';
            console.error('Computation error:', error);
        }
    }
}

// 页面加载完成后初始化
window.addEventListener('DOMContentLoaded', () => {
    new UIController();
});
```

**`web/index.html`** (主页面)
```html
<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>2→2散射动量可视化</title>
    <link rel="stylesheet" href="css/style.css">
</head>
<body>
    <div class="container">
        <aside class="control-panel">
            <h1>散射参数输入</h1>
            <form id="input-form">
                <fieldset>
                    <legend>入射动量 p1 (fm⁻¹)</legend>
                    <input type="number" id="p1x" step="0.01" value="0.5" required>
                    <input type="number" id="p1y" step="0.01" value="0.0" required>
                    <input type="number" id="p1z" step="0.01" value="1.8" required>
                </fieldset>
                
                <fieldset>
                    <legend>入射动量 p2 (fm⁻¹)</legend>
                    <input type="number" id="p2x" step="0.01" value="-0.5" required>
                    <input type="number" id="p2y" step="0.01" value="0.0" required>
                    <input type="number" id="p2z" step="0.01" value="-1.8" required>
                </fieldset>
                
                <fieldset>
                    <legend>粒子质量 (fm⁻¹)</legend>
                    <label>m1: <input type="number" id="m1" step="0.01" value="1.52" required></label>
                    <label>m2: <input type="number" id="m2" step="0.01" value="1.52" required></label>
                    <label>m3: <input type="number" id="m3" step="0.01" value="1.52" required></label>
                    <label>m4: <input type="number" id="m4" step="0.01" value="1.52" required></label>
                </fieldset>
                
                <button type="submit">计算并可视化</button>
            </form>
            
            <div id="status" class="status"></div>
            <div id="physics-info" class="info-panel"></div>
        </aside>
        
        <main id="canvas-container" class="canvas-container"></main>
    </div>
    
    <script type="module" src="js/ui.js"></script>
</body>
</html>
```

**`web/css/style.css`** (样式)
```css
* {
    margin: 0;
    padding: 0;
    box-sizing: border-box;
}

body {
    font-family: 'Arial', sans-serif;
    overflow: hidden;
}

.container {
    display: flex;
    height: 100vh;
}

.control-panel {
    width: 350px;
    padding: 20px;
    background: #f5f5f5;
    overflow-y: auto;
    border-right: 2px solid #ddd;
}

.control-panel h1 {
    font-size: 20px;
    margin-bottom: 20px;
    color: #333;
}

fieldset {
    margin-bottom: 15px;
    padding: 10px;
    border: 1px solid #ccc;
    border-radius: 5px;
}

legend {
    font-weight: bold;
    color: #555;
}

input[type="number"] {
    width: 100%;
    padding: 8px;
    margin: 5px 0;
    border: 1px solid #ddd;
    border-radius: 3px;
}

button {
    width: 100%;
    padding: 12px;
    background: #0066cc;
    color: white;
    border: none;
    border-radius: 5px;
    font-size: 16px;
    cursor: pointer;
    margin-top: 10px;
}

button:hover {
    background: #0052a3;
}

.status {
    margin-top: 15px;
    padding: 10px;
    border-radius: 5px;
    font-weight: bold;
}

.status.loading {
    background: #fff3cd;
    color: #856404;
}

.status.success {
    background: #d4edda;
    color: #155724;
}

.status.error {
    background: #f8d7da;
    color: #721c24;
}

.info-panel {
    margin-top: 15px;
    padding: 10px;
    background: white;
    border-radius: 5px;
    border: 1px solid #ddd;
}

.info-panel h3 {
    margin-bottom: 10px;
    color: #333;
}

.info-panel p {
    margin: 5px 0;
    font-size: 14px;
}

.canvas-container {
    flex: 1;
    position: relative;
}
```

**`ScatteringVisualization.jl`** (Julia后端可视化 - 可选保留用于测试)
```julia
using Plots

# 绘制3D椭球和动量矢量
function plot_scattering_3d(result::ScatteringKinematics; n_samples=200)
    # 采样椭球表面
    points = sample_ellipsoid_surface(n_samples, result.p_star, 
                                      result.affine_A, result.affine_b)
    px = [p[1] for p in points]
    py = [p[2] for p in points]
    pz = [p[3] for p in points]
    
    # 创建3D散点图(椭球)
    p = scatter(px, py, pz, 
                markersize=1, alpha=0.3, color=:lightblue,
                label="可行域椭球", legend=:topright)
    
    # 绘制入射动量(红色箭头)
    quiver!([0], [0], [0], 
            quiver=([result.p1_lab[1]], [result.p1_lab[2]], [result.p1_lab[3]]),
            color=:red, linewidth=2, label="p₁入射")
    quiver!([0], [0], [0],
            quiver=([result.p2_lab[1]], [result.p2_lab[2]], [result.p2_lab[3]]),
            color=:red, linewidth=2, label="p₂入射")
    
    # 绘制出射动量(蓝色箭头)
    quiver!([0], [0], [0],
            quiver=([result.p3_lab[1]], [result.p3_lab[2]], [result.p3_lab[3]]),
            color=:blue, linewidth=2, label="p₃出射")
    quiver!([0], [0], [0],
            quiver=([result.p4_lab[1]], [result.p4_lab[2]], [result.p4_lab[3]]),
            color=:blue, linewidth=2, label="p₄出射")
    
    xlabel!("px [fm⁻¹]")
    ylabel!("py [fm⁻¹]")
    zlabel!("pz [fm⁻¹]")
    title!("2→2散射动量空间 (s=$(round(result.s, digits=3)) fm⁻²)")
    
    return p
end

# 绘制投影视图
function plot_projections(result::ScatteringKinematics; n_samples=200)
    points = sample_ellipsoid_surface(n_samples, result.p_star,
                                      result.affine_A, result.affine_b)
    
    # xy投影
    p1 = scatter([p[1] for p in points], [p[2] for p in points],
                 markersize=1, alpha=0.3, label="椭球投影")
    quiver!([0, 0], [0, 0], 
            quiver=([result.p1_lab[1], result.p2_lab[1]], 
                    [result.p1_lab[2], result.p2_lab[2]]),
            color=:red, linewidth=2, label="入射")
    quiver!([0, 0], [0, 0],
            quiver=([result.p3_lab[1], result.p4_lab[1]],
                    [result.p3_lab[2], result.p4_lab[2]]),
            color=:blue, linewidth=2, label="出射")
    xlabel!("px"); ylabel!("py"); title!("xy投影")
    
    # xz投影
    p2 = scatter([p[1] for p in points], [p[3] for p in points],
                 markersize=1, alpha=0.3, label="")
    quiver!([0, 0], [0, 0],
            quiver=([result.p1_lab[1], result.p2_lab[1]],
                    [result.p1_lab[3], result.p2_lab[3]]),
            color=:red, linewidth=2, label="")
    quiver!([0, 0], [0, 0],
            quiver=([result.p3_lab[1], result.p4_lab[1]],
                    [result.p3_lab[3], result.p4_lab[3]]),
            color=:blue, linewidth=2, label="")
    xlabel!("px"); ylabel!("pz"); title!("xz投影")
    
    # yz投影
    p3 = scatter([p[2] for p in points], [p[3] for p in points],
                 markersize=1, alpha=0.3, label="")
    quiver!([0, 0], [0, 0],
            quiver=([result.p1_lab[2], result.p2_lab[2]],
                    [result.p1_lab[3], result.p2_lab[3]]),
            color=:red, linewidth=2, label="")
    quiver!([0, 0], [0, 0],
            quiver=([result.p3_lab[2], result.p4_lab[2]],
                    [result.p3_lab[3], result.p4_lab[3]]),
            color=:blue, linewidth=2, label="")
    xlabel!("py"); ylabel!("pz"); title!("yz投影")
    
    plot(p1, p2, p3, layout=(1, 3), size=(1200, 400))
end

# 显示物理量信息
function print_kinematics_info(result::ScatteringKinematics)
    println("=" ^ 60)
    println("散射运动学信息")
    println("=" ^ 60)
    println("质心系参数:")
    println("  s = $(result.s) fm⁻² ($(result.s * 197.327^2) MeV²)")
    println("  √s = $(sqrt(result.s)) fm⁻¹ ($(sqrt(result.s) * 197.327) MeV)")
    println("  p* = $(result.p_star) fm⁻¹")
    println("  β = $(norm(result.beta))")
    println("  γ = $(result.gamma)")
    println("\n椭球参数:")
    println("  中心: $(round.(result.ellipsoid.center, digits=4))")
    println("  半轴长: $(round.(result.ellipsoid.half_lengths, digits=4)) fm⁻¹")
    println("\n出射动量:")
    println("  θ* = $(rad2deg(result.theta_star))°, φ* = $(rad2deg(result.phi_star))°")
    println("  p₃ = $(round.(result.p3_lab, digits=4)) fm⁻¹")
    println("  p₄ = $(round.(result.p4_lab, digits=4)) fm⁻¹")
    println("=" ^ 60)
end
```

### 实现优先级

**Phase 1: 后端核心** (必需)
- [ ] `FrameTransformations.jl` - Källen函数、质心系参数、仿射变换
- [ ] `EllipsoidCalculation.jl` - 椭球参数、特征分解、表面采样
- [ ] `MomentumMapping.jl` - 动量映射主逻辑、物理验证
- [ ] `HTTPServer.jl` - REST API、JSON序列化、CORS支持
- [ ] 基础测试: 四动量守恒、质壳条件、boost可逆性

**Phase 2: 前端基础** (必需)
- [ ] `web/index.html` - 页面结构和表单
- [ ] `web/js/api.js` - HTTP通信模块
- [ ] `web/js/visualization.js` - three.js椭球和箭头渲染
- [ ] `web/js/ui.js` - 表单处理和数据更新
- [ ] `web/css/style.css` - 页面布局和样式
- [ ] `server.jl` - 启动脚本

**Phase 3: 交互增强** (按优先级)
- [ ] **P1**: 在椭球上点击选择出射方向
  * 实现射线拾取(Raycaster)
  * 更新p3, p4矢量显示
  * 重新调用API计算物理量
- [ ] **P2**: 实时调整入射动量
  * 添加滑块控件
  * 防抖优化(避免频繁请求)
  * 动态更新可视化
- [ ] **P3**: 洛伦兹变换动画
  * 计算中间帧参数
  * TWEEN.js动画库集成
  * 播放/暂停/重置控制

**Phase 4: 扩展功能** (可选,未来)
- [ ] 并排显示CMS和Lab视图(双画布)
- [ ] 导出功能(截图、参数保存)
- [ ] 与`ScatteringAmplitude`集成显示微分截面

### 单位制与数值精度

**单位制** (与`Constants_PNJL.jl`一致)
- 能量/动量: fm⁻¹ (1 fm⁻¹ = 197.327 MeV)
- 角度: rad
- 转换因子: `ħc_MeV_fm = 197.327`

**数值精度要求**
- 四动量守恒: $|\Delta E| < 10^{-10}$ fm⁻¹, $|\Delta \mathbf{p}| < 10^{-10}$ fm⁻¹
- 质壳条件: $|E^2 - \mathbf{p}^2 - m^2| < 10^{-10}$ fm⁻²
- 阈值判断: $s > (m_3 + m_4)^2 + 10^{-12}$ fm⁻²
- 小速度处理: 当 $\beta^2 < 10^{-14}$ 时令 $A = I$

### 测试策略

**`test_unit/test_frame_transformations.jl`**
- Källen函数数值验证
- 质心系速度和洛伦兹因子计算
- 仿射变换矩阵构造($\beta \to 0$ 极限)
- Boost变换可逆性

**`test_unit/test_ellipsoid_calculation.jl`**
- 特征分解正确性(正交性、正定性)
- 半轴长计算
- 表面采样点验证(所有点在椭球上)
- 退化情况($\beta \to 0$时椭球→球)

**`test_unit/test_momentum_mapping.jl`**
- 四动量守恒(多组随机角度)
- 质壳条件
- 阈值边界测试
- 不同粒子质量组合
- 性能基准: 单次计算 <100μs

### 示例用法

#### 启动Web服务

```bash
# 1. 启动Julia HTTP服务器
julia server.jl

# 输出:
# ============================================================
# 服务器已启动: http://localhost:8080
# API端点: POST http://localhost:8080/compute
# 按 Ctrl+C 停止服务器
# ============================================================

# 2. 打开浏览器访问
# 在浏览器中打开: web/index.html
# 或使用本地HTTP服务器:
python -m http.server 8000 --directory web
# 然后访问: http://localhost:8000
```

#### Julia命令行测试(可选)

```julia
include("src/simulation/MomentumMapping.jl")

# 1. 定义入射态和质量 (uu → uu散射)
p1 = [0.5, 0.0, 1.8]  # fm⁻¹
p2 = [-0.5, 0.0, -1.8]
m_u = 1.52  # fm⁻¹ (约300 MeV)

# 2. 选择质心系方向
theta_star = π/4  # 45度
phi_star = π/6    # 30度

# 3. 计算出射动量
result = calculate_outgoing_momenta(p1, p2, m_u, m_u, m_u, m_u,
                                    theta_star, phi_star)

# 4. 验证物理约束
validate_kinematics(result, m_u, m_u, m_u, m_u)

# 5. 显示信息
print_kinematics_info(result)

# 6. 查看椭球参数
println("椭球中心: ", result.ellipsoid.center)
println("半轴长: ", result.ellipsoid.half_lengths)
```

#### 浏览器使用流程

1. **输入参数**：在左侧面板填写入射动量和粒子质量
2. **计算**：点击"计算并可视化"按钮
3. **查看结果**：
   - 右侧3D视图显示椭球(蓝色半透明)和动量箭头(红色入射/蓝色出射)
   - 左侧显示物理量(s, √s, p*, β, γ)
4. **交互**：
   - 鼠标拖动旋转视角
   - 滚轮缩放
   - 右键平移

### 文件清单

新建文件:
```
src/simulation/
├── FrameTransformations.jl      (~120 lines)
├── EllipsoidCalculation.jl      (~80 lines)
├── MomentumMapping.jl           (~130 lines)
├── HTTPServer.jl                (~100 lines)
└── ScatteringVisualization.jl   (~150 lines, 可选保留用于测试)

test_unit/
├── test_frame_transformations.jl     (~100 lines)
├── test_ellipsoid_calculation.jl     (~80 lines)
└── test_momentum_mapping.jl          (~120 lines)

web/
├── index.html                    (~80 lines)
├── css/
│   └── style.css                (~100 lines)
└── js/
    ├── api.js                   (~40 lines)
    ├── visualization.js         (~150 lines)
    └── ui.js                    (~80 lines)

server.jl                         (~30 lines)
examples/
└── web_demo.md                   (~50 lines, 使用说明)
```

预计总代码量: 
- Julia后端: ~730 lines
- JavaScript前端: ~450 lines
- HTML/CSS: ~180 lines
- **总计: ~1360 lines**
