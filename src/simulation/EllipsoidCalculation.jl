"""
    EllipsoidCalculation

椭球几何计算模块，用于计算出射粒子动量在实验室系中的可行域椭球。

提供功能：
- 椭球参数计算（特征分解）
- 主轴方向和半轴长度
- 椭球表面均匀采样
"""
module EllipsoidCalculation

using LinearAlgebra

export EllipsoidParams, calculate_ellipsoid_parameters, sample_ellipsoid_surface

"""
    EllipsoidParams

椭球参数结构体。

# 字段
- `center::Vector{Float64}`: 椭球中心坐标 [fm⁻¹]
- `axes_directions::Matrix{Float64}`: 主轴方向矩阵（列向量为主轴）
- `half_lengths::Vector{Float64}`: 三个主轴的半轴长 [fm⁻¹]
"""
struct EllipsoidParams
    center::Vector{Float64}            # 中心坐标 b
    axes_directions::Matrix{Float64}   # 主轴方向 (列向量)
    half_lengths::Vector{Float64}      # 半轴长
    
    # 内部构造函数：验证数据有效性
    function EllipsoidParams(center::Vector{Float64}, 
                            axes_directions::Matrix{Float64}, 
                            half_lengths::Vector{Float64})
        @assert length(center) == 3 "Center must be 3-dimensional"
        @assert size(axes_directions) == (3, 3) "Axes directions must be 3×3 matrix"
        @assert length(half_lengths) == 3 "Must have 3 half-lengths"
        @assert all(half_lengths .>= 0) "Half-lengths must be non-negative"
        
        new(center, axes_directions, half_lengths)
    end
end

"""
    calculate_ellipsoid_parameters(A, p_star; center_offset=zeros(3))

计算椭球参数，通过对形状矩阵S = A*A'进行特征分解。

# 参数
- `A`: 仿射变换矩阵（3×3）
- `p_star`: 质心系动量模（椭球在质心系中为半径p*的球）[fm⁻¹]
- `center_offset`: 椭球中心偏移（默认为原点）[fm⁻¹]

# 返回
- `EllipsoidParams`: 椭球参数结构体

# 数学原理
在质心系中，出射粒子3的三动量可行域为球面：|p*| = p_star
经过仿射变换 p_lab = A*p_cms + b，球面变为椭球面。
形状矩阵 S = A*A' 的特征分解给出主轴方向和缩放因子。

# 示例
```julia
ellipsoid = calculate_ellipsoid_parameters(A, cms.p_star, center_offset=b)
```
"""
function calculate_ellipsoid_parameters(A::Matrix{<:Real}, p_star::Real; 
                                       center_offset::Vector{<:Real}=zeros(3))
    # 计算形状矩阵 S = A * A'
    S = A * A'
    
    # 特征分解（使用Symmetric确保对称性）
    eigen_result = eigen(Symmetric(S))
    evals = eigen_result.values
    evecs = eigen_result.vectors
    
    # 检查特征值非负（数值误差可能导致小负值）
    if any(evals .< -1e-10)
        @warn "Negative eigenvalues detected" evals
    end
    evals = max.(evals, 0.0)  # 修正小负值
    
    # 半轴长 = p_star * sqrt(特征值)
    half_lengths = p_star * sqrt.(evals)
    
    # 确保中心偏移是正确维度
    if length(center_offset) != 3
        error("center_offset must be 3-dimensional")
    end
    center = Vector{Float64}(center_offset)
    
    return EllipsoidParams(center, evecs, half_lengths)
end

"""
    sample_ellipsoid_surface(N, p_star, A, b)

在椭球表面均匀采样点（通过质心系角度均匀采样实现）。

# 参数
- `N::Int`: 采样点数量
- `p_star`: 质心系动量模 [fm⁻¹]
- `A`: 仿射变换矩阵（3×3）
- `b`: 偏移向量（3-vector）[fm⁻¹]

# 返回
- `Vector{Vector{Float64}}`: N个采样点的坐标数组

# 原理
在质心系中均匀采样角度(θ, φ)，然后通过仿射变换映射到实验室系。
均匀采样：cosθ ~ U[-1,1], φ ~ U[0,2π]

# 示例
```julia
points = sample_ellipsoid_surface(200, cms.p_star, A, b)
```
"""
function sample_ellipsoid_surface(N::Int, p_star::Real, 
                                  A::Matrix{<:Real}, b::Vector{<:Real})
    @assert N > 0 "Number of samples must be positive"
    
    points = Vector{Vector{Float64}}(undef, N)
    
    for i in 1:N
        # 均匀采样：cosθ ∈ [-1, 1]
        cos_theta = 2 * rand() - 1
        theta = acos(cos_theta)
        
        # 均匀采样：φ ∈ [0, 2π]
        phi = 2π * rand()
        
        # 质心系单位向量
        p_cms = p_star * [sin(theta) * cos(phi), 
                         sin(theta) * sin(phi), 
                         cos(theta)]
        
        # Boost到实验室系
        points[i] = A * p_cms + b
    end
    
    return points
end

"""
    sample_ellipsoid_grid(n_theta, n_phi, p_star, A, b)

在椭球表面生成规则网格点（用于绘制网格线）。

# 参数
- `n_theta`: θ方向网格数
- `n_phi`: φ方向网格数  
- `p_star`: 质心系动量模 [fm⁻¹]
- `A`: 仿射变换矩阵（3×3）
- `b`: 偏移向量（3-vector）[fm⁻¹]

# 返回
- `Matrix{Vector{Float64}}`: n_theta × n_phi 网格点矩阵

# 示例
```julia
grid = sample_ellipsoid_grid(20, 40, cms.p_star, A, b)
```
"""
function sample_ellipsoid_grid(n_theta::Int, n_phi::Int, 
                               p_star::Real, 
                               A::Matrix{<:Real}, b::Vector{<:Real})
    @assert n_theta > 1 && n_phi > 1 "Grid dimensions must be > 1"
    
    grid = Matrix{Vector{Float64}}(undef, n_theta, n_phi)
    
    for (i, theta) in enumerate(range(0, π, length=n_theta))
        for (j, phi) in enumerate(range(0, 2π, length=n_phi))
            # 质心系单位向量
            p_cms = p_star * [sin(theta) * cos(phi), 
                             sin(theta) * sin(phi), 
                             cos(theta)]
            
            # Boost到实验室系
            grid[i, j] = A * p_cms + b
        end
    end
    
    return grid
end

"""
    verify_point_on_ellipsoid(point, ellipsoid; tol=1e-8)

验证点是否在椭球表面上。

# 参数
- `point`: 待验证点坐标 [fm⁻¹]
- `ellipsoid`: EllipsoidParams结构体
- `tol`: 容差（默认1e-8）

# 返回
- `Bool`: 点是否在椭球表面上
- `Float64`: 椭球方程值（应接近1）

# 椭球方程
点P在椭球上当且仅当：
(U'(P - center))' * Λ⁻² * (U'(P - center)) = 1
其中U为主轴方向矩阵，Λ为对角半轴长矩阵

# 示例
```julia
on_surface, value = verify_point_on_ellipsoid(point, ellipsoid)
```
"""
function verify_point_on_ellipsoid(point::Vector{<:Real}, 
                                   ellipsoid::EllipsoidParams; 
                                   tol::Real=1e-8)
    # 将点平移到椭球中心坐标系
    p_centered = point - ellipsoid.center
    
    # 转换到主轴坐标系
    p_principal = ellipsoid.axes_directions' * p_centered
    
    # 归一化坐标
    p_normalized = p_principal ./ ellipsoid.half_lengths
    
    # 椭球方程值（应为1）
    eq_value = dot(p_normalized, p_normalized)
    
    on_surface = abs(eq_value - 1.0) < tol
    
    return on_surface, eq_value
end

end # module
