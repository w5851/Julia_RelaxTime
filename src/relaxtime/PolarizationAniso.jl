"""
# Polarization.jl
极化函数的计算，适用于PNJL模型中的赝标量和标量介子。

需要的函数包括B0和A函数的计算,它们的实现见 `OneLoopIntegrals.jl` 文件。
其中A函数的值A_value可以作为输入直接传入以提高效率,这是因为A函数只依赖于质量、化学势、温度和Polyakov环,与动量无关,在多次调用极化函数时可以复用。
"""
module PolarizationAniso

include("../Constants_PNJL.jl")
include("../relaxtime/OneLoopIntegrals.jl")
include("../relaxtime/OneLoopIntegralsAniso.jl")
using .Constants_PNJL: N_color
using .OneLoopIntegrals: B0, EPS_SEGMENT
using .OneLoopIntegralsAniso: B0_correction

"""
    polarization(channel, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1_value, A2_value, num_s_quark)
计算PNJL模型中两粒子碰撞下的赝标量(P)或标量(S)介子(通过channel标记区分)的极化函数(动量各向异性下)
返回 polarization_real, polarization_imag
"""
function polarization_aniso(channel::Symbol, k0::Float64, k_norm::Float64, m1::Float64, m2::Float64, 
                                  μ1::Float64, μ2::Float64, T::Float64, Φ::Float64, Φbar::Float64, 
                                  ξ::Float64, A1_value::Float64, A2_value::Float64, num_s_quark::Int)
    # 计算总系数
    factor = -N_color / (8π^2)
    # 计算能量-化学势组合参数
    λ = k0 + μ1 - μ2

    # 计算B0函数项
    B0_real, B0_imag = B0(λ, k_norm, m1, μ1, m2, μ2, T; Φ=Φ, Φbar=Φbar)
    # 当num_s_quark为1时,要恢复B0对k0的对称性
    if num_s_quark == 1
        λ_extra = -k0 + μ1 - μ2
        B0_real_extra, B0_imag_extra = B0(λ_extra, k_norm, m1, μ1, m2, μ2, T; Φ=Φ, Φbar=Φbar)
        B0_real = 0.5 * (B0_real + B0_real_extra)
        B0_imag = 0.5 * (B0_imag + B0_imag_extra) 
    end

    # 计算B0函数项的修正
    if ξ>EPS_SEGMENT
        B0_corr_real, B0_corr_imag = B0_correction(λ, k_norm, m1, μ1, m2, μ2, T, ξ; Φ=Φ, Φbar=Φbar)
        # 当num_s_quark为1时,要恢复B0对k0的对称性
        if num_s_quark == 1
            B0_corr_real_extra, B0_corr_imag_extra = B0_correction(λ_extra, k_norm, m1, μ1, m2, μ2, T, ξ; Φ=Φ, Φbar=Φbar)
            B0_corr_real = 0.5 * (B0_corr_real + B0_corr_real_extra)
            B0_corr_imag = 0.5 * (B0_corr_imag + B0_corr_imag_extra)
        end
        B0_real += B0_corr_real
        B0_imag += B0_corr_imag
    end

    real_part = A1_value + A2_value
    imag_part = 0.0
    prefactor = k_norm^2 - λ^2
    # 根据通道选择极化函数的计算公式
    if channel == :P # 赝标量通道
        prefactor += (m1 - m2)^2
    elseif channel == :S # 标量通道
        prefactor += (m1 + m2)^2
    else
        throw(ArgumentError("Unsupported channel: $channel. Use :P or :S."))
    end
    real_part += prefactor * B0_real
    imag_part += prefactor * B0_imag
    return factor * real_part, factor * imag_part
end

end # module Polarization