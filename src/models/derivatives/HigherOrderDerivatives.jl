"""
    HigherOrderDerivatives

高阶导公共工具：
- `nth_derivative(f, x, n)`
- `susceptibility_scale(T_fm, n)`
"""
module HigherOrderDerivatives

using ForwardDiff

export nth_derivative, susceptibility_scale

function nth_derivative(f, x, n::Int)
    n >= 1 || throw(ArgumentError("n must be >= 1, got $n"))
    n == 1 && return ForwardDiff.derivative(f, x)
    inner = y -> nth_derivative(f, y, n - 1)
    return ForwardDiff.derivative(inner, x)
end

function susceptibility_scale(T_fm::Real, n::Int)
    T_fm > 0 || throw(ArgumentError("T_fm must be positive, got $T_fm"))
    n >= 1 || throw(ArgumentError("n must be >= 1, got $n"))
    return T_fm^(n - 4)
end

end # module HigherOrderDerivatives
