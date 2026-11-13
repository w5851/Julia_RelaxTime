module TaylorDiffForwardDiffCompat

using ForwardDiff
using TaylorDiff
import ..PNJLIntegrals

import Base: +, -, *, /, promote_rule

const TDScalar = TaylorDiff.TaylorScalar
const FDDual = ForwardDiff.Dual

# A converged pressure value does not certify its Taylor coefficients. Include
# every coefficient when the radial quadrature chooses panels for a response.
@inline PNJLIntegrals._error_norm(x::TDScalar) = max(
    PNJLIntegrals._error_norm(TaylorDiff.value(x)),
    maximum(PNJLIntegrals._error_norm, TaylorDiff.partials(x)))

# The PNJL Taylor-series gap path evaluates ForwardDiff gradients whose scalar
# values are TaylorDiff polynomials. Keep this narrow interop shim local to the
# Models derivative stack.
ForwardDiff.can_dual(::Type{<:TDScalar}) = true

@inline (::Type{TDScalar{T, P}})(x::TDScalar{T, P}) where {T, P} = x

@inline promote_rule(::Type{TDScalar{TV, P}}, ::Type{FDDual{Tag, V, N}}) where {TV, P, Tag, V, N} =
    FDDual{Tag, promote_type(TDScalar{TV, P}, V), N}

@inline +(a::TDScalar, b::FDDual{Tag, V, N}) where {Tag, V, N} =
    FDDual{Tag}(a + ForwardDiff.value(b), ForwardDiff.partials(b))

@inline +(a::FDDual{Tag, V, N}, b::TDScalar) where {Tag, V, N} =
    FDDual{Tag}(ForwardDiff.value(a) + b, ForwardDiff.partials(a))

@inline -(a::TDScalar, b::FDDual{Tag, V, N}) where {Tag, V, N} =
    FDDual{Tag}(a - ForwardDiff.value(b), -ForwardDiff.partials(b))

@inline -(a::FDDual{Tag, V, N}, b::TDScalar) where {Tag, V, N} =
    FDDual{Tag}(ForwardDiff.value(a) - b, ForwardDiff.partials(a))

@inline *(a::TDScalar, b::FDDual{Tag, V, N}) where {Tag, V, N} =
    FDDual{Tag}(a * ForwardDiff.value(b), a * ForwardDiff.partials(b))

@inline *(a::FDDual{Tag, V, N}, b::TDScalar) where {Tag, V, N} =
    FDDual{Tag}(ForwardDiff.value(a) * b, ForwardDiff.partials(a) * b)

@inline function /(a::TDScalar, b::FDDual{Tag, V, N}) where {Tag, V, N}
    bval = ForwardDiff.value(b)
    return FDDual{Tag}(a / bval, ForwardDiff.partials(b) * (-a / (bval * bval)))
end

@inline /(a::FDDual{Tag, V, N}, b::TDScalar) where {Tag, V, N} =
    FDDual{Tag}(ForwardDiff.value(a) / b, ForwardDiff.partials(a) / b)

end # module TaylorDiffForwardDiffCompat
