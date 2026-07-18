"""
    MixedTaylorJets

Small internal multivariate Taylor jet algebra for PNJL mixed conserved-charge
derivatives. Coefficients are stored as Taylor coefficients, i.e.
`c_alpha = derivative_alpha / alpha!`.
"""
module MixedTaylorJets

using ForwardDiff
using ..PNJLCore: PNJLIntegrals

import Base: +, -, *, /, ^, <, <=, >, >=, abs, convert, exp, float, getindex, inv, isfinite
import Base: isless, log, max, min, one, promote_rule, show, sqrt, zero

export MixedTaylorJet, jet_constant, jet_variable, jet_value
export jet_coefficient, jet_extract_derivative, jet_basis, jet_degree_positions

struct MixedTaylorJet{D, N, L} <: Real
    coeffs::NTuple{L, Float64}
end

const FDDual = ForwardDiff.Dual

function _multiindices_exact_degree!(out, current::Vector{Int}, d::Int, remaining::Int)
    if d == length(current)
        current[d] = remaining
        push!(out, Tuple(current))
        return out
    end
    for n in 0:remaining
        current[d] = n
        _multiindices_exact_degree!(out, current, d + 1, remaining - n)
    end
    return out
end

function _make_multiindices(D::Int, N::Int)
    D >= 1 || error("D must be >= 1")
    N >= 0 || error("N must be >= 0")
    out = Tuple{Vararg{Int}}[]
    current = zeros(Int, D)
    for degree in 0:N
        _multiindices_exact_degree!(out, current, 1, degree)
    end
    return out
end

_basis_length(D::Int, N::Int) = binomial(D + N, D)

@generated function _jet_type(::Val{D}, ::Val{N}) where {D, N}
    L = _basis_length(D, N)
    return :(MixedTaylorJet{$D, $N, $L})
end

@generated function jet_basis(::Type{MixedTaylorJet{D, N, L}}) where {D, N, L}
    basis = Tuple(_make_multiindices(D, N))
    length(basis) == L || error("invalid MixedTaylorJet basis length")
    return QuoteNode(basis)
end

@generated function jet_degree_positions(::Type{MixedTaylorJet{D, N, L}}) where {D, N, L}
    basis = _make_multiindices(D, N)
    levels = map(degree -> Tuple(findall(idx -> sum(basis[idx]) == degree, eachindex(basis))), 0:N)
    return QuoteNode(Tuple(levels))
end

@generated function _mul_targets(::Type{MixedTaylorJet{D, N, L}}) where {D, N, L}
    basis = _make_multiindices(D, N)
    pos = Dict{Tuple{Vararg{Int}}, Int}(basis[i] => i for i in eachindex(basis))
    targets = Int[]
    for α in basis
        for β in basis
            γ = ntuple(i -> α[i] + β[i], D)
            push!(targets, sum(γ) <= N ? pos[γ] : 0)
        end
    end
    return QuoteNode(Tuple(targets))
end

@inline function _position(::Type{MixedTaylorJet{D, N, L}}, idx::NTuple{D, Int}) where {D, N, L}
    sum(idx) <= N || return 0
    basis = jet_basis(MixedTaylorJet{D, N, L})
    @inbounds for i in 1:L
        basis[i] == idx && return i
    end
    return 0
end

@inline function _factorial_product(idx)
    prod = 1
    @inbounds for n in idx
        prod *= factorial(n)
    end
    return prod
end

@inline _tuple_from_vec(v::AbstractVector{<:Real}, ::Val{L}) where {L} =
    ntuple(i -> Float64(v[i]), Val(L))

@inline function _constant_tuple(x::Float64, ::Val{L}) where {L}
    return ntuple(i -> i == 1 ? x : 0.0, Val(L))
end

@inline function _variable_tuple(x::Float64, variable::Int, ::Val{D}, ::Val{N}, ::Val{L}) where {D, N, L}
    1 <= variable <= D || throw(ArgumentError("variable must be in 1:$D, got $variable"))
    e = ntuple(i -> i == variable ? 1 : 0, D)
    pos = _position(MixedTaylorJet{D, N, L}, e)
    pos != 0 || throw(ArgumentError("cannot create variable for order N=$N"))
    return ntuple(i -> i == 1 ? x : (i == pos ? 1.0 : 0.0), Val(L))
end

@generated function _jet_constant(x::Float64, ::Val{D}, ::Val{N}) where {D, N}
    L = _basis_length(D, N)
    return :(MixedTaylorJet{$D, $N, $L}(_constant_tuple(x, Val($L))))
end

@generated function _jet_variable(x::Float64, variable::Int, ::Val{D}, ::Val{N}) where {D, N}
    L = _basis_length(D, N)
    return :(MixedTaylorJet{$D, $N, $L}(_variable_tuple(x, variable, Val($D), Val($N), Val($L))))
end

@inline function jet_constant(x::Real, D::Int, N::Int)
    D >= 1 || throw(ArgumentError("D must be >= 1, got $D"))
    N >= 1 || throw(ArgumentError("N must be >= 1, got $N"))
    return _jet_constant(Float64(x), Val(D), Val(N))
end

@inline function jet_variable(x::Real, variable::Int, D::Int, N::Int)
    D >= 1 || throw(ArgumentError("D must be >= 1, got $D"))
    N >= 1 || throw(ArgumentError("N must be >= 1, got $N"))
    return _jet_variable(Float64(x), variable, Val(D), Val(N))
end

@inline jet_value(x::MixedTaylorJet) = x.coeffs[1]
@inline PNJLIntegrals._primal_float(x::MixedTaylorJet) = jet_value(x)
@inline getindex(x::MixedTaylorJet, i::Int) = x.coeffs[i]

@inline function jet_coefficient(x::MixedTaylorJet{D, N, L}, pos::Int) where {D, N, L}
    1 <= pos <= L || throw(ArgumentError("coefficient position must be in 1:$L, got $pos"))
    return x.coeffs[pos]
end

@inline function jet_coefficient(x::MixedTaylorJet{D, N, L}, idx::NTuple{D, Int}) where {D, N, L}
    pos = _position(MixedTaylorJet{D, N, L}, idx)
    return pos == 0 ? 0.0 : x.coeffs[pos]
end

@inline function jet_extract_derivative(x::MixedTaylorJet{D, N, L}, idx::NTuple{D, Int}) where {D, N, L}
    sum(idx) <= N || throw(ArgumentError("requested derivative order $(sum(idx)) exceeds jet order $N"))
    return jet_coefficient(x, idx) * _factorial_product(idx)
end

@inline zero(::Type{MixedTaylorJet{D, N, L}}) where {D, N, L} =
    MixedTaylorJet{D, N, L}(ntuple(_ -> 0.0, Val(L)))
@inline one(::Type{MixedTaylorJet{D, N, L}}) where {D, N, L} =
    MixedTaylorJet{D, N, L}(_constant_tuple(1.0, Val(L)))
@inline zero(x::MixedTaylorJet{D, N, L}) where {D, N, L} = zero(MixedTaylorJet{D, N, L})
@inline one(x::MixedTaylorJet{D, N, L}) where {D, N, L} = one(MixedTaylorJet{D, N, L})

@inline convert(::Type{MixedTaylorJet{D, N, L}}, x::MixedTaylorJet{D, N, L}) where {D, N, L} = x
@inline convert(::Type{MixedTaylorJet{D, N, L}}, x::Real) where {D, N, L} =
    MixedTaylorJet{D, N, L}(_constant_tuple(Float64(x), Val(L)))
@inline (::Type{MixedTaylorJet{D, N, L}})(x::Real) where {D, N, L} =
    convert(MixedTaylorJet{D, N, L}, x)

@inline promote_rule(::Type{MixedTaylorJet{D, N, L}}, ::Type{<:Real}) where {D, N, L} =
    MixedTaylorJet{D, N, L}
@inline promote_rule(::Type{MixedTaylorJet{D, N, L}}, ::Type{MixedTaylorJet{D, N, L}}) where {D, N, L} =
    MixedTaylorJet{D, N, L}

@inline +(a::MixedTaylorJet{D, N, L}, b::MixedTaylorJet{D, N, L}) where {D, N, L} =
    MixedTaylorJet{D, N, L}(ntuple(i -> a.coeffs[i] + b.coeffs[i], Val(L)))
@inline +(a::MixedTaylorJet, b::Real) = a + convert(typeof(a), b)
@inline +(a::Real, b::MixedTaylorJet) = convert(typeof(b), a) + b

@inline -(a::MixedTaylorJet{D, N, L}, b::MixedTaylorJet{D, N, L}) where {D, N, L} =
    MixedTaylorJet{D, N, L}(ntuple(i -> a.coeffs[i] - b.coeffs[i], Val(L)))
@inline -(a::MixedTaylorJet, b::Real) = a - convert(typeof(a), b)
@inline -(a::Real, b::MixedTaylorJet) = convert(typeof(b), a) - b
@inline -(a::MixedTaylorJet{D, N, L}) where {D, N, L} =
    MixedTaylorJet{D, N, L}(ntuple(i -> -a.coeffs[i], Val(L)))

function *(a::MixedTaylorJet{D, N, L}, b::MixedTaylorJet{D, N, L}) where {D, N, L}
    acc = zeros(Float64, L)
    targets = _mul_targets(MixedTaylorJet{D, N, L})
    @inbounds for i in 1:L
        ai = a.coeffs[i]
        iszero(ai) && continue
        row = (i - 1) * L
        for j in 1:L
            target = targets[row + j]
            target == 0 && continue
            bj = b.coeffs[j]
            iszero(bj) && continue
            acc[target] += ai * bj
        end
    end
    return MixedTaylorJet{D, N, L}(_tuple_from_vec(acc, Val(L)))
end
@inline *(a::MixedTaylorJet, b::Real) = a * convert(typeof(a), b)
@inline *(a::Real, b::MixedTaylorJet) = convert(typeof(b), a) * b

function inv(a::MixedTaylorJet{D, N, L}) where {D, N, L}
    c0 = jet_value(a)
    iszero(c0) && throw(DomainError(c0, "cannot invert a MixedTaylorJet with zero constant coefficient"))
    r = (a - c0) / c0
    term = one(a)
    total = one(a)
    for _ in 1:N
        term = -term * r
        total = total + term
    end
    return total / c0
end

@inline /(a::MixedTaylorJet, b::MixedTaylorJet) = a * inv(b)
@inline /(a::MixedTaylorJet{D, N, L}, b::Real) where {D, N, L} =
    MixedTaylorJet{D, N, L}(ntuple(i -> a.coeffs[i] / Float64(b), Val(L)))
@inline /(a::Real, b::MixedTaylorJet) = convert(typeof(b), a) / b

function ^(a::MixedTaylorJet, n::Integer)
    n == 0 && return one(a)
    n < 0 && return inv(a ^ (-n))
    result = one(a)
    base = a
    k = n
    while k > 0
        if isodd(k)
            result = result * base
        end
        k >>>= 1
        k > 0 && (base = base * base)
    end
    return result
end

@inline function abs(a::MixedTaylorJet)
    return jet_value(a) < 0 ? -a : a
end

function exp(a::MixedTaylorJet{D, N, L}) where {D, N, L}
    c0 = jet_value(a)
    r = a - c0
    term = one(a)
    total = one(a)
    for k in 1:N
        term = term * r / k
        total = total + term
    end
    return Base.exp(c0) * total
end

function log(a::MixedTaylorJet{D, N, L}) where {D, N, L}
    c0 = jet_value(a)
    c0 > 0 || throw(DomainError(c0, "log requires a positive constant coefficient"))
    r = (a - c0) / c0
    term = r
    total = zero(a)
    for k in 1:N
        total = isodd(k) ? total + term / k : total - term / k
        term = term * r
    end
    return Base.log(c0) + total
end

function sqrt(a::MixedTaylorJet{D, N, L}) where {D, N, L}
    c0 = jet_value(a)
    c0 > 0 || throw(DomainError(c0, "sqrt requires a positive constant coefficient"))
    r = (a - c0) / c0
    total = one(a)
    term = one(a)
    coeff = 1.0
    for k in 1:N
        coeff *= (0.5 - (k - 1)) / k
        term = term * r
        total = total + coeff * term
    end
    return Base.sqrt(c0) * total
end

@inline isfinite(a::MixedTaylorJet) = all(isfinite, a.coeffs)
@inline float(a::MixedTaylorJet) = a
@inline isless(a::MixedTaylorJet, b::MixedTaylorJet) = jet_value(a) < jet_value(b)
@inline isless(a::MixedTaylorJet, b::Real) = jet_value(a) < Float64(b)
@inline isless(a::Real, b::MixedTaylorJet) = Float64(a) < jet_value(b)
@inline <(a::MixedTaylorJet, b::MixedTaylorJet) = isless(a, b)
@inline <(a::MixedTaylorJet, b::Real) = isless(a, b)
@inline <(a::Real, b::MixedTaylorJet) = isless(a, b)
@inline <=(a::MixedTaylorJet, b::MixedTaylorJet) = !isless(b, a)
@inline <=(a::MixedTaylorJet, b::Real) = !(b < a)
@inline <=(a::Real, b::MixedTaylorJet) = !(b < a)
@inline >(a::MixedTaylorJet, b::MixedTaylorJet) = isless(b, a)
@inline >(a::MixedTaylorJet, b::Real) = isless(b, a)
@inline >(a::Real, b::MixedTaylorJet) = isless(b, a)
@inline >=(a::MixedTaylorJet, b::MixedTaylorJet) = !isless(a, b)
@inline >=(a::MixedTaylorJet, b::Real) = !(a < b)
@inline >=(a::Real, b::MixedTaylorJet) = !(a < b)
@inline max(a::MixedTaylorJet, b::MixedTaylorJet) = isless(a, b) ? b : a
@inline max(a::MixedTaylorJet, b::Real) = max(a, convert(typeof(a), b))
@inline max(a::Real, b::MixedTaylorJet) = max(convert(typeof(b), a), b)
@inline min(a::MixedTaylorJet, b::MixedTaylorJet) = isless(a, b) ? a : b
@inline min(a::MixedTaylorJet, b::Real) = min(a, convert(typeof(a), b))
@inline min(a::Real, b::MixedTaylorJet) = min(convert(typeof(b), a), b)

function show(io::IO, x::MixedTaylorJet{D, N, L}) where {D, N, L}
    print(io, "MixedTaylorJet{$D,$N}(", x.coeffs, ")")
end

# ForwardDiff interop for gradients whose scalar values are MixedTaylorJet.
ForwardDiff.can_dual(::Type{<:MixedTaylorJet}) = true

@inline promote_rule(::Type{MixedTaylorJet{D, N, L}}, ::Type{FDDual{Tag, V, M}}) where {D, N, L, Tag, V, M} =
    FDDual{Tag, promote_type(MixedTaylorJet{D, N, L}, V), M}

@inline +(a::MixedTaylorJet, b::FDDual{Tag, V, M}) where {Tag, V, M} =
    FDDual{Tag}(a + ForwardDiff.value(b), ForwardDiff.partials(b))
@inline +(a::FDDual{Tag, V, M}, b::MixedTaylorJet) where {Tag, V, M} =
    FDDual{Tag}(ForwardDiff.value(a) + b, ForwardDiff.partials(a))

@inline -(a::MixedTaylorJet, b::FDDual{Tag, V, M}) where {Tag, V, M} =
    FDDual{Tag}(a - ForwardDiff.value(b), -ForwardDiff.partials(b))
@inline -(a::FDDual{Tag, V, M}, b::MixedTaylorJet) where {Tag, V, M} =
    FDDual{Tag}(ForwardDiff.value(a) - b, ForwardDiff.partials(a))

@inline *(a::MixedTaylorJet, b::FDDual{Tag, V, M}) where {Tag, V, M} =
    FDDual{Tag}(a * ForwardDiff.value(b), a * ForwardDiff.partials(b))
@inline *(a::FDDual{Tag, V, M}, b::MixedTaylorJet) where {Tag, V, M} =
    FDDual{Tag}(ForwardDiff.value(a) * b, ForwardDiff.partials(a) * b)

@inline function /(a::MixedTaylorJet, b::FDDual{Tag, V, M}) where {Tag, V, M}
    bval = ForwardDiff.value(b)
    return FDDual{Tag}(a / bval, ForwardDiff.partials(b) * (-a / (bval * bval)))
end
@inline /(a::FDDual{Tag, V, M}, b::MixedTaylorJet) where {Tag, V, M} =
    FDDual{Tag}(ForwardDiff.value(a) / b, ForwardDiff.partials(a) / b)

end # module MixedTaylorJets
