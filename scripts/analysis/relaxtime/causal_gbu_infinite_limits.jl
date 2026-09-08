"""Analytic unitary phase limits; no finite offset is mistaken for the limit."""
module CausalGBUInfiniteLimits

"""For rho(U+u)=A*sqrt(u)+o(sqrt(u)), F(U) is real and nonzero.

The sign of A sets the side of approach. The independent root count is NOT
an input. A supplied numerical error budget protects the real-inverse sign.
"""
function threshold_limit(inverse_at_threshold,A;inverse_error_budget=1e-6)
    all(isfinite,(inverse_at_threshold,A,inverse_error_budget)) && inverse_error_budget>=0 ||
        throw(ArgumentError("finite threshold data and nonnegative error budget required"))
    A!=0 || throw(ArgumentError("vanishing leading cut requires a higher-order limit"))
    abs(imag(inverse_at_threshold))<=inverse_error_budget || throw(ArgumentError("threshold inverse must be real"))
    f=real(inverse_at_threshold)
    abs(f)>inverse_error_budget || throw(ArgumentError("Mott threshold inside the kernel error budget"))
    return (phase=f>0 ? 0. : copysign(Float64(pi),A),inverse_margin=abs(f)-inverse_error_budget,
        leading_cut_sign=sign(A),count_used=false)
end

"""Keep an argument-principle rectangle inside a gap without excluding a found root."""
function gap_window(a,b,roots;maximum_margin=1e-7)
    all(isfinite,(a,b,maximum_margin)) && a<b && 0<maximum_margin<(b-a)/2 &&
        all(r->isfinite(r) && a<r<b,roots) || throw(ArgumentError("ordered gap and strictly interior roots required"))
    ml=min(maximum_margin,minimum((r-a)/8 for r in roots;init=maximum_margin))
    mr=min(maximum_margin,minimum((b-r)/8 for r in roots;init=maximum_margin))
    left,right=a+ml,b-mr
    a<left<right<b && all(r->left<r<right,roots) ||
        throw(ArgumentError("gap/root separation unresolved at this precision"))
    return (left=left,right=right,left_margin=ml,right_margin=mr)
end
end
