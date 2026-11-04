# Verify that the reference comparison itself is not an allocation hot spot.
include(joinpath(@__DIR__,"..","..","scripts","analysis","relaxtime","ThermalWeightedSampling.jl"))
include(joinpath(@__DIR__,"..","..","scripts","analysis","relaxtime","StreamedThermalValidation.jl"))
using .ThermalWeightedSampling, .StreamedThermalValidation

function boxed_comparison(candidate,reference,previous)
    query=let candidate_table=candidate,reference_table=reference,previous_table=previous
        s->begin
            c=candidate_table(s);r=reference_table(s);p=previous_table(s)
            (c,r,abs(c-r),abs(p-r))
        end
    end
    c=r=0.0
    return query
end

function query_loop(query,n)
    result=0.0
    for i in 1:n
        result+=query(1.0+2.0*i/(n+1))[3]
    end
    return result
end

a=WeightedSigmaTable([1.0,2.0,3.0],[0.3,0.5,0.2],0.0,0.0)
b=WeightedSigmaTable([1.0,2.0,3.0],[0.31,0.48,0.21],0.0,0.0)
bad=boxed_comparison(a,b,a)
good=TableComparison(a,b,a)
query_loop(bad,10);query_loop(good,10)
old=@timed query_loop(bad,100_000)
new=@timed query_loop(good,100_000)
old.value==new.value || error("comparison changed its numerical result")
println((queries=100_000,old_seconds=old.time,new_seconds=new.time,
    old_allocated_bytes=old.bytes,new_allocated_bytes=new.bytes,
    old_capture_types=fieldtypes(typeof(bad)),new_capture_types=fieldtypes(typeof(good))))
