module ThermalAngularQuadrature

export AngularGroup, angular_groups

struct AngularGroup
    pairs::Vector{Tuple{Int,Int}}
    cosines::Vector{Float64}
    weights::Vector{Float64}
end

"""Combine identical relative-angle arguments without assuming occupation symmetry.

`cos_relative[k,a,b]` must be the arguments computed by the original quadrature.
Dictionary keys use exact `isequal` equality, including the sign of zero. There
is no angular tolerance, interpolation, or change of the quadrature nodes.
The result is local to the caller, so custom grids cannot reuse stale geometry.
"""
function angular_groups(cos_relative::Array{Float64,3}, phi_weights::Vector{Float64})
    size(cos_relative,1)==length(phi_weights) ||
        throw(ArgumentError("azimuthal geometry and weight lengths differ"))
    all(isfinite,cos_relative) && all(isfinite,phi_weights) ||
        throw(ArgumentError("angular geometry and weights must be finite"))
    groups=AngularGroup[]
    keys=Dict{Tuple,Int}()
    for a in axes(cos_relative,2), b in axes(cos_relative,3)
        key=Tuple(cos_relative[:,a,b])
        index=get(keys,key,0)
        if index==0
            cosines=Float64[]; weights=Float64[]; phi_keys=Dict{Float64,Int}()
            for k in eachindex(phi_weights)
                c=cos_relative[k,a,b]
                p=get(phi_keys,c,0)
                if p==0
                    push!(cosines,c); push!(weights,phi_weights[k])
                    phi_keys[c]=length(cosines)
                else
                    weights[p]+=phi_weights[k]
                end
            end
            push!(groups,AngularGroup([(a,b)],cosines,weights))
            keys[key]=length(groups)
        else
            push!(groups[index].pairs,(a,b))
        end
    end
    return groups
end

end
