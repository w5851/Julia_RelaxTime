"""
    VarSchema

变量 schema：定义 x/theta 的键顺序，用于 named <-> vector 双向转换。
"""
struct VarSchema{NX, NT}
    x_keys::NTuple{NX,Symbol}
    theta_keys::NTuple{NT,Symbol}
end

"""
    SchemaRegistry

按 `(model_kind, spec_tag)` 注册变量 schema。
"""
struct SchemaRegistry
    table::Dict{Tuple{Symbol,Symbol},VarSchema}
end

SchemaRegistry() = SchemaRegistry(Dict{Tuple{Symbol,Symbol},VarSchema}())

@inline function _assert_unique_symbols(keys::Tuple, label::Symbol)
    seen = Set{Symbol}()
    for k in keys
        k isa Symbol || throw(ArgumentError("$label keys must be Symbol, got $(typeof(k))"))
        k in seen && throw(ArgumentError("$label keys contain duplicates: :$k"))
        push!(seen, k)
    end
    return nothing
end

function validate_schema(schema::VarSchema; x_dim::Union{Nothing,Int}=nothing, theta_dim::Union{Nothing,Int}=nothing)
    _assert_unique_symbols(schema.x_keys, :x)
    _assert_unique_symbols(schema.theta_keys, :theta)

    if x_dim !== nothing
        length(schema.x_keys) == x_dim || throw(ArgumentError("x_dim mismatch: expected $x_dim, got $(length(schema.x_keys))"))
    end
    if theta_dim !== nothing
        length(schema.theta_keys) == theta_dim || throw(ArgumentError("theta_dim mismatch: expected $theta_dim, got $(length(schema.theta_keys))"))
    end
    return nothing
end

function register_schema!(reg::SchemaRegistry, model_kind::Symbol, spec_tag::Symbol, schema::VarSchema)
    validate_schema(schema)
    reg.table[(model_kind, spec_tag)] = schema
    return reg
end

function schema_for(reg::SchemaRegistry, model_kind::Symbol, spec_tag::Symbol)
    key = (model_kind, spec_tag)
    haskey(reg.table, key) || throw(ArgumentError("missing schema for $(key)"))
    return reg.table[key]
end

@inline function _keys_for_role(schema::VarSchema, role::Symbol)
    if role === :x
        return schema.x_keys
    elseif role === :theta
        return schema.theta_keys
    end
    throw(ArgumentError("role must be :x or :theta, got $(role)"))
end

function named_to_vec(named::NamedTuple, schema::VarSchema, role::Symbol)
    keys = _keys_for_role(schema, role)
    values = ntuple(i -> begin
        key = keys[i]
        haskey(named, key) || throw(ArgumentError("missing key :$key for role $(role)"))
        named[key]
    end, length(keys))

    T = promote_type(map(typeof, values)...)
    vec = Vector{T}(undef, length(keys))
    @inbounds for i in eachindex(keys)
        vec[i] = convert(T, values[i])
    end
    return vec
end

function vec_to_named(vec::AbstractVector{T}, schema::VarSchema, role::Symbol) where {T}
    keys = _keys_for_role(schema, role)
    length(vec) == length(keys) || throw(ArgumentError("vector length mismatch for role $(role): expected $(length(keys)), got $(length(vec))"))
    return NamedTuple{keys}(Tuple(vec[i] for i in eachindex(keys)))
end
