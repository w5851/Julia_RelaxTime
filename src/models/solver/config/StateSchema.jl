"""
    ModelStateSchema

模型状态向量的 schema 描述：
- `model_kind`：模型类型符号
- `fields`：状态字段顺序（flatten/unflatten 的唯一依据）
"""
struct ModelStateSchema{F<:Tuple}
    model_kind::Symbol
    fields::F
end

@inline state_dim(schema::ModelStateSchema) = length(schema.fields)

@inline function _generic_state_fields(dim::Int)
    dim > 0 || throw(ArgumentError("state dim must be positive, got $dim"))
    return ntuple(i -> Symbol("x", i), dim)
end

@inline function _default_state_fields(model_kind::Symbol, dim::Int)
    if model_kind === :NJL
        return (:phi_u, :phi_d, :phi_s)
    elseif model_kind === :NJL2
        return (:phi_u, :phi_d)
    elseif model_kind === :PNJL || model_kind === :PNJLMagnetic || model_kind === :RPNJL
        return (:phi_u, :phi_d, :phi_s, :Phi, :PhiBar)
    elseif model_kind === :Rotation
        return (:phi_u, :phi_d, :phi_s)
    else
        return _generic_state_fields(dim)
    end
end

"""
    schema_for_model(model_kind::Symbol) -> ModelStateSchema

按模型类型构造状态 schema。
"""
function schema_for_model(model_kind::Symbol)
    # The magnetic route requires a caller-supplied positive eB, so schema
    # construction must not instantiate its model with an invalid zero-field
    # default merely to discover the already-defined five-dimensional shape.
    dim = model_kind === :PNJLMagnetic ? 5 : gap_state_dim(create_model(model_kind))
    fields = _default_state_fields(model_kind, dim)
    length(fields) == dim || throw(ArgumentError("schema fields length mismatch for $model_kind: expected $dim, got $(length(fields))"))
    return ModelStateSchema(model_kind, fields)
end

"""
    flatten_state(schema, st::NamedTuple) -> Vector

按 schema 字段顺序展开状态命名元组。
"""
function flatten_state(schema::ModelStateSchema, st::NamedTuple)
    n = state_dim(schema)
    values = ntuple(n) do i
        field = schema.fields[i]
        haskey(st, field) || throw(ArgumentError("missing state field :$field for model $(schema.model_kind)"))
        st[field]
    end
    T = promote_type((typeof(v) for v in values)...)
    x = Vector{T}(undef, n)
    @inbounds for i in 1:n
        x[i] = convert(T, values[i])
    end
    return x
end

"""
    unflatten_state(schema, x::AbstractVector) -> NamedTuple

将向量按 schema 字段顺序还原为命名元组。
"""
function unflatten_state(schema::ModelStateSchema, x::AbstractVector)
    n = state_dim(schema)
    length(x) == n || throw(ArgumentError("state vector length mismatch for $(schema.model_kind): expected $n, got $(length(x))"))
    return NamedTuple{schema.fields}(Tuple(x[i] for i in 1:n))
end
