"""Optional numeric extensions specialize this direct scalar trait."""
_direct_optional_scalar_type(::Type) = nothing

function _optional_scalar_type(type::Type)
    return _optional_scalar_type(type, Set{Any}())
end

function _optional_scalar_type(type::Type, visited::Set{Any})
    type in visited && return nothing
    push!(visited, type)
    direct = _direct_optional_scalar_type(type)
    direct === nothing || return direct
    unwrapped = Base.unwrap_unionall(type)
    unwrapped isa DataType || return nothing
    for parameter in unwrapped.parameters
        parameter isa Type || continue
        candidate = _optional_scalar_type(parameter, visited)
        candidate === nothing || return candidate
    end
    return nothing
end

_hascomplex_type(::Type{<:Complex}) = true
_hascomplex_type(::Type{<:AbstractArray{T}}) where {T} = _hascomplex_type(T)
_hascomplex_type(::Type{T}) where {T<:Tuple} = any(
    parameter -> parameter isa Type && _hascomplex_type(parameter),
    Base.unwrap_unionall(T).parameters,
)
_hascomplex_type(::Type{NamedTuple{Names,T}}) where {Names,T} =
    _hascomplex_type(T)
_hascomplex_type(::Type) = false

"""
    resolve_T(args...)

Resolve the common scalar used by strict materialized constructors. Core
numbers normalize to `Float64`; optional scalar packages participate through
extension methods for `_optional_scalar_type`.
"""
function resolve_T(args...)
    types = map(typeof, args)
    optional = nothing
    for type in types
        optional = _optional_scalar_type(type)
        optional === nothing || break
    end
    scalar = optional === nothing ? BASE_FLOAT : optional
    return any(_hascomplex_type, types) ? Complex{scalar} : scalar
end

function _coerce_elt_to_T end
_coerce_elt_to_T(value::Number, ::Type{T}) where {T<:AbstractFloat} =
    convert(T, value)
_coerce_elt_to_T(value::Bool, ::Type{T}) where {T<:AbstractFloat} = value
_coerce_elt_to_T(::Nothing, ::Type) = nothing
_coerce_elt_to_T(::Missing, ::Type) = missing
_coerce_elt_to_T(value::Union{Bool,Symbol,String,Function,DataType}, ::Type) = value
_coerce_elt_to_T(value, ::Type) = value

function coerce_to_T end
coerce_to_T(value::T, ::Type{T}) where {T} = value
coerce_to_T(value::Real, ::Type{C}) where {P,C<:Complex{P}} =
    C(coerce_to_T(value, P))
coerce_to_T(value::Complex{P}, ::Type{Complex{P}}) where {P} = value
function coerce_to_T(value::Complex, ::Type{Complex{P}}) where {P}
    return Complex{P}(coerce_to_T(real(value), P), coerce_to_T(imag(value), P))
end
coerce_to_T(value::Complex, ::Type{T}) where {T<:Real} =
    coerce_to_T(real(value), T)
coerce_to_T(value::Number, ::Type{T}) where {T} = _coerce_elt_to_T(value, T)
coerce_to_T(values::AbstractArray{T}, ::Type{T}) where {T} = values
coerce_to_T(values::AbstractArray, ::Type{T}) where {T} =
    map(value -> coerce_to_T(value, T), values)
coerce_to_T(values::Tuple{Vararg{T}}, ::Type{T}) where {T} = values
coerce_to_T(values::Tuple, ::Type{T}) where {T} =
    map(value -> coerce_to_T(value, T), values)
function coerce_to_T(tuple_value::NamedTuple{Names}, ::Type{T}) where {Names,T}
    return NamedTuple{Names}(
        map(value -> coerce_to_T(value, T), values(tuple_value)),
    )
end
coerce_to_T(value, ::Type{T}) where {T} = _coerce_elt_to_T(value, T)
