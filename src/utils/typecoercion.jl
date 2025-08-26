
# --- Type detection for Measurement ---
_hasmeas_type(::Type{<:Measurement}) = true
_hasmeas_type(::Type{<:AbstractArray{S}}) where {S} = _hasmeas_type(S)
_hasmeas_type(::Type{<:Tuple{}}) = false
_hasmeas_type(::Type{T}) where {T<:Tuple} =
    any(_hasmeas_type, Base.unwrap_unionall(T).parameters)
_hasmeas_type(::Type{NamedTuple{N,T}}) where {N,T} = _hasmeas_type(T)
_hasmeas_type(T::Union) = _hasmeas_type(T.a) || _hasmeas_type(T.b)
function _hasmeas_type(T::DataType)
    # FIX: Add guard against recursing into Complex, which is self-contained.
    T <: Complex && return false
    isconcretetype(T) && any(_hasmeas_type, fieldtypes(T))
end
_hasmeas_type(::Type) = false

# --- Type detection for Complex ---
_hascomplex_type(::Type{<:Complex}) = true
_hascomplex_type(::Type{<:AbstractArray{S}}) where {S} = _hascomplex_type(S)
_hascomplex_type(::Type{<:Tuple{}}) = false
_hascomplex_type(::Type{T}) where {T<:Tuple} =
    any(_hascomplex_type, Base.unwrap_unionall(T).parameters)
_hascomplex_type(::Type{NamedTuple{N,T}}) where {N,T} = _hascomplex_type(T)
_hascomplex_type(T::Union) = _hascomplex_type(T.a) || _hascomplex_type(T.b)
function _hascomplex_type(T::DataType)
    # FIX: Add guard against recursing into Measurement, which is self-referential
    # and known not to contain Complex types. This prevents StackOverflowError.
    T <: Measurement && return false
    isconcretetype(T) && any(_hascomplex_type, fieldtypes(T))
end
_hascomplex_type(::Type) = false

# --- Public Type Resolver ---
function resolve_T(args...)
    types = map(typeof, args)
    has_meas = any(_hasmeas_type, types)
    has_complex = any(_hascomplex_type, types)

    if has_meas && has_complex
        return Complex{Measurement{BASE_FLOAT}}
    elseif has_meas
        return Measurement{BASE_FLOAT}
    elseif has_complex
        return Complex{BASE_FLOAT}
    else
        return BASE_FLOAT
    end
end

# ---------- Element-wise Coercion (for leaves) ----------
_meas_inner(::Type{Measurement{S}}) where {S} = S
_coerce_elt_to_T(x::Number, ::Type{R}) where {R<:AbstractFloat} = convert(R, x)
_coerce_elt_to_T(x::Number, ::Type{M}) where {M<:Measurement} = zero(M) + x
_coerce_elt_to_T(m::Measurement, ::Type{M}) where {M<:Measurement} =
    measurement(_meas_inner(M)(value(m)), _meas_inner(M)(uncertainty(m)))
_coerce_elt_to_T(m::Measurement, ::Type{R}) where {R<:AbstractFloat} = convert(R, value(m))
_coerce_elt_to_T(::Nothing, ::Type{T}) where {T} = nothing
_coerce_elt_to_T(::Missing, ::Type{T}) where {T} = missing
_coerce_elt_to_T(x::Bool, ::Type{M}) where {M<:Measurement} = x
_coerce_elt_to_T(x::Bool, ::Type{R}) where {R<:AbstractFloat} = x
_coerce_elt_to_T(x::Union{Symbol,String,Function,DataType}, ::Type{T}) where {T} = x
_coerce_elt_to_T(x, ::Type{T}) where {T} = x

# ---------- Public Coercion API (Containers & Scalars) ----------

# Promote Real to Complex when target is Complex.
coerce_to_T(x::Real, ::Type{C}) where {P,C<:Complex{P}} = C(coerce_to_T(x, P))

# General numbers fall back to element coercion.
coerce_to_T(x::Number, ::Type{T}) where {T} = _coerce_elt_to_T(x, T)

# Coerce a Complex number to a target Complex type.
coerce_to_T(x::Complex, ::Type{C}) where {P,C<:Complex{P}} =
    C(coerce_to_T(real(x), P), coerce_to_T(imag(x), P))

# Coerce a Complex number to a Real type (drops imaginary part).
coerce_to_T(x::Complex, ::Type{R}) where {R<:Real} = coerce_to_T(real(x), R)

# Handle containers recursively.
coerce_to_T(A::AbstractArray, ::Type{T}) where {T} = broadcast(y -> coerce_to_T(y, T), A)
coerce_to_T(t::Tuple, ::Type{T}) where {T} = map(y -> coerce_to_T(y, T), t)
coerce_to_T(nt::NamedTuple, ::Type{T}) where {T} =
    NamedTuple{keys(nt)}(map(v -> coerce_to_T(v, T), values(nt)))

# Fallback for other types.
coerce_to_T(x, ::Type{T}) where {T} = _coerce_elt_to_T(x, T)