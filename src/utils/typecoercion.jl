"""
$(TYPEDSIGNATURES)

Determines whether a `Type` contains or is a `Measurements.Measurement` somewhere in its structure. The check is recursive over arrays, tuples (including variadic tuples), named tuples, and union types. For concrete struct types, the predicate descends into field types. Guards are present to avoid infinite recursion through known self‑contained types (e.g., `Complex`).

# Arguments

- `::Type`: Type to inspect \\[dimensionless\\].

# Returns

- `Bool` indicating whether a `Measurement` occurs anywhere within the type structure.

# Notes

- For `AbstractArray{S}`, only the element type `S` is inspected.
- For `Tuple` and `NamedTuple{N,T}`, the parameters are traversed.
- For `Union`, both branches are inspected.
- Concrete `Complex` types are treated as terminal and are not descended.

# Examples

```julia
using Measurements

$(FUNCTIONNAME)(Float64)                    # false
$(FUNCTIONNAME)(Measurement{Float64})       # true
$(FUNCTIONNAME)(Vector{Measurement{Float64}}) # true
$(FUNCTIONNAME)(Tuple{Int, Float64})        # false
$(FUNCTIONNAME)(Union{Int, Measurement{Float64}}) # true
```
"""
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

"""
$(TYPEDSIGNATURES)

Determines whether a `Type` contains or is a `Complex` number type somewhere in its structure. The check is recursive over arrays, tuples (including variadic tuples), named tuples, and union types. For concrete struct types, the predicate descends into field types. Guards are present to avoid infinite recursion through known self‑referential types (e.g., `Measurements.Measurement`).

# Arguments

- `::Type`: Type to inspect \\[dimensionless\\].

# Returns

- `Bool` indicating whether a `Complex` type occurs anywhere within the type structure.

# Notes

- For `AbstractArray{S}`, only the element type `S` is inspected.
- For `Tuple` and `NamedTuple{N,T}`, the parameters are traversed.
- For `Union`, both branches are inspected.
- Concrete `Measurement` types are treated as terminal and are not descended.

# Examples

```julia
$(FUNCTIONNAME)(Float64)                  # false
$(FUNCTIONNAME)(Complex{Float64})         # true
$(FUNCTIONNAME)(Vector{ComplexF64})       # true
$(FUNCTIONNAME)(Tuple{Int, ComplexF64})   # true
```

# Methods

$(_CLEANMETHODLIST)
"""
function _hascomplex_type end
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

"""
$(TYPEDSIGNATURES)

Resolves the **promotion target type** to be used by constructors and coercion utilities based on the runtime arguments. The decision uses structure‑aware predicates for `Measurement` and `Complex`:

- If any argument contains `Measurement` and any contains `Complex`, returns `Complex{Measurement{BASE_FLOAT}}`.
- Else if any contains `Measurement`, returns `Measurement{BASE_FLOAT}`.
- Else if any contains `Complex`, returns `Complex{BASE_FLOAT}`.
- Otherwise returns `BASE_FLOAT`.

# Arguments

- `args...`: Values whose types will drive the promotion decision \\[dimensionless\\].

# Returns

- A `Type` suitable for numeric promotion in subsequent coercion.

# Examples

```julia
using Measurements

T = $(FUNCTIONNAME)(1.0, 2.0)                       # BASE_FLOAT
T = $(FUNCTIONNAME)(1 + 0im, 2.0)                    # Complex{BASE_FLOAT}
T = $(FUNCTIONNAME)(measurement(1.0, 0.1), 2.0)      # Measurement{BASE_FLOAT}
T = $(FUNCTIONNAME)(measurement(1.0, 0.1), 2 + 0im)  # Complex{Measurement{BASE_FLOAT}}
```
"""
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

"""
$(TYPEDSIGNATURES)

Extracts the real inner type `S` from `Measurement{S}`.

# Arguments

- `::Type{Measurement{S}}`: Measurement type wrapper \\[dimensionless\\].

# Returns

- The inner floating‐point type `S` \\[dimensionless\\].

# Examples

```julia
using Measurements

S = $(FUNCTIONNAME)(Measurement{Float64})  # Float64
```
"""
_meas_inner(::Type{Measurement{S}}) where {S} = S

"""
$(TYPEDSIGNATURES)

Element‑wise coercion kernel. Converts a *single leaf value* to the target type `T` while preserving semantics for `Measurement`, numeric types, and sentinels.

# Arguments

- `x`: Input leaf value \\[dimensionless\\].
- `::Type{T}`: Target type \\[dimensionless\\].

# Returns

- Value coerced to the target, according to the rules below.

# Notes

- `Number → R<:AbstractFloat`: uses `convert(R, x)`.
- `Number → M<:Measurement`: embeds the number as a zero‑uncertainty measurement (i.e., `zero(M) + x`).
- `Measurement → M<:Measurement`: recreates with the target inner type (value and uncertainty cast to `_meas_inner(M)`).
- `Measurement → R<:AbstractFloat`: drops uncertainty and converts the nominal value.
- `nothing` and `missing` pass through unchanged.
- `Bool`, `Symbol`, `String`, `Function`, `DataType`: passed through unchanged for measurement/real targets.
- Fallback: returns `x` unchanged.

# Examples

```julia
using Measurements

$(FUNCTIONNAME)(1.2, Float32)                         # 1.2f0
$(FUNCTIONNAME)(1.2, Measurement{Float64})            # 1.2 ± 0.0
$(FUNCTIONNAME)(measurement(2.0, 0.1), Float32)       # 2.0f0
$(FUNCTIONNAME)(measurement(2.0, 0.1), Measurement{Float32})  # 2.0 ± 0.1 (Float32 inner)
$(FUNCTIONNAME)(missing, Float64)                     # missing
```

# Methods

$(_CLEANMETHODLIST)
"""
function _coerce_elt_to_T end
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

"""
$(TYPEDSIGNATURES)

Public coercion API. Converts scalars and containers to a target type `T`, applying element‑wise coercion recursively. Complex numbers are handled by splitting into real and imaginary parts and coercing each side independently.

# Arguments

- `x`: Input value (scalar or container) \\[dimensionless\\].
- `::Type{T}`: Target type \\[dimensionless\\].

# Returns

- Value coerced to the target type:
  - For `Real → Complex{P}`: constructs `Complex{P}(coerce_to_T(re, P), coerce_to_T(im, P))` (imaginary part from `0`).
  - For `Complex → Real`: discards the imaginary part and coerces the real part.
  - For `AbstractArray`, `Tuple`, `NamedTuple`: coerces each element recursively.
  - For other types: defers to `_coerce_elt_to_T`.

# Examples

```julia
using Measurements

# Scalar
$(FUNCTIONNAME)(1.2, Float32)                         # 1.2f0
$(FUNCTIONNAME)(1.2, Measurement{Float64})            # 1.2 ± 0.0
$(FUNCTIONNAME)(1 + 2im, Complex{Float32})            # 1.0f0 + 2.0f0im
$(FUNCTIONNAME)(1 + 2im, Float64)                     # 1.0

# Containers
$(FUNCTIONNAME)([1.0, 2.0], Measurement{Float64})     # measurement array
$(FUNCTIONNAME)((1.0, 2.0), Float32)                  # (1.0f0, 2.0f0)
$(FUNCTIONNAME)((; a=1.0, b=2.0), Float32)            # (a = 1.0f0, b = 2.0f0)
```
# Methods

$(_CLEANMETHODLIST)

# See also

- [`_coerce_elt_to_T`](@ref)
- [`resolve_T`](@ref)
"""
function coerce_to_T end
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