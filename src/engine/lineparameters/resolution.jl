# Exact decimal engineering floors in native per-metre units. Convert before
# frequency or length scaling; these are not solver-error estimates.
const _LINE_RESOLUTION_DEFAULTS = (R=1//10^10, L=1//10^15, G=1//10^12, C=1//10^16)
const _LINE_RESOLUTION_QUANTITIES = (Z, Y, R, X, L, G, B, C)

function _validate_resolution_atol(atol)
    atol === nothing && return nothing
    limits = atol isa NamedTuple ? values(atol) : (atol,)
    if atol isa NamedTuple
        all(name -> name in (:R,:X,:L,:G,:B,:C), keys(atol)) ||
            throw(ArgumentError("atol keys must be R, X, L, G, B, or C; complex zero uses component cutoffs"))
        for pair in ((:X,:L),(:B,:C))
            all(key -> haskey(atol,key),pair) && throw(ArgumentError(
                "supply only one cutoff for $(join(pair, '/')); their thresholds are linked by 2πf"))
        end
    end
    all(value -> value isa Real && !(value isa Bool) && isfinite(value) && value >= 0, limits) ||
        throw(ArgumentError("atol must be finite and nonnegative"))
    return nothing
end

_resolution_frequencies(source, supplied) = supplied
function _resolution_frequencies(source::LineParameters, supplied)
    supplied === nothing || supplied == source.f || throw(ArgumentError(
        "supplied frequencies must match the stored LineParameters frequency axis"))
    return source.f
end

_cutoff_number(::Type{T}, value::Rational) where {T} = T(value)
_cutoff_number(::Type{T}, value::Real) where {T} = convert(promote_type(T,typeof(float(value))),value)

function _line_resolution_tolerance(quantity, ::Type{T}, f, atol;
        result_basis=:pul, line_length=nothing) where {T}
    name = nameof(quantity)
    atol isa Real && name in (:Z,:Y) && throw(ArgumentError(
        "complex zero requires component cutoffs, for example atol=(R=0, X=0)"))
    overrides = atol === nothing ? (;) : atol isa NamedTuple ? atol : NamedTuple{(name,)}((atol,))
    component = function (selector)
        key=nameof(selector)
        haskey(overrides,key) && return _cutoff_number(T,overrides[key])
        partner = key === :X ? :L : key === :L ? :X : key === :B ? :C : key === :C ? :B : nothing
        linked = key in (:X,:B) || (partner !== nothing && haskey(overrides,partner))
        if linked
            f === nothing && return nothing
            angular = T(2)*T(π) .* T.(nominal.(f))
            other = _line_resolution_tolerance(getfield(@__MODULE__,partner),T,f,overrides;
                result_basis,line_length)
            other === nothing && return nothing
            # L/C at zero frequency are unavailable; never divide by zero.
            return key in (:X,:B) ? angular .* other :
                broadcast((cutoff,w) -> iszero(w) ? zero(cutoff) : cutoff/w,other,angular)
        end
        cutoff=T(getproperty(_LINE_RESOLUTION_DEFAULTS,key))
        result_basis === :pul && return cutoff
        line_length === nothing && throw(ArgumentError(
            "total $key reporting requires retained line length or an explicit total-unit cutoff"))
        length_value=nominal(line_length)
        isfinite(length_value) && length_value>0 || throw(ArgumentError("line length must be finite and positive"))
        return cutoff * length_value
    end
    name === :Z && return (R=component(R),X=component(X))
    name === :Y && return (G=component(G),B=component(B))
    return component(quantity)
end

_resolution_available(value) = false
_resolution_available(value::Real) = isfinite(nominal(value)) && isfinite(uncertainty(value))
_resolution_available(value::Complex) = _resolution_available(real(value)) && _resolution_available(imag(value))
_resolution_unresolved(value::Number, tolerance::Real) =
    _resolution_available(value) && abs(nominal(value)) <= tolerance
_resolution_unresolved(value, tolerance) = false

_aligned_cutoff(tolerance::Real, values) = tolerance
_aligned_cutoff(tolerance::AbstractArray, values::Number) = only(tolerance)
_aligned_cutoff(tolerance::AbstractArray, values::AbstractArray) =
    reshape(tolerance, (ntuple(_ -> 1,ndims(values)-1)...,length(tolerance)))
_cutoff_available(value::Real) = isfinite(value) && value>=0
_cutoff_available(value::AbstractArray) = all(_cutoff_available,value)
_cutoff_available(value::NamedTuple) = all(_cutoff_available,values(value))
_cutoff_available(::Nothing) = false

"""
$(TYPEDSIGNATURES)

Classify original quantities using `abs(nominal(value)) ≤ cutoff`. Complex zero
requires both Cartesian components to satisfy their cutoffs. Availability is
separate; physical standard uncertainty never contributes to a reporting floor.

# Keywords

- `atol`: Component cutoffs in native basis units.
- `frequencies`: Aligned frequencies \\[Hz\\] for linked X/L and B/C cutoffs.
- `result_basis=:pul`: Per-metre or `:total` quantities.
- `line_length`: Physical length \\[m\\] for scaling default total-unit floors.

# Returns

- Applied thresholds, native units, and aligned availability and nominal-zero
  masks. An unknown frequency-dependent threshold is explicitly unassessed.
"""
function observation_resolution(values::Union{Number,AbstractArray}, selector::Function;
        atol=nothing, frequencies=nothing, result_basis::Symbol=:pul, line_length=nothing)
    _validate_resolution_atol(atol)
    _check_basis(result_basis)
    selector in _LINE_RESOLUTION_QUANTITIES ||
        return observation_resolution(nothing,selector;atol,frequencies)
    if frequencies isa Real
        isfinite(frequencies) && frequencies>=0 || throw(ArgumentError("frequency must be finite and nonnegative in Hz"))
    elseif frequencies !== nothing
        frequencies isa AbstractVector && all(f -> f isa Real && isfinite(f) && f>=0,frequencies) ||
            throw(ArgumentError("frequencies must be finite and nonnegative in Hz"))
        (values isa Number ? 1 : size(values,ndims(values))) == length(frequencies) ||
            throw(DimensionMismatch("frequencies must match the quantity depth"))
    end
    T=typeof(float(real(nominal(zero(Base.nonmissingtype(eltype(values)))))))
    tolerance=_line_resolution_tolerance(selector,T,frequencies,atol;result_basis,line_length)
    assessed = tolerance !== nothing && (!(tolerance isa NamedTuple) || all(!isnothing,Base.values(tolerance)))
    available=_resolution_available.(values)
    if selector in (L,C) && frequencies!==nothing
        available = available .& (.!iszero.(_aligned_cutoff(frequencies,values)))
    end
    assessed || return (kind=:unassessed,atol=nothing,unit=Units.native_unit(selector,result_basis),
        unresolved=nothing,available)
    _cutoff_available(tolerance) || throw(ArgumentError("cutoffs must remain finite and nonnegative"))
    unresolved = if selector in (Z,Y)
        first_limit,last_limit=Base.values(tolerance)
        _resolution_unresolved.(real.(values),_aligned_cutoff(first_limit,values)) .&
            _resolution_unresolved.(imag.(values),_aligned_cutoff(last_limit,values))
    else
        _resolution_unresolved.(values,_aligned_cutoff(tolerance,values))
    end
    return (kind=:declared_floor,atol=tolerance,unit=Units.native_unit(selector,result_basis),
        unresolved=unresolved .& available,available)
end

observation_resolution(source::Union{SeriesImpedance,ShuntAdmittance},request::Function;
    atol=nothing,frequencies=nothing) = observation_resolution(source,(request,);atol,frequencies)

# Observation and comparison admit unavailable DC proxies without changing the
# strict raw L/C accessors used by numerical algorithms.
function _line_observation_values(source, request; frequencies=nothing)
    identity=request_identity(request)
    selector=identity isa Tuple ? first(identity) : identity
    indices=request_indices(request)
    if selector in (L,C) && source isa Union{LineParameters,SeriesImpedance,ShuntAdmittance}
        f=_resolution_frequencies(source,frequencies)
        if source isa Union{SeriesImpedance,ShuntAdmittance} && !isempty(indices) && length(indices)>3 && first(indices) isa AbstractVector{<:Real}
            f=first(indices)
            indices=Base.tail(indices)
        end
        f===nothing && throw(ArgumentError("L/C observations require frequencies in Hz"))
        transforms=identity isa Tuple ? Base.tail(identity) : ()
        diagonal=diag in transforms
        proxy=selector===L ? X : B
        values=diagonal ? observe(source,proxy,diag,indices...) : observe(source,proxy,indices...)
        sample=length(indices)==(diagonal ? 2 : 3) ? last(indices) : Colon()
        selected=f[sample]
        aligned=_aligned_cutoff(selected,values)
        T=typeof(float(real(nominal(zero(eltype(values))))))
        V=typeof(zero(eltype(values))/(T(2)*T(π)*one(eltype(f))))
        output=Array{Union{Missing,V}}(undef,values isa Number ? () : size(values))
        broadcast!(output,values,aligned) do value,frequency
            iszero(frequency) && return missing
            value/(T(2)*T(π)*frequency)
        end
        return output
    end
    return request isa Function ? observe(source,request) : observe(source,request...)
end

"""
$(TYPEDSIGNATURES)

Classify an owned line observation in native units, using original complex
components for polar requests. This operation never reads a prepared table.
"""
function observation_resolution(source::Union{AbstractCoreResult,SeriesImpedance,ShuntAdmittance},
        request;atol=nothing,frequencies=nothing)
    _validate_resolution_atol(atol)
    identity=request_identity(request)
    selector=identity isa Tuple ? first(identity) : identity
    selector in _LINE_RESOLUTION_QUANTITIES || return observation_resolution(nothing,request;atol,frequencies)
    f=_resolution_frequencies(source,frequencies)
    indices=request_indices(request)
    if source isa Union{SeriesImpedance,ShuntAdmittance} && selector in (L,C) && !isempty(indices)
        embedded=first(indices)
        f===nothing || f==embedded || throw(ArgumentError("request and supplied frequencies differ"))
        f=embedded
        indices=Base.tail(indices)
    end
    if f!==nothing
        f isa AbstractVector && all(x -> x isa Real && isfinite(x) && x>=0,f) ||
            throw(ArgumentError("frequencies must be finite and nonnegative in Hz"))
        parent=selector in (Z,R,X,L) ? Z : Y
        size(observe(source,parent),3)==length(f) || throw(DimensionMismatch("frequencies must match tensor depth"))
    end
    transforms=identity isa Tuple ? Base.tail(identity) : ()
    polar=any(transform -> transform in (abs,angle),transforms)
    diagonal=diag in transforms
    values = if polar
        diagonal ? observe(source,selector,diag,indices...) : observe(source,selector,indices...)
    else
        _line_observation_values(source,request;frequencies=f)
    end
    sample=length(indices)==(diagonal ? 2 : 3) ? last(indices) : Colon()
    return observation_resolution(values,selector;atol,frequencies=f===nothing ? nothing : f[sample],
        result_basis=basis(source),line_length=_resolution_length(source))
end
