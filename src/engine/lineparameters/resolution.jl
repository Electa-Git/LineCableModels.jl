"Revision of the shared declared line-observable resolution contract."
const OBSERVABLE_RESOLUTION_REVISION = 2

# Native basis units: ohm, H, S, F; per metre for :pul. These are reporting
# cutoffs, not a posteriori forward-error estimates from any solver.
const _LINE_RESOLUTION_DEFAULTS = (R=1e-10, L=1e-15, G=1e-12, C=1e-16)
const _LINE_RESOLUTION_QUANTITIES = (Z, Y, R, X, L, G, B, C)

function _validate_resolution_atol(atol)
    atol === nothing && return nothing
    limits = atol isa NamedTuple ? values(atol) : (atol,)
    if atol isa NamedTuple
        all(name -> name in (:Z, :Y, :R, :X, :L, :G, :B, :C), keys(atol)) ||
            throw(ArgumentError("atol keys must be Z, Y, R, X, L, G, B, or C"))
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

function _line_resolution_tolerance(quantity, ::Type{T}, f, atol) where {T}
    name = nameof(quantity)
    overrides = atol === nothing ? (;) : atol isa NamedTuple ? atol : NamedTuple{(name,)}((atol,))
    haskey(overrides, name) && return convert(T, getproperty(overrides, name))
    limits = merge(_LINE_RESOLUTION_DEFAULTS, overrides)
    name in keys(_LINE_RESOLUTION_DEFAULTS) && return convert(T, getproperty(limits, name))
    f === nothing && return nothing
    # A reporting cutoff uses nominal coordinates, just as its scalar precision
    # comes from nominal values. LEP may retain even exact frequencies in a
    # Measurement-typed vector. This does not alter that axis or its uncertainty.
    # Convert before arithmetic to preserve Float32/BigFloat precision.
    angular = T(2) * T(π) .* T.(nominal.(f))
    name === :X && return angular .* T(limits.L)
    name === :B && return angular .* T(limits.C)
    name === :Z && return T(limits.R) .+ angular .* T(limits.L)
    return T(limits.G) .+ angular .* T(limits.C)
end

_resolution_unresolved(value::Number, tolerance) =
    isfinite(value) && abs(nominal(value)) <= tolerance
_resolution_unresolved(value, tolerance) = false
_resolution_unresolved(value::Number, tolerance, ::Val{:uncertainty}) =
    isfinite(value) && isfinite(uncertainty(value)) && abs(uncertainty(value)) <= tolerance
function _resolution_unresolved(value::Complex, tolerance, ::Val{:uncertainty})
    spread = hypot(uncertainty(real(value)), uncertainty(imag(value)))
    return isfinite(value) && isfinite(spread) && spread <= tolerance
end
_resolution_unresolved(value, tolerance, ::Val{:uncertainty}) = false

"""
$(TYPEDSIGNATURES)

Resolve physical reporting cutoffs for detached numerical projections. The
quantity and native basis must be explicit; X/B/Z/Y also require frequency
context in Hz. Nominal magnitude and standard uncertainty are assessed
independently against the same cutoff. Complex spread uses the root-sum-square
of the real and imaginary standard uncertainties. This does not infer numerical
accuracy from tensor magnitudes.
"""
function observation_resolution(values::Union{Number,AbstractArray}, selector::Function;
        atol=nothing, frequencies=nothing, result_basis::Symbol=:pul)
    _validate_resolution_atol(atol)
    _check_basis(result_basis)
    selector in _LINE_RESOLUTION_QUANTITIES ||
        return observation_resolution(nothing, selector; atol, frequencies)
    if frequencies isa Real
        isfinite(frequencies) && frequencies>=0 ||
            throw(ArgumentError("resolution frequency must be finite and nonnegative in Hz"))
    elseif frequencies !== nothing
        frequencies isa AbstractVector &&
            all(value -> value isa Real && isfinite(value) && value >= 0, frequencies) ||
            throw(ArgumentError("resolution frequencies must be a finite nonnegative vector in Hz"))
        (values isa Number ? 1 : size(values, ndims(values))) == length(frequencies) ||
            throw(DimensionMismatch("resolution frequencies must match the tensor depth"))
    end
    T = typeof(float(real(nominal(zero(Base.nonmissingtype(eltype(values)))))))
    tolerance = _line_resolution_tolerance(selector, T, frequencies, atol)
    tolerance === nothing && return observation_resolution(nothing, selector; atol, frequencies)
    limits = tolerance isa Real ? (tolerance,) : tolerance
    all(value -> isfinite(value) && value >= 0, limits) || throw(ArgumentError(
        "reporting cutoffs must remain finite and nonnegative in the result precision"))
    aligned = tolerance isa AbstractArray ? values isa Number ? only(tolerance) :
        reshape(tolerance, (ntuple(_ -> 1, ndims(values)-1)..., length(tolerance))) : tolerance
    return (kind=:declared_floor, revision=OBSERVABLE_RESOLUTION_REVISION,
        atol=tolerance, unit=Units.native_unit(selector, result_basis),
        unresolved=_resolution_unresolved.(values, aligned),
        uncertainty_unresolved=_resolution_unresolved.(values, aligned, Val(:uncertainty)))
end

# These scientific arrays own their basis; route function requests through
# their native observation method rather than the detached-array contract.
observation_resolution(source::Union{SeriesImpedance,ShuntAdmittance}, request::Function;
        atol=nothing, frequencies=nothing) =
    observation_resolution(source, (request,); atol, frequencies)

"""
$(TYPEDSIGNATURES)

Apply the same physical cutoffs to a retained numerical publication. Frequency
coordinates are converted from their published unit before slicing. No source
result or absent statistic is reconstructed.
"""
function observation_resolution(source::ObservationPublication,request;atol=nothing,frequencies=nothing)
    resolved=observation_request(source,request)
    index=findfirst(selector -> Units.quantity(selector)==resolved.quantity,_LINE_RESOLUTION_QUANTITIES)
    index===nothing && return observation_resolution(nothing,request;atol,frequencies)
    values=request isa Tuple ? observe(source,request...) : observe(source,request)
    f=frequencies
    if f===nothing && haskey(source.columns,:frequency)
        f=unique(source.columns.frequency)
        contract=get(source.metadata.observation_columns,:frequency,nothing)
        contract===nothing || (f=f.*Units.scale_factor(contract.unit,Units.units(:base,:hertz)))
    end
    indices=resolved.indices
    physical=resolved.identity isa Tuple && length(resolved.identity)==3 && !isempty(indices) ? Base.tail(indices) : indices
    f===nothing || isempty(physical) || (f=f[last(physical)])
    return observation_resolution(values,_LINE_RESOLUTION_QUANTITIES[index];
        atol,frequencies=f,result_basis=basis(source))
end

"""
$(TYPEDSIGNATURES)

Resolve line-quantity reporting cutoffs in native basis units: R/X/Z in
\\[Ω/m\\], L in \\[H/m\\], G/B/Y in \\[S/m\\], C in \\[F/m\\] for `:pul`,
and the corresponding total units for `:total`. X/B follow `2πf` times the
L/C cutoffs; Z/Y use the sum of component cutoffs unless directly overridden.
Nominal magnitude and standard uncertainty have separate, aligned masks.
Resolved uncertainty is retained even when the nominal value is unresolved.
For phase requests, the masks assess the underlying complex quantity; its
nominal mask determines whether a phase can be published.
"""
function observation_resolution(source::Union{AbstractCoreResult, SeriesImpedance, ShuntAdmittance},
        request; atol=nothing, frequencies=nothing)
    _validate_resolution_atol(atol)
    identity = request_identity(request)
    selector = identity isa Tuple ? first(identity) : identity
    if !(selector in _LINE_RESOLUTION_QUANTITIES) ||
            (identity isa Tuple && !(last(identity) in (abs, angle, diag)))
        return observation_resolution(nothing, request; atol, frequencies)
    end
    f = _resolution_frequencies(source, frequencies)
    indices = request_indices(request)
    standalone_lc = source isa Union{SeriesImpedance, ShuntAdmittance} && selector in (L, C)
    if standalone_lc && !isempty(indices)
        embedded = first(indices)
        f === nothing || f == embedded || throw(ArgumentError("request and supplied frequencies differ"))
        f = embedded
        indices = Base.tail(indices)
    end
    if f !== nothing
        f isa AbstractVector && all(value -> value isa Real && isfinite(value) && value >= 0, f) ||
            throw(ArgumentError("resolution frequencies must be a finite nonnegative vector in Hz"))
        parent = selector in (Z, R, X, L) ? Z : Y
        size(observe(source, parent), 3) == length(f) ||
            throw(DimensionMismatch("resolution frequencies must match the tensor depth"))
    end
    phase = identity isa Tuple && last(identity) === angle
    values = phase ? observe(source, selector, indices...) :
        request isa Function ? observe(source, request) : observe(source, request...)
    diagonal = identity isa Tuple && last(identity) === diag
    sample_selector = length(indices) == (diagonal ? 2 : 3) ? last(indices) : Colon()
    selected_f = f === nothing ? nothing : f[sample_selector]
    return observation_resolution(values, selector; atol, frequencies=selected_f,
        result_basis=basis(source))
end
