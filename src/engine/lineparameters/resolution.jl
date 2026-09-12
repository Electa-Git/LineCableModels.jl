"Revision of the shared declared line-observable resolution contract."
const OBSERVABLE_RESOLUTION_REVISION = 1

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
    # Convert before arithmetic so Float32 and BigFloat do not silently inherit
    # Float64 arithmetic, and presentation units cannot affect classification.
    angular = T(2) * T(π) .* T.(f)
    name === :X && return angular .* T(limits.L)
    name === :B && return angular .* T(limits.C)
    name === :Z && return T(limits.R) .+ angular .* T(limits.L)
    return T(limits.G) .+ angular .* T(limits.C)
end

_resolution_uncertain(value::Number) = !iszero(uncertainty(value))
_resolution_uncertain(value::Complex) =
    !iszero(uncertainty(real(value))) || !iszero(uncertainty(imag(value)))
_resolution_unresolved(value::Number, tolerance) =
    isfinite(value) && !_resolution_uncertain(value) && abs(nominal(value)) <= tolerance
_resolution_unresolved(value, tolerance) = false

"""
$(TYPEDSIGNATURES)

Resolve line-quantity reporting cutoffs in native basis units: R/X/Z in
\\[Ω/m\\], L in \\[H/m\\], G/B/Y in \\[S/m\\], C in \\[F/m\\] for `:pul`,
and the corresponding total units for `:total`. X/B follow `2πf` times the
L/C cutoffs; Z/Y use the sum of component cutoffs unless directly overridden.
Physical uncertainty is retained, not classified as deterministic residue.
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
    T = typeof(float(real(nominal(zero(eltype(values))))))
    diagonal = identity isa Tuple && last(identity) === diag
    sample_selector = length(indices) == (diagonal ? 2 : 3) ? last(indices) : Colon()
    selected_f = f === nothing ? nothing : f[sample_selector]
    tolerance = _line_resolution_tolerance(selector, T, selected_f, atol)
    tolerance === nothing && return observation_resolution(nothing, request; atol, frequencies)
    limits = tolerance isa Real ? (tolerance,) : tolerance
    all(limit -> isfinite(limit) && limit >= 0, limits) || throw(ArgumentError(
        "declared resolution must remain finite and nonnegative in the result precision"))
    aligned = tolerance isa AbstractArray ?
        reshape(tolerance, (ntuple(_ -> 1, ndims(values) - 1)..., length(tolerance))) : tolerance
    unresolved = _resolution_unresolved.(values, aligned)
    return (kind=:declared_floor, revision=OBSERVABLE_RESOLUTION_REVISION,
        atol=tolerance, unit=Units.native_unit(selector, basis(source)), unresolved)
end
