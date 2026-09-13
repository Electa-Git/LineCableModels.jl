"""
$(TYPEDEF)

Store element-wise absolute and reference-normalized root-mean-square benchmark errors.

Each matrix entry contains the error for the corresponding line-parameter
term over the selected frequency samples. Missing values represent explicit
non-applicability, an empty band, or unavailable relative comparison, with
the explanation retained in `details`. Absolute differences remain measured
when only normalization is unavailable.

$(TYPEDFIELDS)
"""
struct RMSError{T <: Real, D <: NamedTuple}
    "Absolute error in the units of the compared quantity."
    absolute::Matrix{Union{Missing, T}}
    "Reference-normalized error [dimensionless]."
    relative::Matrix{Union{Missing, T}}
    "Requested and selected frequency band, tolerance, applicability, and per-term classification."
    details::D
end

function RMSError{T}(absolute::AbstractMatrix, relative::AbstractMatrix;
        details::NamedTuple = (;)) where {T <: Real}
    size(absolute) == size(relative) ||
        throw(DimensionMismatch("RMS error matrices must match"))
    return RMSError{T, typeof(details)}(absolute, relative, details)
end

function RMSError(absolute::AbstractMatrix{T}, relative::AbstractMatrix{S};
        details::NamedTuple = (;)) where {T <: Real, S <: Real}
    return RMSError{promote_type(T, S)}(absolute, relative; details)
end

"""Return RMS comparison metadata without exposing result storage to consumers."""
details(error::RMSError) = error.details
observe(error::RMSError, ::typeof(absolute_error)) = error.absolute
observe(error::RMSError, ::typeof(relative_error)) = error.relative
observe(error::RMSError, ::typeof(absolute_error), indices...) = getindex(error.absolute, indices...)
observe(error::RMSError, ::typeof(relative_error), indices...) = getindex(error.relative, indices...)

"""
$(TYPEDEF)

Store per-term frequency-domain errors for the series impedance and shunt
admittance of two [`LineParameters`](@ref) objects.

$(TYPEDFIELDS)
"""
struct LineParametersBenchmark{T <: Real, Basis, ZE <: RMSError{T}, YE <: RMSError{T}} <:
       AbstractProblemResult
    "Series-impedance error."
    Z::ZE
    "Shunt-admittance error."
    Y::YE
end

function LineParametersBenchmark(
        impedance::RMSError{T},
        admittance::RMSError{T};
        basis::Symbol = :pul
) where {T <: Real}
    _check_basis(basis)
    return LineParametersBenchmark{T, basis, typeof(impedance), typeof(admittance)}(impedance, admittance)
end

basis(::LineParametersBenchmark{T, Basis}) where {T, Basis} = Basis

function observe(benchmark::LineParametersBenchmark, ::typeof(Z), ::typeof(absolute_error), indices...)
    getindex(benchmark.Z.absolute, indices...)
end
function observe(benchmark::LineParametersBenchmark, ::typeof(Z), ::typeof(relative_error), indices...)
    getindex(benchmark.Z.relative, indices...)
end
function observe(benchmark::LineParametersBenchmark, ::typeof(Y), ::typeof(absolute_error), indices...)
    getindex(benchmark.Y.absolute, indices...)
end
function observe(benchmark::LineParametersBenchmark, ::typeof(Y), ::typeof(relative_error), indices...)
    getindex(benchmark.Y.relative, indices...)
end

function observe(benchmark::LineParametersBenchmark, ::typeof(Z), ::typeof(absolute_error))
    benchmark.Z.absolute
end
function observe(benchmark::LineParametersBenchmark, ::typeof(Z), ::typeof(relative_error))
    benchmark.Z.relative
end
function observe(benchmark::LineParametersBenchmark, ::typeof(Y), ::typeof(absolute_error))
    benchmark.Y.absolute
end
function observe(benchmark::LineParametersBenchmark, ::typeof(Y), ::typeof(relative_error))
    benchmark.Y.relative
end

function observables(::Type{<:LineParametersBenchmark})
    (
        (Z, absolute_error),
        (Z, relative_error),
        (Y, absolute_error),
        (Y, relative_error)
    )
end

function _rms_series(reference::AbstractVector, candidate::AbstractVector,
        normalization::Symbol, tolerance::AbstractVector, candidate_tolerance::AbstractVector)
    difference_norm = norm(reference .- candidate)
    sample_normalizer = sqrt(oftype(difference_norm, length(reference)))
    absolute = difference_norm / sample_normalizer
    unresolved_reference = count(abs.(reference) .<= tolerance)
    unresolved_candidate = count(abs.(candidate) .<= candidate_tolerance)
    counts = (; reference=unresolved_reference, candidate=unresolved_candidate)
    # Eligibility is two-sided and independent of normalization. Never reduce a
    # band's sample population to hide an undefined pairwise relative comparison.
    for (operand, count) in pairs(counts)
        if count > 0
            entire_trace = count == length(reference)
            status = if operand === :reference
                entire_trace ? :reference_below_tolerance : :reference_sample_below_tolerance
            else
                entire_trace ? :candidate_below_tolerance : :candidate_sample_below_tolerance
            end
            reason = "Samples at or below declared resolution: reference " *
                "$(counts.reference)/$(length(reference)), candidate $(counts.candidate)/$(length(reference)); " *
                "relative RMS requires both operands above tolerance at every selected sample; no samples were omitted"
            return (; absolute, relative = missing, status, reason, counts)
        end
    end
    relative = normalization === :pointwise ?
               norm((candidate .- reference) ./ reference) / sample_normalizer :
               difference_norm / norm(reference)
    return (; absolute, relative, status = :compared, reason = nothing, counts)
end

"""
    compare(reference::AbstractArray{<:Number,3}, candidate; normalization=:reference_rms, atol=0)

Measure per-entry absolute and relative RMS differences across the third axis.
The caller must establish equal physical coordinates, units and terminal order.
`atol` is a nonnegative scalar or one tolerance per sample, in the input units.
Relative RMS is `missing` if either operand's magnitude is at or below `atol`
at any selected sample, for either normalization. Absolute RMS retains every
measured difference, and `details` explains the unavailable relative comparison.
"""
function compare(
        reference::AbstractArray{<:Number, 3}, candidate::AbstractArray{<:Number, 3};
        normalization::Symbol = :reference_rms, atol = 0)
    axes(reference) == axes(candidate) ||
        throw(DimensionMismatch("RMS tensor axes must match"))
    isempty(reference) && throw(ArgumentError("RMS tensors cannot be empty"))
    normalization in (:reference_rms, :pointwise) ||
        throw(ArgumentError("unknown RMS normalization"))
    tolerance=atol isa Real ? fill(atol, size(reference, 3)) : collect(atol)
    length(tolerance) == size(reference, 3) ||
        throw(DimensionMismatch("one tolerance per sample is required"))
    all(value -> value isa Real && isfinite(value) && value >= 0, tolerance) ||
        throw(ArgumentError("RMS tolerances must be finite and nonnegative"))
    return _rms_arrays(reference, candidate, normalization, tolerance, tolerance)
end

function _rms_arrays(reference, candidate, normalization, tolerance, candidate_tolerance)
    all(isfinite, reference) && all(isfinite, candidate) ||
        throw(ArgumentError("RMS tensors must be finite"))
    T=promote_type(typeof(float(real(zero(eltype(reference))))), typeof(float(real(zero(eltype(candidate))))))
    errors=[_rms_series(view(reference, row, column, :),
                view(candidate, row, column, :), normalization, tolerance, candidate_tolerance)
            for row in axes(reference, 1), column in axes(reference, 2)]
    return RMSError{T}(getproperty.(errors, :absolute), getproperty.(errors, :relative);
        details = (normalization, atol = tolerance, sample_count = size(reference, 3),
            status = getproperty.(errors, :status), normalization_reason = getproperty.(errors, :reason),
            unresolved_samples = getproperty.(errors, :counts),
            resolution=(revision=OBSERVABLE_RESOLUTION_REVISION, kind=:explicit_floor, unit=nothing)))
end

"""
$(TYPEDSIGNATURES)

Compare two line-parameter results using absolute and reference-normalized
root-mean-square errors for each Z and Y matrix term across frequency.

For the reference series ``A_{ij}`` and candidate series ``B_{ij}`` at one
matrix term, the default `normalization=:reference_rms` uses:

```math
\\mathrm{RMS}_{\\mathrm{abs},ij} =
\\sqrt{\\frac{1}{N_f}\\sum_{k=1}^{N_f}\\left|A_{ij,k}-B_{ij,k}\\right|^2}
```

```math
\\mathrm{RMS}_{\\mathrm{rel},ij} =
\\sqrt{\\frac{\\sum_{k=1}^{N_f}\\left|A_{ij,k}-B_{ij,k}\\right|^2}
{\\sum_{k=1}^{N_f}\\left|A_{ij,k}\\right|^2}}
```

The operands must have identical frequency samples, tensor dimensions, basis,
and domain. Comparison does not reorder conductors, interpolate frequency
samples, convert basis, or apply a reduction.

Both operands must exceed the declared numerical-zero tolerance in magnitude
at every selected sample. Otherwise relative error is `missing`, including for
identical traces. Absolute RMS difference remains available.

Keyword arguments are shared with the single-observable `compare` method:
`normalization`, `band`, `fundamental`, `harmonics`, `atol`, and `unsupported`. Full-band error
is the default; optional sub-bands only slice the stored results. Each Z/Y
error retains its selected samples and numerical-zero classification.

# Arguments

- `reference`: Reference line parameters.
- `candidate`: Candidate line parameters on the same frequency samples.

# Returns

- A [`LineParametersBenchmark`](@ref). Absolute Z errors use the impedance
  units selected by the result basis, and absolute Y errors use the corresponding
  admittance units. Relative errors are dimensionless.

# Errors

- `DimensionMismatch` when Z/Y tensor dimensions differ.
- `ArgumentError` when frequencies, basis, or domain differ.
"""
function compare(reference::LineParameters, candidate::LineParameters; kwargs...)
    return LineParametersBenchmark(compare(reference, candidate, Z; kwargs...),
        compare(reference, candidate, Y; kwargs...); basis = basis(reference))
end

"""
$(TYPEDSIGNATURES)

Compare a selected line-parameter observable over stored frequency samples.
The first operand sets the relative-error normalization, not scientific truth.

# Arguments

- `reference`, `candidate`: Results with identical frequency coordinates, basis,
  domain, and matrix dimensions.
- `quantity`: `Z`, `Y`, `R`, `X`, `L`, `G`, `B`, or `C`. In particular, comparing `G`
  separately prevents displacement current from hiding dielectric-loss differences.

# Keywords

- `normalization`: `:reference_rms` (default) divides the absolute RMS error
  by the reference RMS. `:pointwise` instead computes the RMS of the
  sample-wise relative errors. Both return dimensionless fractions, not percentages.
- `band`: `:all` (default), explicit closed `(lower, upper)` Hz bounds, or
  `:dc` (0.1–100 Hz), `:harmonic` (`fundamental` to `harmonics*fundamental`),
  `:narrow` (1000–1000000 Hz), or `:wide` (strictly above 1000000 Hz).
- `fundamental`: Fundamental frequency in Hz; default 50.
- `harmonics`: Upper harmonic order; default 50. Every stored sample in the
  harmonic band is used, not only samples at integer harmonics.
- `atol`: Declared absolute reporting resolution in the observable's native basis
  units, not a certified floating-point error bound.
  A scalar applies to the requested quantity; a NamedTuple selects tolerances
  by quantity symbol. Defaults are 1e-10 for R, 1e-12 for G, 1e-15 for L,
  and 1e-16 for C, per meter for `:pul` and total units for `:total`.
  Unless overridden directly, X/B use `2πf*atol_L`/`2πf*atol_C`,
  Z uses `atol_R + 2πf*atol_L` and Y uses
  `atol_G + 2πf*atol_C` at each sample. A fixed admittance threshold would
  otherwise treat the same small capacitance differently across the spectrum.
- `unsupported`: NamedTuple of quantity symbols and explanatory strings for
  explicitly unavailable comparisons. Backend declarations may also supply
  this map as `details.comparison_unsupported` in either result.

Closed endpoints snap to the nearest stored frequency by absolute Hz distance;
ties select the lower frequency. Partial overlap uses the available portion.
Disjoint bands and an empty `:wide` band return `missing` errors with
`:no_samples`. No interpolation, extrapolation, weighting, or computation runs
are introduced. A one-sample band is valid.

Absolute RMS always retains the measured difference. Relative RMS requires both
operand magnitudes to exceed `atol` at every selected sample, for either
normalization. An entire trace within tolerance gives `missing` with status
`:reference_below_tolerance` or `:candidate_below_tolerance`; a partially
negligible trace gives `:reference_sample_below_tolerance` or
`:candidate_sample_below_tolerance`. The check is local to the selected band.
No denominator floor is introduced. Source arrays are never modified.

For `normalization=:pointwise`, the relative error is

```math
\\sqrt{\\frac{1}{N_f}\\sum_{k=1}^{N_f}
\\left|\\frac{B_{ij,k}-A_{ij,k}}{A_{ij,k}}\\right|^2}.
```

Samples are never omitted to obtain an eligible subset. Absolute RMS and relative
eligibility are independent of normalization. Each cell retains its explanation
in `details.normalization_reason`.

# Returns

- [`RMSError`](@ref), including actual bounds, sample indices/count, tolerance,
  reason, and a per-term status matrix in `details`.
"""
function compare(reference::AbstractCoreResult, candidate::AbstractCoreResult,
        quantity::Union{typeof(Z), typeof(Y), typeof(R), typeof(X), typeof(L), typeof(G), typeof(B), typeof(C)};
        normalization::Symbol = :reference_rms,
        band = :all, fundamental::Real = 50.0, harmonics::Integer = 50,
        atol = nothing, unsupported::NamedTuple = (;))
    validate(compare; normalization, band, fundamental, harmonics, atol, unsupported)
    left_coordinates=get(details(reference), :coordinates, nothing)
    right_coordinates=get(details(candidate), :coordinates, nothing)
    if left_coordinates !== nothing && right_coordinates !== nothing
        left_coordinates == right_coordinates || throw(ArgumentError("reference and candidate output terminal identities differ"))
    end
    f = frequencies(reference)
    isempty(f) && throw(ArgumentError("reference frequencies cannot be empty"))
    f == frequencies(candidate) || throw(ArgumentError(
        "reference and candidate frequencies must match exactly and in order"))
    basis(reference) === basis(candidate) || throw(ArgumentError("reference and candidate basis must match"))
    domain(reference) === domain(candidate) || throw(ArgumentError("reference and candidate domains must match"))
    issorted(f) || throw(ArgumentError("frequency-band comparison requires ascending stored frequencies"))
    left, right = observe(reference, quantity), observe(candidate, quantity)
    declared = merge(get(details(candidate), :comparison_unsupported, (;)),
        get(details(reference), :comparison_unsupported, (;)), unsupported)
    return compare(left, right, quantity; frequencies=f, result_basis=basis(reference),
        reference_resolution=observation_resolution(reference, quantity; atol, frequencies=f),
        candidate_resolution=observation_resolution(candidate, quantity; atol, frequencies=f),
        normalization, band, fundamental, harmonics, atol, unsupported=declared)
end

"""
$(TYPEDSIGNATURES)

Compare physical tensors on an explicit frequency axis using the same band and
two-sided resolution rules as core results. Arrays use native units in
`result_basis` (`:pul` or `:total`); frequencies use Hz. Owner-supplied resolution
records preserve operand-specific precision. No samples are omitted and no
physical calculation is performed.
"""
function compare(left::AbstractArray{<:Number,3}, right::AbstractArray{<:Number,3},
        quantity::Function; frequencies::AbstractVector, result_basis::Symbol,
        reference_resolution=nothing, candidate_resolution=nothing,
        normalization::Symbol=:reference_rms, band=:all, fundamental::Real=50.0,
        harmonics::Integer=50, atol=nothing, unsupported::NamedTuple=(;))
    validate(compare; normalization, band, fundamental, harmonics, atol, unsupported)
    _check_basis(result_basis)
    quantity in _LINE_RESOLUTION_QUANTITIES || throw(ArgumentError("unsupported physical RMS quantity"))
    f = frequencies
    !isempty(f) && issorted(f) && all(value -> isfinite(value) && value >= 0, f) ||
        throw(ArgumentError("frequency-band comparison requires finite ascending nonnegative frequencies"))
    size(left) == size(right) || throw(DimensionMismatch("reference and candidate quantity dimensions must match"))
    !isempty(left) && size(left, 3) == length(f) ||
        throw(DimensionMismatch("quantity dimensions must match the stored frequencies"))
    requested = if band === :all
        (first(f), last(f))
    elseif band === :dc
        (0.1, 100.0)
    elseif band === :harmonic
        (fundamental, fundamental * harmonics)
    elseif band === :narrow
        (1e3, 1e6)
    elseif band === :wide
        (1e6, Inf)
    elseif band isa Tuple{Real, Real}
        band
    else
        throw(ArgumentError("band must be :all, :dc, :harmonic, :narrow, :wide, or (lower, upper) in Hz"))
    end
    lower, upper = requested
    isfinite(lower) && lower >= 0 && !isnan(upper) && upper >= lower ||
        throw(ArgumentError("frequency bounds must satisfy 0 ≤ lower ≤ upper with finite lower"))
    indices = if band === :wide
        (searchsortedlast(f, 1e6) + 1):length(f)
    elseif upper < first(f) || lower > last(f)
        1:0
    elseif band === :all
        1:length(f)
    else
        first_index = argmin(abs.(f .- lower))
        last_index = isinf(upper) ? length(f) : argmin(abs.(f .- upper))
        first_index:last_index
    end
    T = promote_type(typeof(float(real(zero(eltype(left))))),
        typeof(float(real(zero(eltype(right))))))
    name = Symbol(nameof(quantity))
    resolution = reference_resolution === nothing ?
        observation_resolution(left, quantity; atol, frequencies=f, result_basis) : reference_resolution
    candidate_resolution = candidate_resolution === nothing ?
        observation_resolution(right, quantity; atol, frequencies=f, result_basis) : candidate_resolution
    tolerance = resolution.atol isa Real ? fill(T(resolution.atol), length(indices)) :
        T.(resolution.atol[indices])
    candidate_tolerance = candidate_resolution.atol isa Real ?
        fill(T(candidate_resolution.atol), length(indices)) : T.(candidate_resolution.atol[indices])
    reason = get(unsupported, name, nothing)
    reason === nothing || reason isa AbstractString && !isempty(reason) ||
        throw(ArgumentError("unsupported comparisons require a nonempty explanatory string for $name"))
    status = reason !== nothing ? :unsupported : isempty(indices) ? :no_samples : :compared
    reason === nothing && isempty(indices) &&
        (reason = "No stored samples in the requested frequency band")
    absolute = Matrix{Union{Missing, T}}(missing, size(left, 1), size(left, 2))
    relative = similar(absolute)
    fill!(relative, missing)
    classifications = fill(status, size(absolute))
    normalization_reasons = Matrix{Union{Nothing, String}}(nothing, size(absolute))
    unresolved_samples = fill((reference=0, candidate=0), size(absolute))
    if status === :compared
        error=_rms_arrays(left[:, :, indices], right[:, :, indices], normalization,
            tolerance, candidate_tolerance)
        absolute .= error.absolute
        relative .= error.relative
        classifications .= error.details.status
        normalization_reasons .= error.details.normalization_reason
        unresolved_samples .= error.details.unresolved_samples
    end
    bounds = isempty(indices) ? (missing, missing) :
             (f[first(indices)], f[last(indices)])
    comparison_details = (; quantity = name, normalization, band,
        requested_bounds = requested, actual_bounds = bounds,
        indices, sample_count = length(indices), fundamental, harmonics, atol = tolerance,
        candidate_atol = candidate_tolerance,
        status = classifications, reason, normalization_reason = normalization_reasons,
        unresolved_samples,
        resolution=(; resolution.revision, resolution.kind, resolution.unit))
    # Empty bands and supported bands have the same result type on a Gridspace.
    # Preserve the frequency scalar type while admitting an absent bound/reason.
    detail_types=map(keys(comparison_details)) do key
        key === :actual_bounds ? NTuple{2,Union{Missing,eltype(f)}} :
        key === :reason ? Union{Nothing,String} : typeof(getproperty(comparison_details,key))
    end
    stable_details=NamedTuple{keys(comparison_details),Tuple{detail_types...}}(values(comparison_details))
    return RMSError{T}(absolute, relative; details=stable_details)
end

"""Compare explicitly contextualized detached products without inferring coordinates."""
function compare(reference::NamedTuple{(:result,:metadata)},candidate::NamedTuple{(:result,:metadata)},
        request; kwargs...)
    a,b=reference.metadata,candidate.metadata
    a.port_order==b.port_order && a.frequencies==b.frequencies && a.basis==b.basis && a.domain==b.domain ||
        throw(ArgumentError("detached comparison coordinates, basis or domain differ"))
    identity=request_identity(request)
    identity isa Function || identity isa Tuple && length(identity)==3 ||
        throw(ArgumentError("detached comparisons require a physical quantity or a selected statistical product"))
    prefix=identity isa Tuple ? identity : (identity,)
    left=if reference.result isa AbstractCoreResult || reference.result isa ObservationPublication
        observe(reference.result,prefix...)
    else
        length(reference.result)==1 || throw(ArgumentError("select one UQ point before comparing with a detached publication"))
        observe(reference.result,prefix...,1)
    end
    right=if candidate.result isa AbstractCoreResult || candidate.result isa ObservationPublication
        observe(candidate.result,prefix...)
    else
        length(candidate.result)==1 || throw(ArgumentError("select one UQ point before comparing with a detached publication"))
        observe(candidate.result,prefix...,1)
    end
    quantity=identity isa Tuple ? identity[2] : identity
    error=compare(left,right,quantity;frequencies=a.frequencies,result_basis=a.basis,kwargs...)
    semantics=identity isa Tuple ? (statistical_semantics=(revision=1,representation=:retained_products),) : (;)
    return RMSError{Base.nonmissingtype(eltype(observe(error,absolute_error)))}(
        observe(error,absolute_error),observe(error,relative_error);
        details=merge(details(error),(;request=identity),semantics))
end

"""
$(TYPEDSIGNATURES)

Validate RMS comparison controls before accessing results or starting a calculation.
Frequency bounds and `fundamental` use Hz. Absolute tolerances use the units of
the selected quantities. This method performs no numerical calculation.
"""
function validate(::typeof(compare); normalization = :reference_rms, band = :all,
        fundamental = 50.0, harmonics = 50, atol = nothing, unsupported = (;))
    normalization in (:reference_rms, :pointwise) || throw(ArgumentError(
        "normalization must be :reference_rms or :pointwise"))
    fundamental isa Real && isfinite(fundamental) && fundamental > 0 ||
        throw(ArgumentError("fundamental must be finite and positive Hz"))
    harmonics isa Integer && !(harmonics isa Bool) && harmonics > 0 ||
        throw(ArgumentError("harmonics must be a positive integer"))
    if band isa Tuple{Real, Real}
        lower, upper = band
        isfinite(lower) && lower >= 0 && !isnan(upper) && upper >= lower ||
            throw(ArgumentError("frequency bounds must satisfy 0 ≤ lower ≤ upper with finite lower"))
    else
        band in (:all, :dc, :harmonic, :narrow, :wide) || throw(ArgumentError(
            "band must be :all, :dc, :harmonic, :narrow, :wide, or (lower, upper) in Hz"))
    end
    _validate_resolution_atol(atol)
    unsupported isa NamedTuple && isempty(setdiff(keys(unsupported), (:Z, :Y, :R, :X, :L, :G, :B, :C))) ||
        throw(ArgumentError("unsupported must name Z, Y, R, X, L, G, B, or C"))
    all(reason -> reason isa AbstractString && !isempty(reason), unsupported) ||
        throw(ArgumentError("unsupported comparisons require nonempty explanatory strings"))
    return nothing
end
