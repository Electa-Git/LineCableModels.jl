"""
$(TYPEDEF)

Store element-wise absolute and reference-normalised root-mean-square benchmark errors.

Each matrix entry contains the error for the corresponding line-parameter
term over the selected frequency samples. Missing values represent explicit
non-applicability or an empty band, with the explanation retained in `details`.

$(TYPEDFIELDS)
"""
struct RMSError{T <: Real, D <: NamedTuple}
    "Absolute error in the units of the compared quantity."
    absolute::Matrix{Union{Missing, T}}
    "Reference-normalised error [dimensionless]."
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
        normalization::Symbol = :reference_rms)
    difference_norm = sum(abs2, reference .- candidate)
    reference_norm = sum(abs2, reference)
    absolute = sqrt(difference_norm / length(reference))
    relative = if normalization === :pointwise
        sqrt(sum(eachindex(reference, candidate)) do k
            difference = candidate[k] - reference[k]
            iszero(reference[k]) ?
            (iszero(difference) ? zero(absolute) : oftype(absolute, Inf)) :
            abs2(difference / reference[k])
        end / length(reference))
    elseif iszero(reference_norm)
        iszero(difference_norm) ? zero(absolute) : oftype(absolute, Inf)
    else
        sqrt(difference_norm / reference_norm)
    end
    return (; absolute, relative)
end

"""
$(TYPEDSIGNATURES)

Compare two line-parameter results using absolute and reference-normalised
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

When one reference term has zero norm, its relative error is zero for an
identical candidate term and `Inf` otherwise.

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
- `quantity`: `Z`, `Y`, `R`, `L`, `G`, or `C`. In particular, comparing `G`
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
- `atol`: Absolute numerical-zero tolerance, in the observable's basis units.
  A scalar applies to the requested quantity; a NamedTuple selects tolerances
  by quantity symbol. Defaults are 1e-10 for R, 1e-12 for G, 1e-15 for L,
  and 1e-16 for C, per metre for `:pul` and total units for `:total`.
  Unless overridden directly, Z uses `atol_R + 2πf*atol_L` and Y uses
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

If both complete traces lie within `atol`, their deviation is reported as zero
and classified `:below_tolerance`. Otherwise the ordinary RMS expression is
used without a denominator floor. Source arrays are never modified.

For `normalization=:pointwise`, the relative error is

```math
\\sqrt{\\frac{1}{N_f}\\sum_{k=1}^{N_f}
\\left|\\frac{B_{ij,k}-A_{ij,k}}{A_{ij,k}}\\right|^2}.
```

An exact zero-over-zero sample contributes zero; a nonzero difference over an
exact zero reference contributes `Inf`. Every selected sample remains in the
average. A finite reference-RMS error can therefore coexist with an infinite
pointwise error. The absolute RMS error is independent of normalization.

# Returns

- [`RMSError`](@ref), including actual bounds, sample indices/count, tolerance,
  reason, and a per-term status matrix in `details`.
"""
function compare(reference::LineParameters, candidate::LineParameters,
        quantity::Union{typeof(Z), typeof(Y), typeof(R), typeof(L), typeof(G), typeof(C)};
        normalization::Symbol = :reference_rms,
        band = :all, fundamental::Real = 50.0, harmonics::Integer = 50,
        atol = nothing, unsupported::NamedTuple = (;))
    normalization in (:reference_rms, :pointwise) || throw(ArgumentError(
        "normalization must be :reference_rms or :pointwise"))
    isempty(reference.f) && throw(ArgumentError("reference frequencies cannot be empty"))
    isempty(reference.Z.values) &&
        throw(ArgumentError("reference Z tensor cannot be empty"))
    isempty(reference.Y.values) &&
        throw(ArgumentError("reference Y tensor cannot be empty"))
    size(reference.Z) == size(candidate.Z) || throw(DimensionMismatch(
        "reference and candidate Z dimensions must match",
    ))
    size(reference.Y) == size(candidate.Y) || throw(DimensionMismatch(
        "reference and candidate Y dimensions must match",
    ))
    reference.f == candidate.f || throw(ArgumentError(
        "reference and candidate frequencies must match exactly and in order",
    ))
    basis(reference) === basis(candidate) || throw(ArgumentError(
        "reference and candidate basis must match",
    ))
    domain(reference) === domain(candidate) || throw(ArgumentError(
        "reference and candidate domains must match",
    ))
    isfinite(fundamental) && fundamental > 0 ||
        throw(ArgumentError("fundamental must be finite and positive Hz"))
    harmonics > 0 || throw(ArgumentError("harmonics must be a positive integer"))
    issorted(reference.f) ||
        throw(ArgumentError("frequency-band comparison requires ascending stored frequencies"))
    requested = if band === :all
        (first(reference.f), last(reference.f))
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
        (searchsortedlast(reference.f, 1e6) + 1):length(reference.f)
    elseif upper < first(reference.f) || lower > last(reference.f)
        1:0
    elseif band === :all
        1:length(reference.f)
    else
        first_index = argmin(abs.(reference.f .- lower))
        last_index = isinf(upper) ? length(reference.f) : argmin(abs.(reference.f .- upper))
        first_index:last_index
    end
    left = observe(reference, quantity)
    right = observe(candidate, quantity)
    T = promote_type(typeof(float(real(zero(eltype(left))))),
        typeof(float(real(zero(eltype(right))))))
    name = Symbol(nameof(quantity))
    defaults = (R = 1e-10, L = 1e-15, G = 1e-12, C = 1e-16)
    overrides = atol === nothing ? (;) :
                atol isa NamedTuple ? atol : NamedTuple{(name,)}((atol,))
    isempty(setdiff(keys(overrides), (:Z, :Y, :R, :L, :G, :C))) ||
        throw(ArgumentError("atol keys must be Z, Y, R, L, G, or C"))
    all(v -> v isa Real && isfinite(v) && v >= 0, values(overrides)) ||
        throw(ArgumentError("atol must be finite and nonnegative in $(basis(reference)) units"))
    limits = merge(defaults, overrides)
    tolerance = if haskey(overrides, name)
        fill(convert(T, getproperty(overrides, name)), length(indices))
    elseif name === :Z
        T[limits.R + 2π*f*limits.L for f in reference.f[indices]]
    elseif name === :Y
        T[limits.G + 2π*f*limits.C for f in reference.f[indices]]
    else
        fill(convert(T, getproperty(limits, name)), length(indices))
    end
    reason = get(unsupported, name, nothing)
    for result in (reference, candidate)
        declared = get(result.details, :comparison_unsupported, (;))
        reason === nothing && (reason = get(declared, name, nothing))
    end
    reason === nothing || reason isa AbstractString && !isempty(reason) ||
        throw(ArgumentError("unsupported comparisons require a nonempty explanatory string for $name"))
    status = reason !== nothing ? :unsupported : isempty(indices) ? :no_samples : :compared
    reason === nothing && isempty(indices) &&
        (reason = "No stored samples in the requested frequency band")
    absolute = Matrix{Union{Missing, T}}(missing, size(left, 1), size(left, 2))
    relative = similar(absolute)
    fill!(relative, missing)
    classifications = fill(status, size(absolute))
    if status === :compared
        for row in axes(left, 1), column in axes(left, 2)

            a = @view left[row, column, indices]
            b = @view right[row, column, indices]
            if all(abs.(a) .<= tolerance) && all(abs.(b) .<= tolerance)
                absolute[row, column] = relative[row, column] = zero(T)
                classifications[row, column] = :below_tolerance
            else
                error = _rms_series(a, b, normalization)
                absolute[row, column],
                relative[row, column] = error.absolute, error.relative
            end
        end
    end
    bounds = isempty(indices) ? (missing, missing) :
             (reference.f[first(indices)], reference.f[last(indices)])
    details = (; quantity = name, normalization, band,
        requested_bounds = requested, actual_bounds = bounds,
        indices, sample_count = length(indices), fundamental, harmonics, atol = tolerance,
        status = classifications, reason)
    return RMSError{T}(absolute, relative; details)
end
