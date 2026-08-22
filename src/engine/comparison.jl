"""
$(TYPEDEF)

Store element-wise absolute and reference-normalized root-mean-square errors.

Each matrix entry contains the error for the corresponding line-parameter
term evaluated across its complete frequency vector.

$(TYPEDFIELDS)
"""
struct RMSError{T <: Real}
    "Absolute error in the units of the compared quantity."
    absolute::Matrix{T}
    "Reference-normalized error [dimensionless]."
    relative::Matrix{T}
end

"""
$(TYPEDEF)

Store per-term frequency-domain errors for the series impedance and shunt
admittance of two [`LineParameters`](@ref) objects.

$(TYPEDFIELDS)
"""
struct LineParametersComparison{T <: Real}
    "Series-impedance error."
    Z::RMSError{T}
    "Shunt-admittance error."
    Y::RMSError{T}
end

function _rms_series(reference::AbstractVector, candidate::AbstractVector)
    difference_norm = sum(abs2, reference .- candidate)
    reference_norm = sum(abs2, reference)
    absolute = sqrt(difference_norm / length(reference))
    relative = if iszero(reference_norm)
        iszero(difference_norm) ? zero(absolute) : oftype(absolute, Inf)
    else
        sqrt(difference_norm / reference_norm)
    end
    return (; absolute, relative)
end

function _rms_error(reference::AbstractArray{<:Any, 3}, candidate::AbstractArray{<:Any, 3})
    errors = [_rms_series(@view(reference[row, column, :]), @view(candidate[row, column, :]))
              for row in axes(reference, 1), column in axes(reference, 2)]
    absolute = getproperty.(errors, :absolute)
    relative = getproperty.(errors, :relative)
    T = promote_type(eltype(absolute), eltype(relative))
    return RMSError{T}(Matrix{T}(absolute), Matrix{T}(relative))
end

"""
$(TYPEDSIGNATURES)

Compare two line-parameter results using absolute and reference-normalized
root-mean-square errors for each Z and Y matrix term across frequency.

For the reference series ``A_{ij}`` and candidate series ``B_{ij}`` at one
matrix term:

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
and domain. This operation does not reorder conductors, interpolate frequency
samples, convert basis, or apply a reduction.

When one reference term has zero norm, its relative error is zero for an
identical candidate term and `Inf` otherwise.

# Arguments

- `reference`: Reference line parameters.
- `candidate`: Candidate line parameters on the same frequency samples.

# Returns

- A [`LineParametersComparison`](@ref). Absolute Z errors use the impedance
  units selected by the result basis, and absolute Y errors use the corresponding
  admittance units. Relative errors are dimensionless.

# Errors

- `DimensionMismatch` when Z/Y tensor dimensions differ.
- `ArgumentError` when frequencies, basis, or domain differ.
"""
function compare(reference::LineParameters, candidate::LineParameters)
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
    impedance = _rms_error(reference.Z.values, candidate.Z.values)
    admittance = _rms_error(reference.Y.values, candidate.Y.values)
    T = promote_type(eltype(impedance.absolute), eltype(admittance.absolute))
    return LineParametersComparison{T}(
        RMSError{T}(Matrix{T}(impedance.absolute), Matrix{T}(impedance.relative)),
        RMSError{T}(Matrix{T}(admittance.absolute), Matrix{T}(admittance.relative))
    )
end
