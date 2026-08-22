"""
$(TYPEDEF)

Store the absolute and reference-normalized root-mean-square errors for one
line-parameter tensor.

$(TYPEDFIELDS)
"""
struct RMSError{T <: Real}
    "Absolute root-mean-square error."
    absolute::T
    "Root-mean-square error normalized by the reference norm."
    relative::T
end

"""
$(TYPEDEF)

Store root-mean-square errors for the series-impedance and shunt-admittance
tensors of two [`LineParameters`](@ref) objects.

$(TYPEDFIELDS)
"""
struct LineParametersComparison{T <: Real}
    "Series-impedance error."
    Z::RMSError{T}
    "Shunt-admittance error."
    Y::RMSError{T}
end

function _rms_error(reference::AbstractArray, candidate::AbstractArray)
    difference_norm = sum(abs2, reference .- candidate)
    reference_norm = sum(abs2, reference)
    absolute = sqrt(difference_norm / length(reference))
    relative = if iszero(reference_norm)
        iszero(difference_norm) ? zero(absolute) : oftype(absolute, Inf)
    else
        sqrt(difference_norm / reference_norm)
    end
    T = promote_type(typeof(absolute), typeof(relative))
    return RMSError{T}(convert(T, absolute), convert(T, relative))
end

"""
$(TYPEDSIGNATURES)

Compare two line-parameter results using absolute and reference-normalized
root-mean-square errors for their complete Z and Y tensors.

The operands must have identical frequency samples, tensor dimensions, basis,
and domain. This operation does not reorder conductors, interpolate frequency
samples, convert basis, or apply a reduction.

When a reference tensor has zero norm, the relative error is zero for an
identical candidate and `Inf` otherwise.

# Returns

A [`LineParametersComparison`](@ref).

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
    T = promote_type(typeof(impedance.absolute), typeof(admittance.absolute))
    return LineParametersComparison{T}(
        RMSError{T}(convert(T, impedance.absolute), convert(T, impedance.relative)),
        RMSError{T}(convert(T, admittance.absolute), convert(T, admittance.relative))
    )
end
