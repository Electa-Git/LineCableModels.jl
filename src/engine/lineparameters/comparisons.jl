"""
$(TYPEDSIGNATURES)

Compare one reference with every result while retaining the problem and
formulation axes. The returned ParametricResult contains per-term RMSError
values in the same order as the input results. No formulation is evaluated.
"""
function compare(reference::AbstractCoreResult, result::ParametricResult,
        quantity::Union{Function,Tuple}; kwargs...)
    isempty(result.axes) && throw(ArgumentError("comparison requires retained problem/formulation axes"))
    errors=[compare(reference, value, quantity::Union{Function,Tuple}; kwargs...) for value in result]
    return ParametricResult(result.formulation, errors, result.axes, ComputationDetails())
end

"""
$(TYPEDSIGNATURES)

Compare two retained result spaces with an explicit reference index for every
result point. `pairing` is a vector of `(reference_index, result_index)`
pairs. Each result must occur exactly once. Use the declared pairs even when the lengths are equal.
"""
function compare(reference::ParametricResult, result::ParametricResult,
        quantity::Union{Function,Tuple}; pairing=nothing, kwargs...)
    pairing === nothing && throw(ArgumentError("two result spaces require explicit pairing"))
    length(pairing) == length(result) &&
        sort(last.(pairing)) == collect(eachindex(result.values)) ||
        throw(ArgumentError("pairing must select every result exactly once"))
    foreach(pair -> validate(pair, compare, reference), pairing)
    ordered=sort(collect(pairing); by=last)
    errors=[compare(reference[i], result[j], quantity::Union{Function,Tuple}; kwargs...) for (i,j) in ordered]
    return ParametricResult(result.formulation, errors, result.axes, ComputationDetails())
end

# One pair of a comparison pairing: a reference index within the reference result, and a
# result index. Bools are not indices.
function validate(
        pair::Tuple{Integer, Integer},
        ::typeof(compare),
        reference::ParametricResult
)
    1 <= first(pair) <= length(reference) || throw(ArgumentError(
        "pairing contains an invalid reference/result index; got $pair"))
    return pair
end
function validate(
        pair::Union{Tuple{Bool, Integer}, Tuple{Integer, Bool}},
        ::typeof(compare),
        reference::ParametricResult
)
    throw(ArgumentError("pairing contains an invalid reference/result index; got $pair"))
end
function validate(pair, ::typeof(compare), reference::ParametricResult)
    throw(ArgumentError(
        "pairing contains an invalid reference/result index; got $(repr(pair))"))
end
