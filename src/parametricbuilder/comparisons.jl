"""
$(TYPEDSIGNATURES)

Compare one reference with every candidate while retaining the problem and
formulation axes. The returned ParametricResult contains per-term RMSError
values in the same order as the candidate results. No formulation is evaluated.
"""
function compare(reference::Grammar.AbstractCoreResult, candidate::ParametricResult,
        quantity; kwargs...)
    isempty(candidate.axes) && throw(ArgumentError("comparison requires retained problem/formulation axes"))
    errors=[compare(reference, value, quantity; kwargs...) for value in candidate]
    return ParametricResult(candidate.formulation, errors, candidate.axes, (;))
end

"""
$(TYPEDSIGNATURES)

Compare two retained result spaces with an explicit reference index for every
candidate point. `pairing` is a vector of `(reference_index, candidate_index)`
pairs; each candidate must occur exactly once. Equal lengths do not imply pairing.
"""
function compare(reference::ParametricResult, candidate::ParametricResult,
        quantity; pairing=nothing, kwargs...)
    pairing === nothing && throw(ArgumentError("two result spaces require explicit pairing"))
    length(pairing) == length(candidate) &&
        sort(last.(pairing)) == collect(eachindex(candidate.values)) ||
        throw(ArgumentError("pairing must select every candidate exactly once"))
    all(pair -> pair isa Tuple{Integer,Integer} &&
        !(first(pair) isa Bool) && !(last(pair) isa Bool) &&
        1 <= first(pair) <= length(reference), pairing) ||
        throw(ArgumentError("pairing contains an invalid reference/candidate index"))
    ordered=sort(collect(pairing); by=last)
    errors=[compare(reference[i], candidate[j], quantity; kwargs...) for (i,j) in ordered]
    return ParametricResult(candidate.formulation, errors, candidate.axes, (;))
end
