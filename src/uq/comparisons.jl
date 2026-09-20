"""
$(TYPEDSIGNATURES)

Compare the selected statistic of two UQ result spaces without pooling points
or trials. `pairing` lists one reference/candidate pair per candidate, in any
order. A singleton reference can be broadcast; multiple references require
explicit pairing. Return RMS results in candidate order. Physical cutoffs,
frequency bands and normalization are owned by Engine.
"""
function Engine.compare(reference::Union{MonteCarloResult, LinearErrorResult},
        candidate::Union{MonteCarloResult, LinearErrorResult}, request::Tuple;
        pairing = nothing, kwargs...)
    identity = request_identity(request)
    identity isa Tuple && length(identity) == 3 && first(identity) === statistics ||
        throw(ArgumentError("UQ comparison requires (statistics, quantity, statistic)"))
    isempty(request_indices(request)) || throw(ArgumentError(
        "comparison point selection uses pairing; tensor entries are compared individually"))
    observation_request(reference, request)
    observation_request(candidate, request)
    selected = pairing === nothing ?
               (length(reference) == 1 ? [(1, index) for index in eachindex(candidate)] :
                throw(ArgumentError("multiple reference points require explicit pairing"))) :
               collect(pairing)
    length(selected) == length(candidate) &&
    all(
        pair -> pair isa Tuple{Integer, Integer} && !(first(pair) isa Bool) &&
                !(last(pair) isa Bool) && 1 <= first(pair) <= length(reference),
        selected) &&
    sort(last.(selected)) == collect(eachindex(candidate)) ||
        throw(ArgumentError("pairing must reference valid points and select each candidate exactly once"))
    sort!(selected; by = last)
    for (left_index, right_index) in selected
        left, right = reference[left_index], candidate[right_index]
        basis(left) === basis(right) && Engine.domain(left) === Engine.domain(right) ||
            throw(ArgumentError("UQ point basis or domain differs"))
        frequencies(left) == frequencies(right) ||
            throw(ArgumentError("UQ point frequencies differ"))
        left_ports = get(details(left).data, :coordinates, nothing)
        right_ports = get(details(right).data, :coordinates, nothing)
        left_ports === nothing || right_ports === nothing || left_ports == right_ports ||
            throw(ArgumentError("UQ point output terminals differ"))
    end
    return map(selected) do (left_index, right_index)
        left=reference[left_index]
        a = observe(reference, identity..., left_index)
        b = observe(candidate, identity..., right_index)
        error = Engine.compare(a, b, identity[2]; frequencies = frequencies(left),
            result_basis = basis(left), kwargs...)
        return Engine.RMSError{Base.nonmissingtype(eltype(observe(error, Engine.absolute_error)))}(
            observe(error, Engine.absolute_error), observe(error, Engine.relative_error);
            details = ComputationDetails(merge(details(error).data, (; request = identity,
                estimators=(
                    reference=reference isa MonteCarloResult ? :empirical : :first_order,
                    candidate=candidate isa MonteCarloResult ? :empirical : :first_order)))))
    end
end
