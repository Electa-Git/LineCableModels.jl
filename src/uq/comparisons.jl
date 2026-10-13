"""
$(TYPEDSIGNATURES)

Compare the selected statistic of two UQ result spaces without pooling points
or trials. `pairing` lists one reference-result pair per result, in any
order. A singleton reference can be broadcast. Multiple references require
explicit pairing. Return RMS results in result order. Physical cutoffs,
frequency bands and normalization are owned by Engine.
"""
function Engine.compare(reference::Union{MonteCarloResult, LinearErrorResult},
        result::Union{MonteCarloResult, LinearErrorResult}, request::Tuple;
        pairing = nothing, kwargs...)
    identity = request_identity(request)
    identity isa Tuple && length(identity) == 3 && first(identity) === statistics ||
        throw(ArgumentError("UQ comparison requires (statistics, quantity, statistic)"))
    isempty(request_indices(request)) || throw(ArgumentError(
        "comparison point selection uses pairing; tensor entries are compared individually"))
    observation_request(reference, request)
    observation_request(result, request)
    selected = pairing === nothing ?
               (length(reference) == 1 ? [(1, index) for index in eachindex(result)] :
                throw(ArgumentError("multiple reference points require explicit pairing"))) :
               collect(pairing)
    for pair in selected
        validate(pair, Engine.compare, reference)
    end
    length(selected) == length(result) && sort(last.(selected)) == collect(eachindex(result)) ||
        throw(ArgumentError("pairing must reference valid points and select each result exactly once"))
    sort!(selected; by = last)
    for (left_index, right_index) in selected
        left, right = reference[left_index], result[right_index]
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
        b = observe(result, identity..., right_index)
        error = Engine.compare(a, b, identity[2]; frequencies = frequencies(left),
            result_basis = basis(left), kwargs...)
        return Engine.RMSError{Base.nonmissingtype(eltype(observe(error, Engine.absolute_error)))}(
            observe(error, Engine.absolute_error), observe(error, Engine.relative_error);
            details = ComputationDetails(merge(details(error).data, (; request = identity,
                estimators=(
                    reference=reference isa MonteCarloResult ? :empirical : :first_order,
                    result=result isa MonteCarloResult ? :empirical : :first_order)))))
    end
end

# Each pairing pair holds integer reference and result positions, not Bools, and its
# reference position lies within the reference.
function validate(pair::Tuple{Integer, Integer}, ::typeof(Engine.compare),
        reference::Union{MonteCarloResult, LinearErrorResult})
    1 <= first(pair) <= length(reference) || throw(ArgumentError(
        "pairing must reference valid points and select each result exactly once; got $pair"))
    return pair
end
function validate(pair::Union{Tuple{Bool, Integer}, Tuple{Integer, Bool}},
        ::typeof(Engine.compare), reference::Union{MonteCarloResult, LinearErrorResult})
    throw(ArgumentError(
        "pairing must reference valid points and select each result exactly once; got $pair"))
end
function validate(pair, ::typeof(Engine.compare),
        reference::Union{MonteCarloResult, LinearErrorResult})
    throw(ArgumentError("pairing must reference valid points and select each result exactly " *
                        "once; got $(repr(pair))"))
end
