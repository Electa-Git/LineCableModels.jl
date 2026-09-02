# Plotting-only request normalization and detached data preparation. Scientific
# values continue to be produced by the public observation API; this file owns
# only the Makie recipe's convenience policy.

_line_request_family(request) = Units.family(request_quantity(request))

function _diagonal_request(request)
    identity = request_identity(request)
    return identity isa Tuple && last(identity) === diag
end

_family_parent(::Val{:series}) = Z
_family_parent(::Val{:shunt}) = Y

function _validate_plot_requests(object, requests)
    requests isa Tuple || throw(ArgumentError("requests must be a tuple"))
    isempty(requests) && throw(ArgumentError(
        "at least one explicit observable request is required",
    ))
    all(request -> request isa Tuple, requests) || throw(ArgumentError(
        "line plots accept explicit observable request tuples",
    ))
    validate_observables(object, requests)
    all(requests) do scientific_request
        expected = _diagonal_request(scientific_request) ? 2 : 3
        length(request_indices(scientific_request)) == expected
    end || throw(ArgumentError(
        "line plots require mode/frequency indices for diagonal requests and " *
        "row/column/frequency indices otherwise",
    ))
    return requests
end

function _frequency_observation(values, target)
    native = Units.units(:base, :hertz)
    factor = Units.scale_factor(native, target)
    return (;
        values = map(value -> value * factor, values),
        quantity = Units.Quantity{:frequency}(),
        unit = target
    )
end

function _published_frequency(object, input, selector)
    target = Units.units(input.freq_unit, :hertz)
    if !(object isa LineParameters)
        return _frequency_observation(input.frequencies[selector], target)
    end
    published = only(observables(
        object,
        ((frequencies, selector),);
        units = (target,),
        clip = input.clip
    ))
    if input.frequencies !== nothing
        supplied = _frequency_observation(input.frequencies[selector], target)
        supplied.values == published.values || throw(ArgumentError(
            "supplied frequencies do not match the LineParameters frequency axis",
        ))
    end
    return published
end

function _publish_request(object, request, target, clip::Bool)
    return only(observables(object, (request,); units = (target,), clip))
end

function _request_coordinates(object, request)
    dimensions = object isa LineParameters ? size(Z(object)) : size(object)
    if _diagonal_request(request)
        mode, frequency = request_indices(request)
        modes = observation_indices(mode, dimensions[1])
        samples = observation_indices(frequency, dimensions[3])
        return modes, [1], samples
    end
    row, column, frequency = request_indices(request)
    rows = observation_indices(row, dimensions[1])
    columns = observation_indices(column, dimensions[2])
    frequency_count = object isa LineParameters ? nfrequencies(object) : size(object, 3)
    samples = observation_indices(frequency, frequency_count)
    return rows, columns, samples
end

function _materialized_line_request(object, input, request)
    rows, columns, samples = _request_coordinates(object, request)
    identity = observation_request(object, request).identity
    prefix = identity isa Tuple ? identity : (identity,)
    _diagonal_request(request) && return (prefix..., rows, samples)
    if object isa SeriesImpedance && identity === L
        prefix = (L, input.frequencies)
    elseif object isa ShuntAdmittance && identity === C
        prefix = (C, input.frequencies)
    end
    return (prefix..., rows, columns, samples)
end

function _publish_line_source(object, input, requests)
    coordinates = map(request -> _request_coordinates(object, request), requests)
    sample_indices = last.(coordinates)
    all(==(first(sample_indices)), sample_indices) || throw(DimensionMismatch(
        "all requests on one line dashboard must select the same frequency indices",
    ))
    frequency = _published_frequency(object, input, first(sample_indices))
    targets = unit_targets(
        requests,
        basis(object);
        length_prefix = input.length_unit,
        overrides = input.quantity_units
    )
    observations = map(requests, targets, coordinates) do request, target, indices
        observation = _publish_request(
            object,
            _materialized_line_request(object, input, request),
            target,
            input.clip
        )
        _diagonal_request(request) || return observation
        rows, _, samples = indices
        values = reshape(observation.values, length(rows), 1, length(samples))
        return merge(observation, (; values))
    end
    all(observation -> size(observation.values, 3) == length(frequency.values), observations) ||
        throw(DimensionMismatch("frequency count does not match line-parameter samples"))
    return (; frequency, observations, coordinates)
end

function _prepare_line_observations(
        object::Union{LineParameters, SeriesImpedance, ShuntAdmittance};
        frequencies = nothing,
        requests,
        freq_unit = :base,
        length_unit = :kilo,
        quantity_units = nothing,
        clip::Bool = true
)
    _validate_plot_requests(object, requests)
    supplied = frequencies === nothing ? nothing : collect(frequencies)
    if object isa Union{SeriesImpedance, ShuntAdmittance}
        supplied === nothing && throw(ArgumentError(
            "frequencies are required for SeriesImpedance and ShuntAdmittance",
        ))
        length(supplied) == size(object, 3) || throw(DimensionMismatch(
            "frequency vector length does not match the parameter depth",
        ))
    end
    if supplied !== nothing
        all(isfinite, supplied) || throw(ArgumentError("frequencies must be finite"))
        any(request -> request_identity(request) in (L, C), requests) &&
            any(iszero, supplied) &&
            throw(DomainError(
                supplied,
                "inductance and capacitance are undefined at zero frequency"
            ))
    end
    input = (;
        frequencies = supplied,
        freq_unit,
        length_unit,
        quantity_units,
        clip
    )
    published = _publish_line_source(object, input, requests)
    length(published.frequency.values) <= 1 &&
        @warn "Frequency vector has $(length(published.frequency.values)) sample(s); nothing to plot."
    return published
end

function _supports_log_values(samples)
    found = false
    samples === nothing && return false
    for sample in samples
        found = true
        value = nominal(sample)
        uncertainty_value = abs(uncertainty(sample))
        value isa Real && isfinite(value) && isfinite(uncertainty_value) &&
        value - uncertainty_value > 0 || return false
    end
    return found
end

_axis_scales(values) = _supports_log_values(values) ? (:linear, :log10) : (:linear,)

_line_families(::SeriesImpedance) = (Val(:series),)
_line_families(::ShuntAdmittance) = (Val(:shunt),)
_line_families(::LineParameters) = (Val(:series), Val(:shunt))

function _line_axis_groups(requests, selected_indices = eachindex(requests))
    keys = Any[]
    groups = Vector{Vector{Int}}()
    for index in selected_indices
        key = request_quantity(requests[index])
        group_index = findfirst(candidate -> candidate == key, keys)
        if group_index === nothing
            push!(keys, key)
            push!(groups, Int[index])
        else
            push!(groups[group_index], index)
        end
    end
    return Tuple(Tuple(group) for group in groups)
end
