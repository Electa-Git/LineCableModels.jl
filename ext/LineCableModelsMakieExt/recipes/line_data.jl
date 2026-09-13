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

function _validate_plot_ydata(object, ydata)
    ydata isa Tuple || throw(ArgumentError("ydata must be a tuple"))
    isempty(ydata) && throw(ArgumentError(
        "at least one explicit observable request is required",
    ))
    all(request -> request isa Tuple, ydata) || throw(ArgumentError(
        "line plots accept explicit observable request tuples",
    ))
    validate_observables(object, ydata)
    all(ydata) do scientific_request
        expected = _diagonal_request(scientific_request) ? 2 : 3
        length(request_indices(scientific_request)) == expected
    end || throw(ArgumentError(
        "line plots require mode/frequency indices for diagonal requests and " *
        "row/column/frequency indices otherwise",
    ))
    return ydata
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

function _publish_request(object, request, target, input)
    publication = observables(object, (request,); units=(target,),
        clip=input.clip, atol=input.atol, frequencies=input.frequencies)
    observation = only(publication)
    contract = getproperty(publication.metadata.observation_columns,
        Symbol(Units.symbol(observation.quantity)))
    return (; observation, resolution=contract.resolution)
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
    selector = identity isa Tuple ? first(identity) : identity
    if (object isa SeriesImpedance && selector === L) ||
            (object isa ShuntAdmittance && selector === C)
        prefix = (prefix..., input.frequencies)
    end
    _diagonal_request(request) && return (prefix..., rows, samples)
    return (prefix..., rows, columns, samples)
end

function _publish_line_source(object, input, ydata)
    coordinates = map(request -> _request_coordinates(object, request), ydata)
    sample_indices = last.(coordinates)
    all(==(first(sample_indices)), sample_indices) || throw(DimensionMismatch(
        "all requests on one line dashboard must select the same frequency indices",
    ))
    frequency = _published_frequency(object, input, first(sample_indices))
    targets = unit_targets(
        ydata,
        basis(object);
        length_prefix = input.length_unit,
        overrides = input.quantity_units
    )
    publications = map(ydata, targets, coordinates) do request, target, indices
        publication = _publish_request(
            object,
            _materialized_line_request(object, input, request),
            target,
            input
        )
        _diagonal_request(request) || return publication
        observation = publication.observation
        rows, _, samples = indices
        values = reshape(observation.values, length(rows), 1, length(samples))
        return merge(publication, (observation=merge(observation, (; values)),))
    end
    observations = map(publication -> publication.observation, publications)
    resolutions = map(publication -> publication.resolution, publications)
    all(observation -> size(observation.values, 3) == length(frequency.values), observations) ||
        throw(DimensionMismatch("frequency count does not match line-parameter samples"))
    return (; frequency, observations, coordinates, resolutions)
end

function _prepare_line_observations(
        object::Union{LineParameters, SeriesImpedance, ShuntAdmittance};
        frequencies = nothing,
        ydata,
        freq_unit = :base,
        length_unit = :kilo,
        quantity_units = nothing,
        clip::Bool = true,
        atol = nothing
)
    _validate_plot_ydata(object, ydata)
    atol isa Real && length(unique(request_identity.(ydata))) > 1 && throw(ArgumentError(
        "a scalar atol requires one plotted quantity; use keyed native-unit tolerances"))
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
        any(request -> request_identity(request) in (L, C), ydata) &&
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
        clip,
        atol
    )
    published = _publish_line_source(object, input, ydata)
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

function _prepare_line_observations(source::Union{LineCableModels.AbstractUncertaintyResult,ObservationPublication};
        point::Integer,ydata,freq_unit=:base,length_unit=:kilo,quantity_units=nothing,
        clip::Bool=true,atol=nothing,frequencies=nothing,sample_indices=nothing)
    retained=source isa ObservationPublication
    retained && point!=1 && throw(ArgumentError("a retained publication contains one selected point"))
    dimensions=retained ? size(first(source).values) : size(observe(source[point],Z))
    coordinates=map(ydata) do request
        indices=request_indices(request)
        if length(indices)==4
            first(indices) isa Colon || first(indices)==point || throw(ArgumentError(
                "statistical request point differs from the selected comparison; select its pair explicitly"))
            indices=Base.tail(indices)
        end
        isempty(indices) && (indices=(Colon(),Colon(),Colon()))
        length(indices)==3 || throw(ArgumentError("UQ matrix plots require row, column and frequency indices"))
        selected=map(observation_indices,indices,dimensions)
        sample_indices===nothing && return selected
        frequency=filter(in(sample_indices),last(selected))
        isempty(frequency) && throw(ArgumentError("the selected band and request have no common frequency samples"))
        return (selected[1],selected[2],frequency)
    end
    all(indices -> last(indices)==last(first(coordinates)),coordinates) || throw(DimensionMismatch("UQ plot frequency selections differ"))
    targets=unit_targets(ydata,basis(source);length_prefix=length_unit,overrides=quantity_units)
    publications=map(ydata,targets,coordinates) do request,target,indices
        identity=request_identity(request)
        identity isa Tuple && length(identity)==3 && first(identity)===LineCableModels.statistics ||
            throw(ArgumentError("UQ matrix plots require an explicit selected statistic"))
        publication=observables(source,((identity...,point,indices...),);units=(target,),clip,atol)
        observation=only(publication)
        eltype(observation.values) <: Real || throw(ArgumentError(
            "complex statistical plots require explicit real-valued quantities, such as R/X or G/B"))
        contract=getproperty(publication.metadata.observation_columns,Symbol(Units.symbol(observation.quantity)))
        (;observation,resolution=contract.resolution)
    end
    f=if retained
        selected=unique(source.columns.frequency)
        contract=get(source.metadata.observation_columns,:frequency,nothing)
        contract===nothing ? selected : selected.*Units.scale_factor(contract.unit,Units.units(:base,:hertz))
    else
        LineCableModels.frequencies(source[point])
    end
    frequencies===nothing || frequencies==f || throw(ArgumentError("supplied UQ frequencies differ"))
    frequency=_frequency_observation(f[last(first(coordinates))],Units.units(freq_unit,:hertz))
    return (;frequency,observations=map(value -> value.observation,publications),coordinates,
        resolutions=map(value -> value.resolution,publications))
end

function _axis_scales(values; signed_log::Bool=false)
    _supports_log_values(values) && return (:linear, :log10)
    # Benchmark matrices commonly contain negative mutual terms. Preserve their
    # signs and zeros instead of suppressing the page-level log-y control.
    signed_log && !isempty(values) && all(value -> begin
        sample=nominal(value)
        sample isa Real && isfinite(sample) && isfinite(abs(uncertainty(value)))
    end,values) && return (:linear, :pseudolog10)
    return (:linear,)
end
