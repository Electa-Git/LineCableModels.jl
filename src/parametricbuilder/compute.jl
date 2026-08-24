struct CableConstantsMaterializer end

function (::CableConstantsMaterializer)(design, separation, earth_resistivity)
    return Engine.CableConstantsProblem(design; separation, earth_resistivity)
end

function Engine.CableConstantsProblem(
        design::Gridspace{<:DataModel.CableDesign};
        separation = nothing,
        earth_resistivity = 100.0,
        combine::Symbol = :product
)
    return Gridspace{Engine.CableConstantsProblem}(
        CableConstantsMaterializer(),
        (
            _gridspace_axis(design),
            _gridspace_axis(separation),
            _gridspace_axis(earth_resistivity)
        ),
        (:design, :separation, :earth_resistivity);
        combine
    )
end

mutable struct _ManifestState
    automatic_keys::IdDict{Any, Int}
end

_ManifestState() = _ManifestState(IdDict{Any, Int}())

function _manifest_tree(value, state::_ManifestState)
    if value isa AutomaticGridKey
        label = get!(state.automatic_keys, value.token) do
            length(state.automatic_keys) + 1
        end
        return (automatic_grid = label,)
    elseif value isa NamedGridKey
        return (named_grid = _manifest_tree(value.value, state),)
    elseif value isa Type
        return string(value)
    elseif value isa Function
        names = fieldnames(typeof(value))
        fields = map(name -> _manifest_tree(getfield(value, name), state), names)
        return (function_type = string(typeof(value)), fields = NamedTuple{names}(fields))
    elseif value isa NamedTuple
        return NamedTuple{keys(value)}(map(item -> _manifest_tree(item, state), values(value)))
    elseif value isa AbstractDict
        entries = [(_manifest_tree(key, state), _manifest_tree(item, state))
                   for (key, item) in pairs(value)]
        sort!(entries; by = entry -> repr(first(entry)))
        return (type = string(typeof(value)), entries = tuple(entries...))
    elseif value isa AbstractSet
        entries = [_manifest_tree(item, state) for item in value]
        sort!(entries; by = repr)
        return (type = string(typeof(value)), entries = tuple(entries...))
    elseif value isa Tuple
        return map(item -> _manifest_tree(item, state), value)
    elseif value isa AbstractArray
        return (
            type = string(typeof(value)),
            size = size(value),
            values = map(item -> _manifest_tree(item, state), Tuple(value))
        )
    elseif value isa Union{Nothing, Missing, Bool, Number, Symbol, AbstractString, Char}
        return value
    elseif isstructtype(typeof(value))
        names = fieldnames(typeof(value))
        fields = map(name -> _manifest_tree(getfield(value, name), state), names)
        return (type = string(typeof(value)), fields = NamedTuple{names}(fields))
    end
    return (type = string(typeof(value)), representation = repr(value))
end

_manifest_tree(value) = _manifest_tree(value, _ManifestState())

function _configuration_failure(configuration, exception)
    return (
        coordinate = configuration_manifest(configuration),
        exception,
        message = sprint(showerror, exception)
    )
end

function _skippable_configuration_error(exception)
    exception isa Union{
        ArgumentError, AssertionError, DimensionMismatch, DomainError
    }
end

_append_result(::Nothing, value) = typeof(value)[value]

function _append_result(values::Vector{T}, value) where {T}
    value isa T && (push!(values, value); return values)
    W = typejoin(T, typeof(value))
    W === Any && throw(ArgumentError(
        "configuration results do not share one primitive result grammar",
    ))
    widened = Vector{W}(undef, length(values) + 1)
    copyto!(widened, values)
    widened[end] = value
    return widened
end

function _details(; failures, manifest)
    return Dict{Symbol, NamedTuple}(
        :failures => (entries = failures,),
        :manifest => manifest
    )
end

function _coupling_manifest(binding_sets)
    state = _ManifestState()
    entries = map(binding_sets) do bindings
        map(bindings) do binding
            (
                key = _manifest_tree(binding.key, state),
                index = binding.index,
                cardinality = binding.cardinality
            )
        end
    end
    return (entries = entries,)
end

function _manifest(formulation, options, binding_sets; random = nothing)
    return (
        backend = _manifest_tree(formulation),
        options = _manifest_tree(options),
        random,
        coupling = _coupling_manifest(binding_sets)
    )
end

function _traverse(problem::ParametricProblem, formulation, result_type)
    values = nothing
    resolved = NamedTuple[]
    resolved_bindings = Tuple[]
    failures = NamedTuple[]
    for configuration in configurations(problem.space)
        parameterization = configuration_manifest(configuration)
        try
            materialized = materialize(configuration)
            values = _append_result(
                values,
                compute(materialized, formulation.inner; options = problem.options)
            )
            push!(resolved, parameterization)
            push!(resolved_bindings, configuration.bindings)
        catch exception
            exception isa InterruptException && rethrow()
            (formulation.invalid === :error ||
             !_skippable_configuration_error(exception)) && rethrow()
            push!(failures, _configuration_failure(configuration, exception))
        end
    end
    typed_values = values === nothing ? AbstractProblemResult[] : values
    manifest_value = _manifest(formulation.inner, problem.options, resolved_bindings)
    return result_type(
        formulation,
        typed_values,
        resolved,
        _details(; failures, manifest = manifest_value)
    )
end

function compute(problem::ParametricProblem, formulation::Combinatorial)
    _traverse(problem, formulation, ParametricResult)
end
