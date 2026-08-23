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

function CalculationManifest(
        resolved_parameterization,
        problem_assumptions,
        formulation,
        execution_policy,
        calculation_options
)
    solver = string(typeof(formulation))
    payload = (
        resolved_parameterization = resolved_parameterization,
        problem_assumptions = problem_assumptions,
        formulation = _manifest_tree(formulation),
        solver,
        execution_policy = _manifest_tree(execution_policy),
        calculation_options = _manifest_tree(calculation_options)
    )
    digest = bytes2hex(sha256(_stable_bytes(payload)))
    return CalculationManifest(
        payload.resolved_parameterization,
        payload.problem_assumptions,
        payload.formulation,
        payload.solver,
        payload.execution_policy,
        payload.calculation_options,
        digest
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
        sort!(entries; by = entry -> String(_stable_bytes(first(entry))))
        return (type = string(typeof(value)), entries = tuple(entries...))
    elseif value isa AbstractSet
        entries = [_manifest_tree(item, state) for item in value]
        sort!(entries; by = entry -> String(_stable_bytes(entry)))
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

function _stable_bytes(value)
    io = IOBuffer()
    _write_stable(io, value)
    return take!(io)
end

function _write_stable(io::IO, value)
    if value isa NamedTuple
        print(io, "named{")
        for key in keys(value)
            print(io, repr(key), '=')
            _write_stable(io, getproperty(value, key))
            print(io, ';')
        end
        print(io, '}')
    elseif value isa Tuple
        print(io, "tuple[")
        for item in value
            _write_stable(io, item)
            print(io, ';')
        end
        print(io, ']')
    elseif value isa AbstractArray
        print(io, "array", repr(size(value)), '[')
        for item in value
            _write_stable(io, item)
            print(io, ';')
        end
        print(io, ']')
    elseif value isa AbstractFloat
        print(io, string(typeof(value)), ':', bitstring(value))
    elseif value isa Complex
        print(io, "complex(")
        _write_stable(io, real(value))
        print(io, ',')
        _write_stable(io, imag(value))
        print(io, ')')
    elseif value isa Number
        print(io, string(typeof(value)), ':', repr(value))
    else
        print(io, string(typeof(value)), ':', repr(value))
    end
    return io
end

function _configuration_failure(index, configuration, exception)
    return ConfigurationFailure(
        index,
        configuration_manifest(configuration),
        string(typeof(exception)),
        sprint(showerror, exception)
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

function _details(; failures, manifest_value, samples_value = nothing,
        histograms_value = nothing, random_value = nothing)
    return Dict{Symbol, NamedTuple}(
        :failures => (values = failures,),
        :samples => (values = samples_value,),
        :histograms => (values = histograms_value,),
        :random => random_value === nothing ? (value = nothing,) : random_value,
        :manifest => (value = manifest_value,)
    )
end

function _manifest(problem, formulation, resolved, policy, options)
    return CalculationManifest(
        tuple(resolved...),
        _manifest_tree(problem.space),
        formulation,
        policy,
        options
    )
end

function _traverse(problem::ParametricProblem, formulation, result_type)
    values = nothing
    resolved = NamedTuple[]
    failures = ConfigurationFailure[]
    for (index, configuration) in enumerate(configurations(problem.space))
        parameterization = configuration_manifest(configuration)
        try
            materialized = materialize(configuration)
            values = _append_result(
                values,
                compute(materialized, formulation.inner; options = problem.options)
            )
            push!(resolved, parameterization)
        catch exception
            exception isa InterruptException && rethrow()
            (formulation.invalid === :error ||
             !_skippable_configuration_error(exception)) && rethrow()
            push!(failures, _configuration_failure(index, configuration, exception))
        end
    end
    typed_values = values === nothing ? AbstractProblemResult[] : values
    manifest_value = _manifest(
        problem, formulation.inner, resolved, formulation, problem.options)
    return result_type(
        formulation,
        typed_values,
        problem.space,
        _details(; failures, manifest_value)
    )
end

function compute(problem::ParametricProblem, formulation::Combinatorial)
    _traverse(problem, formulation, ParametricResult)
end
