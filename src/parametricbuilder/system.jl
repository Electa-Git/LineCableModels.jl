function build(
        ::Type{DataModel.LineCableSystem},
        designs,
        placements;
        connections = nothing,
        environment = nothing,
        system_id = "line-cable-system",
        line_length = 1,
        combine::Symbol = :product
)
    values = (
        designs,
        placements,
        connections,
        environment,
        system_id,
        line_length
    )
    caller = (selected...) -> build(DataModel.LineCableSystem, selected...)
    return parameterize(DataModel.LineCableSystem, caller, values; combine)
end

function _placed_system(
        placements,
        environment,
        system_id,
        line_length
)
    placements isa Union{Tuple, AbstractVector} && !isempty(placements) || throw(
        ArgumentError("placed cable declarations must be a nonempty collection")
    )
    declarations = Any[]
    for placement in placements
        if placement isa Union{Tuple, AbstractVector}
            append!(declarations, placement)
        else
            push!(declarations, placement)
        end
    end
    isempty(declarations) && throw(ArgumentError(
        "placed cable declarations must be a nonempty collection"
    ))
    all(declarations) do placement
        placement isa NamedTuple &&
            keys(placement) == (:design, :pose, :connections) &&
            placement.design isa DataModel.CableDesign &&
            placement.pose isa DataModel.Pose2
    end || throw(ArgumentError(
        "each placed cable requires design, pose, and connections"
    ))
    return build(
        DataModel.LineCableSystem,
        getproperty.(declarations, :design),
        getproperty.(declarations, :pose),
        getproperty.(declarations, :connections),
        environment,
        system_id,
        line_length
    )
end

function build(
        ::Type{DataModel.LineCableSystem},
        placements;
        environment = nothing,
        system_id = "line-cable-system",
        line_length = 1,
        combine::Symbol = :product
)
    values = (placements, environment, system_id, line_length)
    return parameterize(
        DataModel.LineCableSystem, _placed_system, values; combine
    )
end

function _line_problem(system, temperature, earth, frequencies)
    return Engine.LineParametersProblem(
        system;
        temperature,
        earth_props = earth,
        frequencies
    )
end

function Engine.LineParametersProblem(
        system::Gridspace{DataModel.LineCableSystem};
        temperature = 20,
        earth_props,
        frequencies = [50],
        combine::Symbol = :product
)
    values = (system, temperature, earth_props, frequencies)
    grids = map(values) do value
        value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
    end
    return Gridspace{Engine.LineParametersProblem}(_line_problem, grids; combine)
end

function Engine.LineParametersProblem(
        system::Gridspace{DataModel.LineCableSystem},
        earth_props::Union{AbstractGrid, Gridspace};
        temperature = 20,
        frequencies = [50],
        combine::Symbol = :product
)
    return Engine.LineParametersProblem(
        system;
        temperature,
        earth_props,
        frequencies,
        combine
    )
end

function _placed_line_problem(
        placements,
        environment,
        system_id,
        line_length,
        temperature,
        earth_props,
        frequencies
)
    system = _placed_system(placements, environment, system_id, line_length)
    return Engine.LineParametersProblem(
        system;
        temperature,
        earth_props,
        frequencies
    )
end

function Engine.LineParametersProblem(
        placements::Union{Tuple, AbstractVector, Gridspace};
        environment = nothing,
        system_id = "line-cable-system",
        line_length = 1,
        temperature = 20,
        earth_props,
        frequencies = [50],
        combine::Symbol = :product
)
    values = (
        placements, environment, system_id, line_length, temperature,
        earth_props, frequencies
    )
    return parameterize(
        Engine.LineParametersProblem, _placed_line_problem, values; combine
    )
end

function Engine.LineParametersProblem(
        system::DataModel.LineCableSystem,
        earth_props::Union{AbstractGrid, Gridspace};
        temperature = 20,
        frequencies = [50],
        combine::Symbol = :product
)
    values = (system, temperature, earth_props, frequencies)
    grids = map(values) do value
        value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
    end
    return Gridspace{Engine.LineParametersProblem}(_line_problem, grids; combine)
end

function Engine.LineParametersProblem(
        designs,
        placements,
        connections,
        environment,
        system_id,
        line_length,
        temperature,
        earth_props,
        frequencies;
        combine::Symbol = :product
)
    values = (
        designs,
        placements,
        connections,
        environment,
        system_id,
        line_length,
        temperature,
        earth_props,
        frequencies
    )
    any(value -> value isa Union{AbstractGrid, Gridspace}, values) || throw(
        MethodError(Engine.LineParametersProblem, values)
    )
    sources = map(values) do value
        value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
    end
    return Gridspace{Engine.LineParametersProblem}(
        _declared_line_problem, sources; combine
    )
end

function _declared_line_problem(
        designs,
        placements,
        connections,
        environment,
        system_id,
        line_length,
        temperature,
        earth_props,
        frequencies
)
    return Engine.LineParametersProblem(
        designs,
        placements,
        connections,
        environment,
        system_id,
        line_length,
        temperature,
        earth_props,
        frequencies
    )
end

function Engine.LineParametersProblem(
        designs,
        placements;
        connections = nothing,
        environment = nothing,
        system_id = "line-cable-system",
        line_length = 1,
        temperature = 20,
        earth_props,
        frequencies = [50],
        combine::Symbol = :product
)
    return Engine.LineParametersProblem(
        designs,
        placements,
        connections,
        environment,
        system_id,
        line_length,
        temperature,
        earth_props,
        frequencies;
        combine
    )
end
