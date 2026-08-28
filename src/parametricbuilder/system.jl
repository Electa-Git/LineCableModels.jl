function _earth_model(rho, eps_r, mu_r, thickness)
    EarthProps.EarthModel(rho, eps_r, mu_r; thickness)
end
_earth_model(rho, eps_r, mu_r, ::Nothing) = EarthProps.EarthModel(rho, eps_r, mu_r)

"Construct static earth properties directly or as an explicit finite space."
function Earth(;
        rho,
        eps_r = 1.0,
        mu_r = 1.0,
        thickness = nothing,
        combine::Symbol = :product
)
    values = (rho, eps_r, mu_r, thickness)
    return _construction(EarthProps.EarthModel, _earth_model, values; combine)
end

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
    return _construction(DataModel.LineCableSystem, caller, values; combine)
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
    all(placements) do placement
        placement isa NamedTuple &&
            keys(placement) == (:design, :pose, :connections) &&
            placement.design isa DataModel.CableDesign &&
            placement.pose isa DataModel.Pose2
    end || throw(ArgumentError(
        "each placed cable requires design, pose, and connections"
    ))
    return build(
        DataModel.LineCableSystem,
        getproperty.(placements, :design),
        getproperty.(placements, :pose),
        getproperty.(placements, :connections),
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
    return _construction(
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
    return _construction(
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
    caller = (selected...) -> Engine.LineParametersProblem(selected...)
    return Gridspace{Engine.LineParametersProblem}(caller, sources; combine)
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
