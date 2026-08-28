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
    combine in (:product, :zip) || throw(ArgumentError(
        "combine must be :product or :zip"
    ))
    values = (rho, eps_r, mu_r, thickness)
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        return Gridspace{EarthProps.EarthModel}(_earth_model, grids; combine)
    end
    return _earth_model(values...)
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
    combine in (:product, :zip) || throw(ArgumentError(
        "combine must be :product or :zip"
    ))
    values = (
        designs,
        placements,
        connections,
        environment,
        system_id,
        line_length
    )
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        caller = (selected...) -> build(DataModel.LineCableSystem, selected...)
        return Gridspace{DataModel.LineCableSystem}(caller, grids; combine)
    end
    return build(DataModel.LineCableSystem, values...; combine)
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
