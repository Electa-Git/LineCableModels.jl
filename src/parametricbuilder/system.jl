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

function _system_space(
        designs::Tuple;
        positions,
        connections,
        environment,
        system_id,
        line_length,
        combine
)
    count = length(designs)
    count > 0 || throw(ArgumentError("a line system requires one cable design"))
    values = (
        designs...,
        positions,
        connections,
        environment,
        String(system_id),
        line_length
    )
    build = function (resolved...)
        selected = DataModel.CableDesign[resolved[index] for index in 1:count]
        offset = count
        return DataModel.LineCableSystem(
            selected;
            positions = resolved[offset + 1],
            connections = resolved[offset + 2],
            environment = resolved[offset + 3],
            system_id = resolved[offset + 4],
            line_length = resolved[offset + 5]
        )
    end
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        return Gridspace{DataModel.LineCableSystem}(build, grids; combine)
    end
    return build(values...)
end

function LineCableSystem(
        designs::Tuple;
        positions,
        connections = nothing,
        environment = nothing,
        system_id::AbstractString = "line-cable-system",
        line_length = 1,
        combine::Symbol = :product
)
    all(design -> design isa Union{DataModel.CableDesign,
            Gridspace{DataModel.CableDesign}}, designs) || throw(ArgumentError(
        "designs must contain CableDesign values or their Gridspaces"
    ))
    return _system_space(
        designs;
        positions,
        connections,
        environment,
        system_id,
        line_length,
        combine
    )
end

function LineCableSystem(
        design::Gridspace{DataModel.CableDesign};
        position = DataModel.Pose2(0, -1, 0),
        connections = nothing,
        kwargs...
)
    return LineCableSystem(
        (design,);
        positions = position,
        connections,
        kwargs...
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
