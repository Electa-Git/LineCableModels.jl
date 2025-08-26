"""
$(TYPEDSIGNATURES)

Inspects all numerical data within a `LineParametersProblem` and determines the
common floating-point type. If any value (frequencies, geometric properties,
material properties, or earth properties) is a `Measurement`, the function
returns `Measurement{Float64}`. Otherwise, it returns `Float64`.
"""
function _find_common_type(problem::LineParametersProblem)
    # Check frequencies
    any(x -> x isa Measurement, problem.frequencies) && return Measurement{Float64}

    # Check cable system properties
    for cable in problem.system.cables
        (cable.horz isa Measurement || cable.vert isa Measurement) && return Measurement{Float64}
        for component in cable.design_data.components
            if any(x -> x isa Measurement, (
                component.conductor_group.radius_in, component.conductor_group.radius_ext,
                component.insulator_group.radius_in, component.insulator_group.radius_ext,
                component.conductor_props.rho, component.conductor_props.mu_r, component.conductor_props.eps_r,
                component.insulator_props.rho, component.insulator_props.mu_r, component.insulator_props.eps_r,
                component.insulator_group.shunt_capacitance, component.insulator_group.shunt_conductance
            ))
                return Measurement{Float64}
            end
        end
    end

    # Check earth model properties
    if !isnothing(problem.earth_props)
        for layer in problem.earth_props.layers
            if any(x -> x isa Measurement, (layer.rho_g, layer.mu_g, layer.eps_g))
                return Measurement{Float64}
            end
        end
    end

    if !isnothing(problem.temperature)
        if problem.temperature isa Measurement
            return Measurement{Float64}
        end

    end

    return Float64
end

function _get_earth_data(formulation::AbstractEHEMFormulation, earth_model::EarthModel, frequencies::Vector{<:REALSCALAR}, T::DataType)
    return formulation(earth_model, frequencies, T)
end

"""
Default method for when no EHEM formulation is provided. It processes the
original multi-layer earth model, promoting its property vectors to the
target numeric type `T` and returning it in the standard "Array of Structs" format.
"""
function _get_earth_data(::Nothing, earth_model::EarthModel, frequencies::Vector{<:REALSCALAR}, T::DataType)
    return [
        (
            rho_g=T.(layer.rho_g),
            eps_g=T.(layer.eps_g),
            mu_g=T.(layer.mu_g)
        )
        for layer in earth_model.layers
    ]
end