import DataFrames: DataFrame

"""
$(TYPEDSIGNATURES)

Convert a cable design to a tabular description.

# Arguments

- `design`: Materialised cable design.
- `format`: Table layout. Use `:components` for equivalent component values,
  `:detailed` for one column per cable layer, or `:baseparams` for compact
  solver-input fields. Default: `:components`.

# Returns

- A `DataFrame` in the selected layout.

# Errors

- Throws an error when `format` is unsupported.
"""
function DataFrame(design::CableDesign, format::Symbol = :components)::DataFrame
    if format == :baseparams
        return DataFrame(_base_parameters(design))
    elseif format == :components
        # Component-level properties
        properties = [
            :radius_in_con,
            :radius_ext_con,
            :rho_con,
            :alpha_con,
            :mu_con,
            :radius_ext_ins,
            :eps_ins,
            :mu_ins,
            :loss_factor_ins
        ]

        # Initialise the DataFrame.
        data = DataFrame(property = properties)

        # Add one column per cable component.
        for component in design.components
            # Use the component identifier as the column name.
            col = component.id

            # Map the current component structure to the retained row names.
            # Calculate the dielectric loss factor.
            ω = 2 * π * component.insulator_group.reference_frequency
            C_eq = component.insulator_group.shunt_capacitance
            G_eq = component.insulator_group.shunt_conductance
            loss_factor = G_eq / (ω * C_eq)

            # Collect values in the retained row order.
            new_col = [
                component.conductor_group.r_in,               # radius_in_con
                component.conductor_group.r_ex,              # radius_ext_con
                component.conductor_props.rho,                     # rho_con
                component.conductor_props.alpha,                   # alpha_con
                component.conductor_props.mu_r,                    # mu_con
                component.insulator_group.r_ex,              # radius_ext_ins
                component.insulator_props.eps_r,                   # eps_ins
                component.insulator_props.mu_r,                    # mu_ins
                loss_factor                                       # loss_factor_ins
            ]

            # Add to DataFrame
            data[!, col] = new_col
        end

    elseif format == :detailed
        # Detailed part-by-part breakdown
        properties = [
            "type",
            "r_in",
            "r_ex",
            "diam_in",
            "diam_ext",
            "thickness",
            "cross_section",
            "num_wires",
            "resistance",
            "alpha",
            "gmr",
            "gmr/radius",
            "shunt_capacitance",
            "shunt_conductance"
        ]

        # Initialise the DataFrame.
        data = DataFrame(property = properties)

        # Process each component
        for component in design.components
            # Handle conductor group layers
            for (i, part) in enumerate(component.conductor_group.layers)
                # Column name with component ID and layer number
                col = lowercase(component.id) * ", cond. layer " * string(i)

                # Collect values for each property
                new_col = _extract_part_properties(part, properties)

                # Add to DataFrame
                data[!, col] = new_col
            end

            # Handle insulator group layers
            for (i, part) in enumerate(component.insulator_group.layers)
                # Column name with component ID and layer number
                col = lowercase(component.id) * ", ins. layer " * string(i)

                # Collect values for each property
                new_col = _extract_part_properties(part, properties)

                # Add to DataFrame
                data[!, col] = new_col
            end
        end
    else
        Base.error(
            "Unsupported format: $format. Use :components or :detailed",
        )
    end

    return data
end

"""
$(TYPEDSIGNATURES)

Render already-computed cable constants without performing a calculation.
"""
function DataFrame(constants::CableConstants)::DataFrame
    values = (R(constants), L(constants), C(constants))
    return DataFrame(
        parameter = ["R", "L", "C"],
        value = collect(values),
        unit = ["Ω/m", "H/m", "F/m"]
    )
end

"""
$(TYPEDSIGNATURES)

Return the fixed row values used by the `:detailed` cable-design table.

# Arguments

- `part`: Cable part to inspect.
- `properties`: Retained layout argument. The current row order is fixed.

# Returns

- Values for type, radii, diameters, thickness, cross-section, wire count,
  resistance, temperature coefficient, GMR, GMR-to-radius ratio, capacitance,
  and conductance. A field absent from `part` produces `missing`.
"""
function _extract_part_properties(part, properties)
    return [
        lowercase(string(typeof(part))),  # type
        hasfield(typeof(part), :r_in) ?
        getproperty(part, :r_in) : missing,
        hasfield(typeof(part), :r_ex) ?
        getproperty(part, :r_ex) : missing,
        hasfield(typeof(part), :r_in) ?
        2 * getproperty(part, :r_in) : missing,
        hasfield(typeof(part), :r_ex) ?
        2 * getproperty(part, :r_ex) : missing,
        hasfield(typeof(part), :r_ex) &&
        hasfield(typeof(part), :r_in) ?
        (getproperty(part, :r_ex) - getproperty(part, :r_in)) :
        missing,
        hasfield(typeof(part), :cross_section) ?
        getproperty(part, :cross_section) : missing,
        hasfield(typeof(part), :num_wires) ?
        getproperty(part, :num_wires) : missing,
        hasfield(typeof(part), :resistance) ?
        getproperty(part, :resistance) : missing,
        hasfield(typeof(part), :alpha) ||
        (
            hasfield(typeof(part), :material_props) &&
            hasfield(typeof(getproperty(part, :material_props)), :alpha)
        ) ?
        (
            hasfield(typeof(part), :alpha) ?
            getproperty(part, :alpha) :
            getproperty(getproperty(part, :material_props), :alpha)
        ) : missing,
        hasfield(typeof(part), :gmr) ?
        getproperty(part, :gmr) : missing,
        hasfield(typeof(part), :gmr) &&
        hasfield(typeof(part), :r_ex) ?
        (getproperty(part, :gmr) / getproperty(part, :r_ex)) : missing,
        hasfield(typeof(part), :shunt_capacitance) ?
        getproperty(part, :shunt_capacitance) : missing,
        hasfield(typeof(part), :shunt_conductance) ?
        getproperty(part, :shunt_conductance) : missing
    ]
end
