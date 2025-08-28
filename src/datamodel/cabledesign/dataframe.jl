import DataFrames: DataFrame


"""
$(TYPEDSIGNATURES)

Extracts and displays data from a [`CableDesign`](@ref).

# Arguments

- `design`: A [`CableDesign`](@ref) object to extract data from.
- `format`: Symbol indicating the level of detail:
  - `:baseparams`: Basic RLC parameters with nominal value comparison (default).
  - `:components`: Component-level equivalent properties.
  - `:detailed`: Individual cable part properties with layer-by-layer breakdown.
- `S`: Separation distance between cables \\[m\\] (only used for `:baseparams` format). Default: outermost cable diameter.
- `rho_e`: Resistivity of the earth \\[Ω·m\\] (only used for `:baseparams` format). Default: 100.

# Returns

- A `DataFrame` containing the requested cable data in the specified format.

# Examples

```julia
# Get basic RLC parameters
data = DataFrame(design)  # Default is :baseparams format

# Get component-level data
comp_data = DataFrame(design, :components)

# Get detailed part-by-part breakdown
detailed_data = DataFrame(design, :detailed)

# Specify earth parameters for core calculations
core_data = DataFrame(design, :baseparams, S=0.5, rho_e=150)
```

# See also

- [`CableDesign`](@ref)
- [`calc_tubular_resistance`](@ref)
- [`calc_inductance_trifoil`](@ref)
- [`calc_shunt_capacitance`](@ref)
"""
function DataFrame(
    design::CableDesign,
    format::Symbol=:baseparams;
    S::Union{Nothing,Number}=nothing,
    rho_e::Number=100.0,
)::DataFrame
    if format == :baseparams
        # Core parameters calculation
        # Get components from the vector
        if length(design.components) < 2
            throw(ArgumentError("At least two components are required for :baseparams format."))
        end

        cable_core = design.components[1]
        cable_shield = design.components[2]
        cable_outer = design.components[end]

        # Determine separation distance if not provided
        S =
            S === nothing ?
            (
                # Check if we need to use insulator or conductor radius
                isnan(cable_outer.insulator_group.radius_ext) ?
                2 * cable_outer.conductor_group.radius_ext :
                2 * cable_outer.insulator_group.radius_ext
            ) : S

        # Compute R, L, and C using given formulas - mapped to new data structure
        # Cable core resistance
        R =
            calc_tubular_resistance(
                cable_core.conductor_group.radius_in,
                cable_core.conductor_group.radius_ext,
                cable_core.conductor_props.rho,
                0.0, 20.0, 20.0,
            ) * 1e3

        # Inductance calculation
        L =
            calc_inductance_trifoil(
                cable_core.conductor_group.radius_in,
                cable_core.conductor_group.radius_ext,
                cable_core.conductor_props.rho,
                cable_core.conductor_props.mu_r,
                cable_shield.conductor_group.radius_in,
                cable_shield.conductor_group.radius_ext,
                cable_shield.conductor_props.rho,
                cable_shield.conductor_props.mu_r,
                S,
                rho_e=rho_e,
            ) * 1e6

        # Capacitance calculation
        C =
            calc_shunt_capacitance(
                cable_core.conductor_group.radius_ext,
                cable_core.insulator_group.radius_ext,
                cable_core.insulator_props.eps_r,
            ) * 1e6 * 1e3

        # Prepare nominal values from CableDesign
        nominals = [
            design.nominal_data.resistance,
            design.nominal_data.inductance,
            design.nominal_data.capacitance,
        ]

        # Calculate differences
        diffs = map(zip([R, L, C], nominals)) do (computed, nominal)
            if isnothing(nominal)
                return missing
            else
                return abs(nominal - computed) / nominal * 100
            end
        end

        # Compute the comparison DataFrame
        data = DataFrame(
            parameter=["R [Ω/km]", "L [mH/km]", "C [μF/km]"],
            computed=[R, L, C],
            nominal=nominals,
        )

        # Add percent_diff column only for rows with non-nothing nominal values
        data[!, "percent_diff"] = diffs

        # Handle measurement bounds if present
        has_error_bounds = !(isnan(to_lower(R)) || isnan(to_upper(R)))
        if has_error_bounds
            data[!, "lower"] = [to_lower(R), to_lower(L), to_lower(C)]
            data[!, "upper"] = [to_upper(R), to_upper(L), to_upper(C)]

            # Add compliance column only for rows with non-nothing nominal values
            data[!, "in_range?"] = map(zip(data.nominal, data.lower, data.upper)) do (nom, low, up)
                isnothing(nom) ? missing : (nom >= low && nom <= up)
            end
        end

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
            :loss_factor_ins,
        ]

        # Initialize the DataFrame
        data = DataFrame(property=properties)

        # Process each component - now using vector 
        for component in design.components
            # Use component ID as column name
            col = component.id

            # For each component, we need to map new structure to old column names
            # Calculate loss factor from resistivity
            ω = 2 * π * f₀  # Using default frequency
            C_eq = component.insulator_group.shunt_capacitance
            G_eq = component.insulator_group.shunt_conductance
            loss_factor = G_eq / (ω * C_eq)

            # Collect values for each property - mapping from new structure to old property names
            new_col = [
                component.conductor_group.radius_in,               # radius_in_con
                component.conductor_group.radius_ext,              # radius_ext_con
                component.conductor_props.rho,                     # rho_con
                component.conductor_props.alpha,                   # alpha_con
                component.conductor_props.mu_r,                    # mu_con
                component.insulator_group.radius_ext,              # radius_ext_ins
                component.insulator_props.eps_r,                   # eps_ins
                component.insulator_props.mu_r,                    # mu_ins
                loss_factor,                                       # loss_factor_ins
            ]

            # Add to DataFrame
            data[!, col] = new_col
        end

    elseif format == :detailed
        # Detailed part-by-part breakdown
        properties = [
            "type",
            "radius_in",
            "radius_ext",
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
            "shunt_conductance",
        ]

        # Initialize the DataFrame
        data = DataFrame(property=properties)

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
        error("Unsupported format: $format. Use :baseparams, :components, or :detailed")
    end

    return data
end

"""
$(TYPEDSIGNATURES)

Helper function to extract properties from a part for detailed format.

# Arguments

- `part`: An instance of [`AbstractCablePart`](@ref) from which to extract properties.
- `properties`: A vector of symbols indicating which properties to extract (not used in the current implementation).

# Returns

- A vector containing the extracted properties in the following order:
  - `type`: The lowercase string representation of the part's type.
  - `radius_in`: The inner radius of the part, if it exists, otherwise `missing`.
  - `radius_ext`: The outer radius of the part, if it exists, otherwise `missing`.
  - `diameter_in`: The inner diameter of the part (2 * radius_in), if `radius_in` exists, otherwise `missing`.
  - `diameter_ext`: The outer diameter of the part (2 * radius_ext), if `radius_ext` exists, otherwise `missing`.
  - `thickness`: The difference between `radius_ext` and `radius_in`, if both exist, otherwise `missing`.
  - `cross_section`: The cross-sectional area of the part, if it exists, otherwise `missing`.
  - `num_wires`: The number of wires in the part, if it exists, otherwise `missing`.
  - `resistance`: The resistance of the part, if it exists, otherwise `missing`.
  - `alpha`: The temperature coefficient of resistivity of the part or its material, if it exists, otherwise `missing`.
  - `gmr`: The geometric mean radius of the part, if it exists, otherwise `missing`.
  - `gmr_ratio`: The ratio of `gmr` to `radius_ext`, if both exist, otherwise `missing`.
  - `shunt_capacitance`: The shunt capacitance of the part, if it exists, otherwise `missing`.
  - `shunt_conductance`: The shunt conductance of the part, if it exists, otherwise `missing`.

# Notes

This function is used to create a standardized format for displaying detailed information about cable parts.

# Examples

```julia
part = Conductor(...)
properties = [:radius_in, :radius_ext, :resistance]  # Example of properties to extract
extracted_properties = _extract_part_properties(part, properties)
println(extracted_properties)
```
"""
function _extract_part_properties(part, properties)
    return [
        lowercase(string(typeof(part))),  # type
        hasfield(typeof(part), :radius_in) ?
        getfield(part, :radius_in) : missing,
        hasfield(typeof(part), :radius_ext) ?
        getfield(part, :radius_ext) : missing,
        hasfield(typeof(part), :radius_in) ?
        2 * getfield(part, :radius_in) : missing,
        hasfield(typeof(part), :radius_ext) ?
        2 * getfield(part, :radius_ext) : missing,
        hasfield(typeof(part), :radius_ext) &&
        hasfield(typeof(part), :radius_in) ?
        (getfield(part, :radius_ext) - getfield(part, :radius_in)) :
        missing,
        hasfield(typeof(part), :cross_section) ?
        getfield(part, :cross_section) : missing,
        hasfield(typeof(part), :num_wires) ?
        getfield(part, :num_wires) : missing,
        hasfield(typeof(part), :resistance) ?
        getfield(part, :resistance) : missing,
        hasfield(typeof(part), :alpha) ||
        (
            hasfield(typeof(part), :material_props) &&
            hasfield(typeof(getfield(part, :material_props)), :alpha)
        ) ?
        (
            hasfield(typeof(part), :alpha) ?
            getfield(part, :alpha) :
            getfield(getfield(part, :material_props), :alpha)
        ) : missing,
        hasfield(typeof(part), :gmr) ?
        getfield(part, :gmr) : missing,
        hasfield(typeof(part), :gmr) &&
        hasfield(typeof(part), :radius_ext) ?
        (getfield(part, :gmr) / getfield(part, :radius_ext)) : missing,
        hasfield(typeof(part), :shunt_capacitance) ?
        getfield(part, :shunt_capacitance) : missing,
        hasfield(typeof(part), :shunt_conductance) ?
        getfield(part, :shunt_conductance) : missing,
    ]
end