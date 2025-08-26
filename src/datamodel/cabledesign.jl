import ..LineCableModels: add!

"""
$(TYPEDEF)

Stores the nominal electrical and geometric parameters for a cable design, with attributes:

$(TYPEDFIELDS)
"""
struct NominalData
    "Cable designation as per DIN VDE 0271/0276."
    designation_code::Union{Nothing,String}
    "Rated phase-to-earth voltage \\[kV\\]."
    U0::Union{Nothing,Number}
    "Rated phase-to-phase voltage \\[kV\\]."
    "Cross-sectional area of the conductor \\[mm²\\]."
    U::Union{Nothing,Number}
    conductor_cross_section::Union{Nothing,Number}
    "Cross-sectional area of the screen \\[mm²\\]."
    screen_cross_section::Union{Nothing,Number}
    "Cross-sectional area of the armor \\[mm²\\]."
    armor_cross_section::Union{Nothing,Number}
    "Base (DC) resistance of the cable core \\[Ω/km\\]."
    resistance::Union{Nothing,Number}
    "Capacitance of the main insulation \\[μF/km\\]."
    capacitance::Union{Nothing,Number}
    "Inductance of the cable (trifoil formation) \\[mH/km\\]."
    inductance::Union{Nothing,Number}

    @doc """
    $(TYPEDSIGNATURES)

    Initializes a [`NominalData`](@ref) object with optional default values.

    # Arguments
    - `designation_code`: Cable designation (default: `nothing`).
    - `U0`: Phase-to-earth voltage rating \\[kV\\] (default: `nothing`).
    - `U`: Phase-to-phase voltage rating \\[kV\\] (default: `nothing`).
    - `conductor_cross_section`: Conductor cross-section \\[mm²\\] (default: `nothing`).
    - `screen_cross_section`: Screen cross-section \\[mm²\\] (default: `nothing`).
    - `armor_cross_section`: Armor cross-section \\[mm²\\] (default: `nothing`).
    - `resistance`: Cable resistance \\[Ω/km\\] (default: `nothing`).
    - `capacitance`: Cable capacitance \\[μF/km\\] (default: `nothing`).
    - `inductance`: Cable inductance (trifoil) \\[mH/km\\] (default: `nothing`).

    # Returns
    An instance of [`NominalData`](@ref) with the specified nominal properties.

    # Examples
    ```julia
    nominal_data = $(FUNCTIONNAME)(
    	conductor_cross_section=1000,
    	resistance=0.0291,
    	capacitance=0.39,
    )
    println(nominal_data.conductor_cross_section)
    println(nominal_data.resistance)
    println(nominal_data.capacitance)
    ```
    """
    function NominalData(;
        designation_code::Union{Nothing,String}=nothing,
        U0::Union{Nothing,Number}=nothing,
        U::Union{Nothing,Number}=nothing,
        conductor_cross_section::Union{Nothing,Number}=nothing,
        screen_cross_section::Union{Nothing,Number}=nothing,
        armor_cross_section::Union{Nothing,Number}=nothing,
        resistance::Union{Nothing,Number}=nothing,
        capacitance::Union{Nothing,Number}=nothing,
        inductance::Union{Nothing,Number}=nothing,
    )
        return new(
            designation_code,
            U0,
            U,
            conductor_cross_section,
            screen_cross_section,
            armor_cross_section,
            resistance,
            capacitance,
            inductance,
        )
    end
end

"""
$(TYPEDEF)

Represents a [`CableComponent`](@ref), i.e. a group of [`AbstractCablePart`](@ref) objects, with the equivalent geometric and material properties:

$(TYPEDFIELDS)

!!! info "Definition & application"
	Cable components operate as containers for multiple cable parts, allowing the calculation of effective electromagnetic (EM) properties (``\\sigma, \\varepsilon, \\mu``). This is performed by transforming the physical objects within the [`CableComponent`](@ref) into one equivalent coaxial homogeneous structure comprised of one conductor and one insulator, each one represented by effective [`Material`](@ref) types stored in `conductor_props` and `insulator_props` fields.

	The effective properties approach is widely adopted in EMT-type simulations, and involves locking the internal and external radii of the conductor and insulator parts, respectively, and calculating the equivalent EM properties in order to match the previously determined values of R, L, C and G [916943](@cite) [1458878](@cite).

	In applications, the [`CableComponent`](@ref) type is mapped to the main cable structures described in manufacturer datasheets, e.g., core, sheath, armor and jacket.
"""
mutable struct CableComponent
    "Cable component identification (e.g. core/sheath/armor)."
    id::String
    "The conductor group containing all conductive parts."
    conductor_group::ConductorGroup
    "Effective properties of the equivalent coaxial conductor."
    conductor_props::Material
    "The insulator group containing all insulating parts."
    insulator_group::InsulatorGroup
    "Effective properties of the equivalent coaxial insulator."
    insulator_props::Material

    @doc """
    $(TYPEDSIGNATURES)

    Initializes a [`CableComponent`](@ref) object based on its constituent conductor and insulator groups. The constructor performs the following sequence of steps:

    1.  Validate that the conductor and insulator groups have matching radii at their interface.
    2.  Obtain the lumped-parameter values (R, L, C, G) from the conductor and insulator groups, which are computed within their respective constructors.
    3.  Calculate the correction factors and equivalent electromagnetic properties of the conductor and insulator groups:


    | Quantity | Symbol | Function |
    |----------|--------|----------|
    | Resistivity (conductor) | ``\\rho_{con}`` | [`calc_equivalent_rho`](@ref) |
    | Permeability (conductor) | ``\\mu_{con}`` | [`calc_equivalent_mu`](@ref) |
    | Resistivity (insulator) | ``\\rho_{ins}`` | [`calc_sigma_lossfact`](@ref) |
    | Permittivity (insulation) | ``\\varepsilon_{ins}`` | [`calc_equivalent_eps`](@ref) |
    | Permeability (insulation) | ``\\mu_{ins}`` | [`calc_solenoid_correction`](@ref) |

    # Arguments

    - `id`: Cable component identification (e.g. core/sheath/armor).
    - `conductor_group`: The conductor group containing all conductive parts.
    - `insulator_group`: The insulator group containing all insulating parts.

    # Returns

    - A [`CableComponent`](@ref) instance with calculated equivalent properties.

    # Examples

    ```julia
    conductor_group = ConductorGroup(...)
    insulator_group = InsulatorGroup(...)
    cable = $(FUNCTIONNAME)("component_id", conductor_group, insulator_group)  # Create cable component with base parameters @ 50 Hz
    ```

    # See also

    - [`calc_equivalent_rho`](@ref)
    - [`calc_equivalent_mu`](@ref)
    - [`calc_equivalent_eps`](@ref)
    - [`calc_sigma_lossfact`](@ref)
    - [`calc_solenoid_correction`](@ref)
    """
    function CableComponent(
        id::String,
        conductor_group::ConductorGroup,
        insulator_group::InsulatorGroup)
        # Validate the geometry
        if conductor_group.radius_ext != insulator_group.radius_in
            error("Conductor outer radius must match insulator inner radius.")
        end

        # Use pre-calculated values from the groups
        radius_in_con = conductor_group.radius_in
        radius_ext_con = conductor_group.radius_ext
        radius_ext_ins = insulator_group.radius_ext

        # Equivalent conductor properties
        rho_con =
            calc_equivalent_rho(conductor_group.resistance, radius_ext_con, radius_in_con)
        mu_con = calc_equivalent_mu(conductor_group.gmr, radius_ext_con, radius_in_con)
        alpha_con = conductor_group.alpha
        conductor_props =
            Material(rho_con, 0.0, mu_con, conductor_group.layers[1].temperature, alpha_con)

        # Insulator properties
        C_eq = insulator_group.shunt_capacitance
        G_eq = insulator_group.shunt_conductance
        eps_ins = calc_equivalent_eps(C_eq, radius_ext_ins, radius_ext_con)
        rho_ins = 1 / calc_sigma_lossfact(G_eq, radius_ext_con, radius_ext_ins)
        correction_mu_ins = calc_solenoid_correction(
            conductor_group.num_turns,
            radius_ext_con,
            radius_ext_ins,
        )
        mu_ins = correction_mu_ins
        insulator_props =
            Material(rho_ins, eps_ins, mu_ins, insulator_group.layers[1].temperature, 0.0)

        # Initialize object
        return new(
            id,
            conductor_group,
            conductor_props,
            insulator_group,
            insulator_props,
        )
    end
end

"""
$(TYPEDEF)

Represents the design of a cable, including its unique identifier, nominal data, and components.

$(TYPEDFIELDS)
"""
mutable struct CableDesign
    "Unique identifier for the cable design."
    cable_id::String
    "Informative reference data."
    nominal_data::NominalData
    "Vector of cable components."
    components::Vector{CableComponent}

    @doc """
    $(TYPEDSIGNATURES)

    Constructs a [`CableDesign`](@ref) instance.

    # Arguments

    - `cable_id`: Unique identifier for the cable design.
    - `component`: Initial [`CableComponent`](@ref) for the design.
    - `nominal_data`: Reference data for the cable design. Default: `NominalData()`.

    # Returns

    - A [`CableDesign`](@ref) object with the specified properties.

    # Examples

    ```julia
    conductor_group = ConductorGroup(central_conductor)
    insulator_group = InsulatorGroup(main_insulator)
    component = CableComponent(conductor_group, insulator_group)
    design = $(FUNCTIONNAME)("example", component)
    ```

    # See also

    - [`CableComponent`](@ref)
    - [`ConductorGroup`](@ref)
    - [`InsulatorGroup`](@ref)
    """
    function CableDesign(
        cable_id::String,
        component::CableComponent;
        nominal_data::NominalData=NominalData(),
    )
        return new(cable_id, nominal_data, [component])
    end
end

"""
$(TYPEDSIGNATURES)

Constructs a [`CableDesign`](@ref) instance from conductor and insulator groups.

# Arguments

- `cable_id`: Unique identifier for the cable design.
- `conductor_group`: The [`ConductorGroup`](@ref) for the component.
- `insulator_group`: The [`InsulatorGroup`](@ref) for the component.
- `component_id`: ID for the cable component. Default: "core".
- `nominal_data`: Reference data for the cable design. Default: `NominalData()`.

# Returns

- A [`CableDesign`](@ref) object with the specified properties.

# Examples

```julia
conductor_group = ConductorGroup(central_conductor)
insulator_group = InsulatorGroup(main_insulator)
design = $(FUNCTIONNAME)("example", conductor_group, insulator_group)
```
"""
function CableDesign(
    cable_id::String,
    conductor_group::ConductorGroup,
    insulator_group::InsulatorGroup;
    component_id::String="component1",
    nominal_data::NominalData=NominalData(),
)
    component = CableComponent(component_id, conductor_group, insulator_group)
    return CableDesign(cable_id, nominal_data, [component])
end


"""
$(TYPEDSIGNATURES)

Adds a cable component to an existing [`CableDesign`](@ref).

# Arguments

- `design`: A [`CableDesign`](@ref) object where the component will be added.
- `component`: A [`CableComponent`](@ref) to add to the design.

# Returns

- The modified [`CableDesign`](@ref) object.

# Notes

If a component with the same ID already exists, it will be overwritten, and a warning will be logged.

# Examples

```julia
conductor_group = ConductorGroup(wire_array)
insulator_group = InsulatorGroup(insulation)
component = CableComponent("sheath", conductor_group, insulator_group)
$(FUNCTIONNAME)(design, component)
```

# See also

- [`CableDesign`](@ref)
- [`CableComponent`](@ref)
"""
function add!(
    design::CableDesign,
    component::CableComponent,
)
    # Check for existing component with same ID
    existing_idx = findfirst(comp -> comp.id == component.id, design.components)

    if !isnothing(existing_idx)
        @warn "Component with ID '$(component.id)' already exists in the CableDesign and will be overwritten."
        design.components[existing_idx] = component
    else
        # Add new component to the vector
        push!(design.components, component)
    end

    return design
end

"""
$(TYPEDSIGNATURES)

Adds a cable component to an existing [`CableDesign`](@ref) using separate conductor and insulator groups. Performs as a convenience wrapper to construct the [`CableComponent`](@ref) object with reduced boilerplate.

# Arguments

- `design`: A [`CableDesign`](@ref) object where the component will be added.
- `component_id`: ID for the new component.
- `conductor_group`: A [`ConductorGroup`](@ref) for the component.
- `insulator_group`: An [`InsulatorGroup`](@ref) for the component.

# Returns

- The modified [`CableDesign`](@ref) object.

# Examples

```julia
conductor_group = ConductorGroup(wire_array)
insulator_group = InsulatorGroup(insulation)
$(FUNCTIONNAME)(design, "shield", conductor_group, insulator_group)
```

# See also

- [`CableDesign`](@ref)
- [`CableComponent`](@ref)
"""
function add!(
    design::CableDesign,
    component_id::String,
    conductor_group::ConductorGroup,
    insulator_group::InsulatorGroup,
)
    # Create new component
    component = CableComponent(component_id, conductor_group, insulator_group)

    # Call the main function
    return add!(design, component)
end