import ..LineCableModels: add!

"""
$(TYPEDEF)

Represents a composite conductor group assembled from multiple conductive layers or stranded wires.

This structure serves as a container for different [`AbstractConductorPart`](@ref) elements 
(such as wire arrays, strips, and tubular conductors) arranged in concentric layers. 
The `ConductorGroup` aggregates these individual parts and provides equivalent electrical 
properties that represent the composite behavior of the entire assembly.

# Attributes

$(TYPEDFIELDS)
"""
mutable struct ConductorGroup <: AbstractConductorPart
    "Inner radius of the conductor group \\[m\\]."
    radius_in::Number
    "Outer radius of the conductor group \\[m\\]."
    radius_ext::Number
    "Cross-sectional area of the entire conductor group \\[m²\\]."
    cross_section::Number
    "Number of individual wires in the conductor group \\[dimensionless\\]."
    num_wires::Number
    "Number of turns per meter of each wire strand \\[1/m\\]."
    num_turns::Number
    "DC resistance of the conductor group \\[Ω\\]."
    resistance::Number
    "Temperature coefficient of resistance \\[1/°C\\]."
    alpha::Number
    "Geometric mean radius of the conductor group \\[m\\]."
    gmr::Number
    "Vector of conductor layer components."
    layers::Vector{AbstractConductorPart}

    @doc """
    $(TYPEDSIGNATURES)

    Constructs a [`ConductorGroup`](@ref) instance initializing with the central conductor part.

    # Arguments

    - `central_conductor`: An [`AbstractConductorPart`](@ref) object located at the center of the conductor group.

    # Returns

    - A [`ConductorGroup`](@ref) object initialized with geometric and electrical properties derived from the central conductor.

    # Examples

    ```julia
    material_props = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
    central_strip = Strip(0.01, 0.002, 0.05, 10, material_props)
    conductor_group = $(FUNCTIONNAME)(central_strip)
    println(conductor_group.layers)      # Output: [central_strip]
    println(conductor_group.resistance)  # Output: Resistance in \\[Ω\\]
    ```
    """
    function ConductorGroup(central_conductor::AbstractConductorPart)
        num_wires = 0
        num_turns = 0.0

        if central_conductor isa WireArray || central_conductor isa Strip
            num_wires = central_conductor isa WireArray ? central_conductor.num_wires : 1
            num_turns =
                central_conductor.pitch_length > 0 ?
                1 / central_conductor.pitch_length : 0.0
        end

        # Initialize object
        return new(
            central_conductor.radius_in,
            central_conductor.radius_ext,
            central_conductor.cross_section,
            num_wires,
            num_turns,
            central_conductor.resistance,
            central_conductor.material_props.alpha,
            central_conductor.gmr,
            [central_conductor],
        )
    end
end

"""
$(TYPEDSIGNATURES)

Adds a new part to an existing [`ConductorGroup`](@ref) object and updates its equivalent electrical parameters.

# Arguments

- `group`: [`ConductorGroup`](@ref) object to which the new part will be added.
- `part_type`: Type of conductor part to add ([`AbstractConductorPart`](@ref)).
- `args...`: Positional arguments specific to the constructor of the `part_type` ([`AbstractConductorPart`](@ref)) \\[various\\].
- `kwargs...`: Named arguments for the constructor including optional values specific to the constructor of the `part_type` ([`AbstractConductorPart`](@ref)) \\[various\\].

# Returns

- The function modifies the [`ConductorGroup`](@ref) instance in place and does not return a value.

# Notes

- Updates `gmr`, `resistance`, `alpha`, `radius_ext`, `cross_section`, and `num_wires` to account for the new part.
- The `temperature` of the new part defaults to the temperature of the first layer if not specified.
- The `radius_in` of the new part defaults to the external radius of the existing conductor if not specified.

!!! warning "Note"
	- When an [`AbstractCablePart`](@ref) is provided as `radius_in`, the constructor retrieves its `radius_ext` value, allowing the new cable part to be placed directly over the existing part in a layered cable design.
	- In case of uncertain measurements, if the added cable part is of a different type than the existing one, the uncertainty is removed from the radius value before being passed to the new component. This ensures that measurement uncertainties do not inappropriately cascade across different cable parts.

# Examples

```julia
material_props = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
conductor = ConductorGroup(Strip(0.01, 0.002, 0.05, 10, material_props))
$(FUNCTIONNAME)(conductor, WireArray, 0.02, 0.002, 7, 15, material_props, temperature = 25)
```

# See also

- [`ConductorGroup`](@ref)
- [`WireArray`](@ref)
- [`Strip`](@ref)
- [`Tubular`](@ref)
- [`calc_equivalent_gmr`](@ref)
- [`calc_parallel_equivalent`](@ref)
- [`calc_equivalent_alpha`](@ref)
"""
function add!(
    group::ConductorGroup,
    part_type::Type{T},  # The type of conductor part (WireArray, Strip, Tubular)
    args...;  # Arguments specific to the part type
    kwargs...,
) where {T<:AbstractConductorPart}
    # Infer default properties
    radius_in = get(kwargs, :radius_in, group.radius_ext)

    # Create a new named tuple with default temperature
    default_kwargs = (temperature=group.layers[1].temperature,)

    # Create a merged kwargs dictionary
    merged_kwargs = Dict{Symbol,Any}()

    # Add defaults first
    for (key, value) in pairs(default_kwargs)
        merged_kwargs[key] = value
    end

    # Then override with user-provided values
    for (key, value) in kwargs
        merged_kwargs[key] = value
    end

    # Create the new part
    new_part = T(radius_in, args...; merged_kwargs...)

    # Update the Conductor with the new part
    group.gmr = calc_equivalent_gmr(group, new_part)
    group.alpha = calc_equivalent_alpha(
        group.alpha,
        group.resistance,
        new_part.material_props.alpha,
        new_part.resistance,
    )

    group.resistance = calc_parallel_equivalent(group.resistance, new_part.resistance)
    group.radius_ext += (new_part.radius_ext - new_part.radius_in)
    group.cross_section += new_part.cross_section

    # For WireArray, update the number of wires and turns
    if new_part isa WireArray || new_part isa Strip
        cum_num_wires = group.num_wires
        cum_num_turns = group.num_turns
        new_wires = new_part isa WireArray ? new_part.num_wires : 1
        new_turns = new_part.pitch_length > 0 ? 1 / new_part.pitch_length : 0.0
        group.num_wires += new_wires
        group.num_turns =
            (cum_num_wires * cum_num_turns + new_wires * new_turns) / (group.num_wires)
    end

    push!(group.layers, new_part)
    group
end
