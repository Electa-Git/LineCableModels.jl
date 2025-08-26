import ..LineCableModels: add!

"""
$(TYPEDEF)

Represents a composite coaxial insulator group assembled from multiple insulating layers.

This structure serves as a container for different [`AbstractInsulatorPart`](@ref) elements
(such as insulators and semiconductors) arranged in concentric layers.
The `InsulatorGroup` aggregates these individual parts and provides equivalent electrical
properties that represent the composite behavior of the entire assembly, stored in the attributes:

$(TYPEDFIELDS)
"""
mutable struct InsulatorGroup <: AbstractInsulatorPart
    "Inner radius of the insulator group \\[m\\]."
    radius_in::Number
    "Outer radius of the insulator group \\[m\\]."
    radius_ext::Number
    "Cross-sectional area of the entire insulator group \\[m²\\]."
    cross_section::Number
    "Shunt capacitance per unit length of the insulator group \\[F/m\\]."
    shunt_capacitance::Number
    "Shunt conductance per unit length of the insulator group \\[S·m\\]."
    shunt_conductance::Number
    "Vector of insulator layer components."
    layers::Vector{AbstractInsulatorPart}

    @doc """
    $(TYPEDSIGNATURES)

    Constructs an [`InsulatorGroup`](@ref) instance initializing with the initial insulator part.

    # Arguments

    - `initial_insulator`: An [`AbstractInsulatorPart`](@ref) object located at the innermost position of the insulator group.

    # Returns

    - An [`InsulatorGroup`](@ref) object initialized with geometric and electrical properties derived from the initial insulator.

    # Examples

    ```julia
    material_props = Material(1e10, 3.0, 1.0, 20.0, 0.0)
    initial_insulator = Insulator(0.01, 0.015, material_props)
    insulator_group = $(FUNCTIONNAME)(initial_insulator)
    println(insulator_group.layers)           # Output: [initial_insulator]
    println(insulator_group.shunt_capacitance) # Output: Capacitance in [F/m]
    ```
    """
    function InsulatorGroup(initial_insulator::AbstractInsulatorPart)
        # Initialize object
        return new(
            initial_insulator.radius_in,
            initial_insulator.radius_ext,
            initial_insulator.cross_section,
            initial_insulator.shunt_capacitance,
            initial_insulator.shunt_conductance,
            [initial_insulator],
        )
    end
end

"""
$(TYPEDSIGNATURES)

Adds a new part to an existing [`InsulatorGroup`](@ref) object and updates its equivalent electrical parameters.

# Arguments

- `group`: [`InsulatorGroup`](@ref) object to which the new part will be added.
- `part_type`: Type of insulator part to add ([`AbstractInsulatorPart`](@ref)).
- `args...`: Positional arguments specific to the constructor of the `part_type` ([`AbstractInsulatorPart`](@ref)) \\[various\\].
- `kwargs...`: Named arguments for the constructor including optional values specific to the constructor of the `part_type` ([`AbstractInsulatorPart`](@ref)) \\[various\\].

# Returns

- The function modifies the [`InsulatorGroup`](@ref) instance in place and does not return a value.

# Notes

- Updates `shunt_capacitance`, `shunt_conductance`, `radius_ext`, and `cross_section` to account for the new part.
- The `radius_in` of the new part defaults to the external radius of the existing insulator group if not specified.

!!! warning "Note"
	- When an [`AbstractCablePart`](@ref) is provided as `radius_in`, the constructor retrieves its `radius_ext` value, allowing the new cable part to be placed directly over the existing part in a layered cable design.
	- In case of uncertain measurements, if the added cable part is of a different type than the existing one, the uncertainty is removed from the radius value before being passed to the new component. This ensures that measurement uncertainties do not inappropriately cascade across different cable parts.

# Examples

```julia
material_props = Material(1e10, 3.0, 1.0, 20.0, 0.0)
insulator_group = InsulatorGroup(Insulator(0.01, 0.015, material_props))
$(FUNCTIONNAME)(insulator_group, Semicon, 0.015, 0.018, material_props)
```

# See also

- [`InsulatorGroup`](@ref)
- [`Insulator`](@ref)
- [`Semicon`](@ref)
- [`calc_parallel_equivalent`](@ref)
"""
function add!(
    group::InsulatorGroup,
    part_type::Type{T},  # The type of insulator part (Insulator, Semicon)
    args...;  # Arguments specific to the part type
    f::Number=f₀,  # Add the f parameter with default value f₀
    kwargs...,
) where {T<:AbstractInsulatorPart}
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


    # For admittances (parallel combination)
    ω = 2 * π * f
    Y_group = Complex(group.shunt_conductance, ω * group.shunt_capacitance)
    Y_newpart = Complex(new_part.shunt_conductance, ω * new_part.shunt_capacitance)
    Y_equiv = calc_parallel_equivalent(Y_group, Y_newpart)
    group.shunt_capacitance = imag(Y_equiv) / ω
    group.shunt_conductance = real(Y_equiv)

    # Update geometric properties
    group.radius_ext += (new_part.radius_ext - new_part.radius_in)
    group.cross_section += new_part.cross_section

    # Add to layers
    push!(group.layers, new_part)
    group
end



