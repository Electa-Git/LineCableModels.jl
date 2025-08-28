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
mutable struct InsulatorGroup{T<:REALSCALAR} <: AbstractInsulatorPart{T}
    "Inner radius of the insulator group \\[m\\]."
    radius_in::T
    "Outer radius of the insulator group \\[m\\]."
    radius_ext::T
    "Cross-sectional area of the entire insulator group \\[m²\\]."
    cross_section::T
    "Shunt capacitance per unit length of the insulator group \\[F/m\\]."
    shunt_capacitance::T
    "Shunt conductance per unit length of the insulator group \\[S·m\\]."
    shunt_conductance::T
    "Vector of insulator layer components."
    layers::Vector{AbstractInsulatorPart{T}}

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
    function InsulatorGroup{T}(
        radius_in::T,
        radius_ext::T,
        cross_section::T,
        shunt_capacitance::T,
        shunt_conductance::T,
        layers::Vector{AbstractInsulatorPart{T}},
    ) where {T}
        return new{T}(radius_in, radius_ext, cross_section,
            shunt_capacitance, shunt_conductance, layers)
    end

    function InsulatorGroup{T}(initial_insulator::AbstractInsulatorPart{T}) where {T}
        return new{T}(
            initial_insulator.radius_in,
            initial_insulator.radius_ext,
            initial_insulator.cross_section,
            initial_insulator.shunt_capacitance,
            initial_insulator.shunt_conductance,
            AbstractInsulatorPart{T}[initial_insulator],
        )
    end
end

# Convenience outer
InsulatorGroup(ins::AbstractInsulatorPart{T}) where {T} = InsulatorGroup{T}(ins)

"""
$(TYPEDSIGNATURES)

Adds a new part to an existing [`InsulatorGroup`](@ref) object and updates its equivalent electrical parameters.

# Behavior:

1. Apply part-level keyword defaults (from `Validation.keyword_defaults`).
2. Default `radius_in` to `group.radius_ext` if absent.
3. Compute `Tnew = resolve_T(group, radius_in, args..., values(kwargs)..., f)`.
4. If `Tnew === T`, mutate in place; else `coerce_to_T(group, Tnew)` then mutate and **return the promoted group**.

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
    group::InsulatorGroup{T},
    part_type::Type{C},
    args...;
    f::Number=f₀,
    kwargs...
) where {T,C<:AbstractInsulatorPart}

    # 1) Merge declared keyword defaults for this part type
    kwv = _with_kwdefaults(C, (; kwargs...))

    # 2) Default stacking: inner radius = current outer radius unless overridden
    rin = get(kwv, :radius_in, group.radius_ext)
    kwv = haskey(kwv, :radius_in) ? kwv : merge(kwv, (; radius_in=rin))

    # 3) Decide target numeric type using *current group + raw inputs + f*
    Tnew = resolve_T(group, rin, args..., values(kwv)..., f)

    if Tnew === T
        # 4a) Fast path: mutate in place
        return _do_add!(group, C, args...; f, kwv...)
    else
        @warn """
        Adding a `$Tnew` part to an `InsulatorGroup{$T}` returns a **promoted** group.
        Capture the result:  group = add!(group, $C, …)
        """
        promoted = coerce_to_T(group, Tnew)
        return _do_add!(promoted, C, args...; f, kwv...)
    end
end

"""
$(TYPEDSIGNATURES)

Do the actual insertion for `InsulatorGroup` with the group already at the
correct scalar type. Validates/parses the part, coerces to the group’s `T`,
constructs the strict numeric core, and updates geometry and admittances at the
provided frequency.

Returns the mutated group (same object).
"""
function _do_add!(
    group::InsulatorGroup{Tg},
    C::Type{<:AbstractInsulatorPart},
    args...;
    f::Number=f₀,
    kwargs...
) where {Tg}

    # Materialize keyword args into a NamedTuple
    kw = (; kwargs...)

    # Validate + parse with the part’s own pipeline (proxies resolved here)
    ntv = Validation.validate!(C, kw.radius_in, args...; kw...)

    # Build argument order and coerce validated values to group’s T
    order = (Validation.required_fields(C)..., Validation.keyword_fields(C)...)
    coerced = _coerced_args(C, ntv, Tg, order)   # respects coercive_fields(C)
    new_part = C(coerced...)                      # call strict numeric core

    # Parallel admittances at frequency f
    ω = Tg(2π) * coerce_to_T(f, Tg)
    Yg = Complex(group.shunt_conductance, ω * group.shunt_capacitance)
    Yp = Complex(new_part.shunt_conductance, ω * new_part.shunt_capacitance)
    Ye = calc_parallel_equivalent(Yg, Yp)
    group.shunt_conductance = real(Ye)
    group.shunt_capacitance = imag(Ye) / ω

    # Update geometry
    group.radius_ext += new_part.radius_ext - new_part.radius_in
    group.cross_section += new_part.cross_section

    push!(group.layers, new_part)
    return group
end

include("insulatorgroup/base.jl")

