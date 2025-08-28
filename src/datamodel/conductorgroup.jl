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
mutable struct ConductorGroup{T<:REALSCALAR} <: AbstractConductorPart{T}
    "Inner radius of the conductor group \\[m\\]."
    radius_in::T
    "Outer radius of the conductor group \\[m\\]."
    radius_ext::T
    "Cross-sectional area of the entire conductor group \\[m²\\]."
    cross_section::T
    "Number of individual wires in the conductor group \\[dimensionless\\]."
    num_wires::Int
    "Number of turns per meter of each wire strand \\[1/m\\]."
    num_turns::T
    "DC resistance of the conductor group \\[Ω\\]."
    resistance::T
    "Temperature coefficient of resistance \\[1/°C\\]."
    alpha::T
    "Geometric mean radius of the conductor group \\[m\\]."
    gmr::T
    "Vector of conductor layer components."
    layers::Vector{AbstractConductorPart{T}}

    @doc """
    $(TYPEDSIGNATURES)

    Constructs a [`ConductorGroup`](@ref) instance initializing with the central conductor part.

    # Arguments

    - `central_conductor`: An [`AbstractConductorPart`](@ref) object located at the center of the conductor group.

    # Returns

    - A [`ConductorGroup`](@ref) object initialized with geometric and electrical properties derived from the central conductor.
    """
    function ConductorGroup{T}(
        radius_in::T,
        radius_ext::T,
        cross_section::T,
        num_wires::Int,
        num_turns::T,
        resistance::T,
        alpha::T,
        gmr::T,
        layers::Vector{AbstractConductorPart{T}},
    ) where {T}
        return new{T}(radius_in, radius_ext, cross_section, num_wires, num_turns,
            resistance, alpha, gmr, layers)
    end

    function ConductorGroup{T}(central::AbstractConductorPart{T}) where {T}
        num_wires::Int = 0
        num_turns::T = zero(T)

        # only touch fields that exist inside the guarded branches
        if central isa WireArray{T}
            num_wires = central.num_wires
            num_turns = central.pitch_length > zero(T) ? one(T) / central.pitch_length : zero(T)
        elseif central isa Strip{T}
            num_wires = 1
            num_turns = central.pitch_length > zero(T) ? one(T) / central.pitch_length : zero(T)
        end

        return new{T}(
            central.radius_in,
            central.radius_ext,
            central.cross_section,
            num_wires,
            num_turns,
            central.resistance,
            central.material_props.alpha,
            central.gmr,
            AbstractConductorPart{T}[central],
        )
    end
end

# Outer helper that infers T from the central part
ConductorGroup(con::AbstractConductorPart{T}) where {T} = ConductorGroup{T}(con)

"""
$(TYPEDSIGNATURES)

Add a new conductor part to a [`ConductorGroup`](@ref), validating raw inputs,
normalizing proxies, and **promoting** the group’s numeric type if required.

# Behavior:

1. Apply part-level keyword defaults.
2. Default `radius_in` to `group.radius_ext` if absent.
3. Compute `Tnew = resolve_T(group, radius_in, args..., values(kwargs)...)`.
4. If `Tnew === T`, mutate in place; else `coerce_to_T(group, Tnew)` then mutate and **return the promoted group**.

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
    group::ConductorGroup{T},
    part_type::Type{C},
    args...;
    kwargs...
) where {T,C<:AbstractConductorPart}

    # 1) Merge declared keyword defaults for this part type
    kwv = _with_kwdefaults(C, (; kwargs...))

    # 2) Default stacking: inner radius = current outer radius unless overridden
    rin = get(kwv, :radius_in, group.radius_ext)
    kwv = haskey(kwv, :radius_in) ? kwv : merge(kwv, (; radius_in=rin))

    # 3) Decide target numeric type using *current group + raw inputs*
    Tnew = resolve_T(group, rin, args..., values(kwv)...)

    if Tnew === T
        # 4a) Fast path: mutate in place
        return _do_add!(group, C, args...; kwv...)
    else
        @warn """
        Adding a `$Tnew` part to a `ConductorGroup{$T}` returns a **promoted** group.
        Capture the result:  group = add!(group, $C, …)
        """
        promoted = coerce_to_T(group, Tnew)
        return _do_add!(promoted, C, args...; kwv...)
    end
end

"""
$(TYPEDSIGNATURES)

Internal, in-place insertion (no promotion logic). Assumes `:radius_in` was materialized.
Runs Validation → parsing, then coerces fields to the group’s `T` and updates
equivalent properties and book-keeping.
"""
function _do_add!(
    group::ConductorGroup{Tg},
    C::Type{<:AbstractConductorPart},
    args...;
    kwargs...
) where {Tg}
    # Materialize keyword args into a NamedTuple (never poke Base.Pairs internals)
    kw = (; kwargs...)

    # Validate + parse with the part’s own pipeline (proxies resolved here)
    ntv = Validation.validate!(C, kw.radius_in, args...; kw...)

    # Coerce validated values to group’s T and call strict numeric core
    order = (Validation.required_fields(C)..., Validation.keyword_fields(C)...)
    coerced = _coerced_args(C, ntv, Tg, order)      # respects coercive_fields(C)
    new_part = C(coerced...)

    # Update equivalent properties
    group.gmr = calc_equivalent_gmr(group, new_part)
    group.alpha = calc_equivalent_alpha(group.alpha, group.resistance,
        new_part.material_props.alpha,
        new_part.resistance)
    group.resistance = calc_parallel_equivalent(group.resistance, new_part.resistance)
    group.radius_ext += (new_part.radius_ext - new_part.radius_in)
    group.cross_section += new_part.cross_section

    # WireArray / Strip bookkeeping
    if new_part isa WireArray || new_part isa Strip
        old_wires = group.num_wires
        old_turns = group.num_turns
        nw = new_part isa WireArray ? new_part.num_wires : 1
        nt = new_part.pitch_length > 0 ? inv(new_part.pitch_length) : zero(Tg)
        group.num_wires += nw
        group.num_turns = (old_wires * old_turns + nw * nt) / group.num_wires
    end

    push!(group.layers, new_part)
    return group
end

include("conductorgroup/base.jl")