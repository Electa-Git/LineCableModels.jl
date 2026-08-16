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

mutable struct ConductorGroup{T <: REALSCALAR} <: AbstractConductorPart{T}
    "Inner radius of the conductor group \\[m\\]."
    r_in::T
    "Outer radius of the conductor group \\[m\\]."
    r_ex::T
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
            r_in::T,
            r_ex::T,
            cross_section::T,
            num_wires::Int,
            num_turns::T,
            resistance::T,
            alpha::T,
            gmr::T,
            layers::Vector{AbstractConductorPart{T}}
    ) where {T}
        return new{T}(r_in, r_ex, cross_section, num_wires, num_turns,
            resistance, alpha, gmr, layers)
    end

    function ConductorGroup{T}(central::AbstractConductorPart{T}) where {T}
        num_wires::Int = 0
        num_turns::T = zero(T)

        # only touch fields that exist inside the guarded branches
        if central isa CircStrands{T}
            num_wires = central.num_wires
            num_turns = central.pitch_length > zero(T) ? one(T) / central.pitch_length :
                        zero(T)
        elseif central isa Strip{T}
            num_wires = 1
            num_turns = central.pitch_length > zero(T) ? one(T) / central.pitch_length :
                        zero(T)
        end

        return new{T}(
            central.r_in,
            central.r_ex,
            central.cross_section,
            num_wires,
            num_turns,
            central.resistance,
            central.material_props.alpha,
            central.gmr,
            AbstractConductorPart{T}[central]
        )
    end
end

# Outer helper that infers T from the central part
ConductorGroup(con::AbstractConductorPart{T}) where {T} = ConductorGroup{T}(con)

"""
$(TYPEDSIGNATURES)

Add a new conductor part to a [`ConductorGroup`](@ref), validate its numeric
inputs, and promote the group’s numeric type when required.

This updates `gmr`, `resistance`, `alpha`, `r_ex`, `cross_section`, and
`num_wires`. The new part temperature defaults to the temperature of the
first layer, and its inner radius defaults to the external radius of the
existing group.

# Behavior

1. Apply part-level keyword defaults.
2. Default `r_in` to `group.r_ex` if absent.
3. Compute `Tnew = resolve_T(group, r_in, args..., values(kwargs)...)`.
4. If `Tnew === T`, mutate in place; else `coerce_to_T(group, Tnew)` then mutate and **return the promoted group**.

# Arguments

- `group`: [`ConductorGroup`](@ref) object to which the new part will be added.
- `part_type`: Type of conductor part to add ([`AbstractConductorPart`](@ref)).
- `args...`: Positional arguments specific to the constructor of the `part_type` ([`AbstractConductorPart`](@ref)) \\[various\\].
- `kwargs...`: Named arguments for the constructor including optional values specific to the constructor of the `part_type` ([`AbstractConductorPart`](@ref)) \\[various\\].

# Returns

- The updated [`ConductorGroup`](@ref). This is `group` when its numeric type
  is unchanged, or a promoted group when the new part requires another numeric
  type.

!!! warning "Note"
    - The current outer radius is supplied automatically as `r_in`. For radial
      parts, pass exactly one of `radius` or `thickness`; do not pass a layer
      object as a constructor input.
    - For uncertain geometries, the existing outer-radius derivative graph is
      retained. Adjacent layers therefore share one physical boundary and
      cumulative-radius covariance is preserved across part types.

# Examples

```jldoctest
using LineCableModels.DataModel: CircStrands, ConductorGroup
using LineCableModels.Materials: Material

material = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
conductor = ConductorGroup(CircStrands(0.0, 0.001, 1, 0.0, material))
conductor = $(FUNCTIONNAME)(conductor, CircStrands, 0.001, 6, 15.0, material)
@assert length(conductor.layers) == 2
# output
```

"""
function add!(
        group::ConductorGroup{T},
        part_type::Type{C},
        args...;
        kwargs...
) where {T, C <: AbstractConductorPart}

    # 1) Merge declared keyword defaults for this part type
    kwv = _with_kwdefaults(C, (; kwargs...))

    # 2) Default stacking: inner radius = current outer radius unless overridden
    rin = get(kwv, :r_in, group.r_ex)
    kwv = haskey(kwv, :r_in) ? kwv : merge(kwv, (; r_in = rin))

    part_args, kwv = _resolve_group_radial_args(C, rin, args, kwv)

    # 3) Decide target numeric type using *current group + raw inputs*
    Tnew = resolve_T(group, rin, part_args..., values(kwv)...)

    if Tnew === T
        # 4a) Fast path: mutate in place
        return _do_add!(group, C, part_args...; kwv...)
    else
        @warn """
          Adding a `$Tnew` part to a `ConductorGroup{$T}` returns a **promoted** group.
          Capture the result:  group = add!(group, $C, …)
          """
        promoted = coerce_to_T(group, Tnew)
        return _do_add!(promoted, C, part_args...; kwv...)
    end
end

"""
$(TYPEDSIGNATURES)

Internal, in-place insertion (no promotion logic). Assumes `:r_in` was materialized.
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

    # Validate and parse with the part's own numeric input sequence.
    ntv = Validation.validate!(C, kw.r_in, args...; kw...)

    # Coerce validated values to group’s T and call strict numeric core
    order = (Validation.required_fields(C)..., Validation.keyword_fields(C)...)
    coerced = _coerced_args(C, ntv, Tg, order)      # respects coercive_fields(C)
    new_part = C(coerced...)

    return _append_conductor!(group, new_part)
end

function _append_conductor!(
    group::ConductorGroup{T},
    new_part::AbstractConductorPart{T},
) where {T}
    isapprox(new_part.r_in, group.r_ex) || throw(ArgumentError(
        "new conductor layer must start at the group's current outer radius " *
        "$(group.r_ex); got $(new_part.r_in)",
    ))

    # Update equivalent properties
    group.gmr = calc_equivalent_gmr(group, new_part)
    group.alpha = calc_equivalent_alpha(group.alpha, group.resistance,
        new_part.material_props.alpha,
        new_part.resistance)
    group.resistance = calc_parallel_equivalent(group.resistance, new_part.resistance)
    group.r_ex += (new_part.r_ex - new_part.r_in)
    group.cross_section += new_part.cross_section

    # CircStrands / Strip bookkeeping
    if new_part isa CircStrands || new_part isa Strip
        old_wires = group.num_wires
        old_turns = group.num_turns
        nw = new_part isa CircStrands ? new_part.num_wires : 1
        nt = new_part.pitch_length > 0 ? inv(new_part.pitch_length) : zero(Tg)
        group.num_wires += nw
        group.num_turns = (old_wires * old_turns + nw * nt) / group.num_wires
    end

    push!(group.layers, new_part)
    return group
end

function add!(
    group::ConductorGroup{T},
    new_part::AbstractConductorPart{U},
) where {T,U}
    target_type = promote_type(T, U)
    if target_type === T
        return _append_conductor!(group, coerce_to_T(new_part, T))
    end
    promoted = coerce_to_T(group, target_type)
    return _append_conductor!(promoted, coerce_to_T(new_part, target_type))
end

include("conductorgroup/base.jl")
