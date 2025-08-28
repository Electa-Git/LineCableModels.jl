import ..LineCableModels: add!

"""
$(TYPEDEF)

Represents a physically defined cable with position and phase mapping within a system.

$(TYPEDFIELDS)
"""
struct CablePosition{T<:REALSCALAR}
    "The [`CableDesign`](@ref) object assigned to this cable position."
    design_data::CableDesign{T}
    "Horizontal coordinate \\[m\\]."
    horz::T
    "Vertical coordinate \\[m\\]."
    vert::T
    "Phase mapping vector (aligned with design_data.components)."
    conn::Vector{Int}

    @doc """
    $(TYPEDSIGNATURES)

    Constructs a [`CablePosition`](@ref) instance with specified cable design, coordinates, and phase mapping.

    # Arguments

    - `cable`: A [`CableDesign`](@ref) object defining the cable structure.
    - `horz`: Horizontal coordinate \\[m\\].
    - `vert`: Vertical coordinate \\[m\\].
    - `conn`: A dictionary mapping component names to phase indices, or `nothing` for default mapping.

    # Returns

    - A [`CablePosition`](@ref) object with the assigned cable design, coordinates, and phase mapping.

    #=
    !!! note "Phase mapping"
         The `conn` argument is a `Dict` that maps the cable components to their respective phases. The values (1, 2, 3) represent the phase numbers (A, B, C) in a three-phase system. Components mapped to phase 0 will be Kron-eliminated (grounded). Components set to the same phase will be bundled into an equivalent phase.
    =#

    # Examples

    ```julia
    cable_design = CableDesign("example", nominal_data, components_dict)
    xa, ya = 0.0, -1.0  # Coordinates in meters

    # With explicit phase mapping
    cablepos1 = $(FUNCTIONNAME)(cable_design, xa, ya, Dict("core" => 1))

    # With default phase mapping (first component to phase 1, others to 0)
    default_cablepos = $(FUNCTIONNAME)(cable_design, xa, ya)
    ```

    # See also

    - [`CableDesign`](@ref)
    """
    function CablePosition{T}(
        cable::CableDesign{T},
        horz::T,
        vert::T,
        conn::Vector{Int},
    ) where {T<:REALSCALAR}
        # Validate: cable not empty
        @assert !isempty(cable.components) "CableDesign must contain at least one component"

        # Find outermost radius (last component)
        last_comp = cable.components[end]
        r_cond = last_comp.conductor_group.radius_ext
        r_ins = last_comp.insulator_group.radius_ext
        r_max = max(r_cond, r_ins)

        # Validate vertical position
        if iszero(vert)
            throw(ArgumentError("Vertical position cannot be exactly at the air/earth interface (z=0)"))
        end
        if abs(vert) < r_max
            throw(ArgumentError("Vertical position |$vert| must be ≥ cable's outer radius $r_max to avoid crossing z=0"))
        end

        return new{T}(cable, horz, vert, conn)
    end
end

"""
$(TYPEDSIGNATURES)

**Weakly-typed constructor** that infers `T` from the `cable` and coordinates, builds/validates the phase mapping, coerces inputs to `T`, and calls the typed kernel.
"""
function CablePosition(
    cable::Union{CableDesign,Nothing},
    horz::Number,
    vert::Number,
    conn::Union{Dict{String,Int},Nothing}=nothing,
)
    @assert !isnothing(cable) "A valid CableDesign must be provided"
    @assert !isempty(cable.components) "CableDesign must contain at least one component"

    # Build phase mapping vector aligned to component order
    names = [comp.id for comp in cable.components]
    conn_vector = if isnothing(conn)
        [i == 1 ? 1 : 0 for i in 1:length(names)]   # default: first component → phase 1, others grounded
    else
        [get(conn, name, 0) for name in names]
    end

    # Validate provided mapping keys exist (only when conn was given)
    if conn !== nothing
        for component_id in keys(conn)
            if !(component_id in names)
                throw(ArgumentError("Component ID '$component_id' not found in the cable design."))
            end
        end
    end

    # Warn if all grounded
    !all(iszero, conn_vector) || @warn("At least one component should be assigned to a non-zero phase.")

    # Resolve scalar type and coerce — with identity-preserving pass-through
    T = resolve_T(cable, horz, vert)
    cableT = coerce_to_T(cable, T)
    horzT = (horz isa T) ? horz : coerce_to_T(horz, T)
    vertT = (vert isa T) ? vert : coerce_to_T(vert, T)

    return CablePosition{T}(cableT, horzT, vertT, conn_vector)
end

"""
$(TYPEDEF)

Represents a cable system configuration, defining the physical structure, cables, and their positions.

$(TYPEDFIELDS)
"""
mutable struct LineCableSystem{T<:REALSCALAR}
    "Unique identifier for the system."
    system_id::String
    "Length of the cable system \\[m\\]."
    line_length::T
    "Number of cables in the system."
    num_cables::Int
    "Number of actual phases in the system."
    num_phases::Int
    "Cross-section cable positions."
    cables::Vector{CablePosition{T}}

    @doc """
    $(TYPEDSIGNATURES)

    Constructs a [`LineCableSystem`](@ref) with an initial cable position and system parameters.

    # Arguments

    - `system_id`: Identifier for the cable system.
    - `line_length`: Length of the cable system \\[m\\].
    - `cable`: Initial [`CablePosition`](@ref) object defining a cable position and phase mapping.

    # Returns

    - A [`LineCableSystem`](@ref) object initialized with a single cable position.

    # Examples

    ```julia
    cable_design = CableDesign("example", nominal_data, components_dict)
    cablepos1 = CablePosition(cable_design, 0.0, 0.0, Dict("core" => 1))

    cable_system = $(FUNCTIONNAME)("test_case_1", 1000.0, cablepos1)
    println(cable_system.num_phases)  # Prints number of unique phase assignments
    ```

    # See also

    - [`CablePosition`](@ref)
    - [`CableDesign`](@ref)
    """
    @inline function LineCableSystem{T}(
        system_id::String,
        line_length::T,
        cable::CablePosition{T},
    ) where {T<:REALSCALAR}
        # phase accounting from this single position
        conn = cable.conn
        # count unique non-zero phases
        nph = count(x -> x > 0, unique(conn))
        return new{T}(system_id, line_length, 1, nph, CablePosition{T}[cable])
    end

    @doc """
    $(TYPEDSIGNATURES)

    **Strict numeric kernel**. Builds a typed `LineCableSystem{T}` from a vector of `CablePosition{T}`.
    """
    @inline function LineCableSystem{T}(
        system_id::String,
        line_length::T,
        cables::Vector{CablePosition{T}},
    ) where {T<:REALSCALAR}
        @assert !isempty(cables) "At least one CablePosition must be provided"
        # flatten & count phases
        assigned = unique(vcat((cp.conn for cp in cables)...))
        nph = count(x -> x > 0, assigned)
        return new{T}(system_id, line_length, length(cables), nph, cables)
    end
end

"""
$(TYPEDSIGNATURES)

Weakly-typed constructor. Infers scalar type `T` from `line_length` and the `cable` (or its design), coerces as needed, and calls the strict kernel.
"""
function LineCableSystem(
    system_id::String,
    line_length::Number,
    cable::CablePosition,
)
    T = resolve_T(line_length, cable)
    return LineCableSystem{T}(
        system_id,
        coerce_to_T(line_length, T),
        coerce_to_T(cable, T),
    )
end

"""
$(TYPEDSIGNATURES)

Weakly-typed convenience constructor. Builds a `CablePosition` from a `CableDesign` and coordinates, then constructs the system.
"""
function LineCableSystem(
    system_id::String,
    line_length::Number,
    cable::CableDesign,
    horz::Number,
    vert::Number,
    conn::Union{Dict{String,Int},Nothing}=nothing,
)
    pos = CablePosition(cable, horz, vert, conn)
    return LineCableSystem(system_id, line_length, pos)
end

# Outer (maximum) radius of the last component of a position's design
@inline function _outer_radius(cp::CablePosition)
    comp = cp.design_data.components[end]
    return max(comp.conductor_group.radius_ext, comp.insulator_group.radius_ext)
end

# True if two cable disks overlap (strictly), evaluated in a common scalar type T
@inline function _overlaps(a::CablePosition, b::CablePosition, ::Type{T}) where {T}
    x1 = coerce_to_T(a.horz, T)
    y1 = coerce_to_T(a.vert, T)
    r1 = coerce_to_T(_outer_radius(a), T)
    x2 = coerce_to_T(b.horz, T)
    y2 = coerce_to_T(b.vert, T)
    r2 = coerce_to_T(_outer_radius(b), T)
    d = hypot(x1 - x2, y1 - y2)
    return d < (r1 + r2)   # strict overlap; grazing contact (==) allowed
end

"""
$(TYPEDSIGNATURES)

Adds a new cable position to an existing [`LineCableSystem`](@ref), updating its phase mapping and cable count. If adding the position introduces a different numeric scalar type, the system is **promoted** and the promoted system is returned. Otherwise, mutation happens in place.

# Arguments

- `system`: Instance of [`LineCableSystem`](@ref) to which the cable will be added.
- `cable`: A [`CableDesign`](@ref) object defining the cable structure.
- `horz`: Horizontal coordinate \\[m\\].
- `vert`: Vertical coordinate \\[m\\].
- `conn`: Dictionary mapping component names to phase indices, or `nothing` for automatic assignment.

# Returns

- The modified [`LineCableSystem`](@ref) object with the new cable added.

# Examples

```julia
cable_design = CableDesign("example", nominal_data, components_dict)

# Define coordinates for two cables
xa, ya = 0.0, -1.0
xb, yb = 1.0, -2.0

# Create initial system with one cable
cablepos1 = CablePosition(cable_design, xa, ya, Dict("core" => 1))
cable_system = LineCableSystem("test_case_1", 1000.0, cablepos1)

# Add second cable to system
$(FUNCTIONNAME)(cable_system, cable_design, xb, yb, Dict("core" => 2))

println(cable_system.num_cables)  # Prints: 2
```

# See also

- [`LineCableSystem`](@ref)
- [`CablePosition`](@ref)
- [`CableDesign`](@ref)
"""
function add!(system::LineCableSystem{T}, pos::CablePosition) where {T}
    # Decide the common numeric type first
    Tnew = resolve_T(system, pos)

    # Geometric guard once, in a common type (no mutation, no allocation)
    for cp in system.cables
        if _overlaps(cp, pos, Tnew)
            throw(ArgumentError("Cable position overlaps an existing cable (disks intersect)."))
        end
    end

    if Tnew === T
        posT = coerce_to_T(pos, T)  # identity if already T
        push!(system.cables, posT)
        system.num_cables += 1
        assigned = unique(vcat((cp.conn for cp in system.cables)...))
        system.num_phases = count(x -> x > 0, assigned)
        return system
    else
        @warn """
        Adding a `$Tnew` position to a `LineCableSystem{$T}` returns a **promoted** system.
        Capture the result:  system = add!(system, position)
        """
        sysT = coerce_to_T(system, Tnew)
        posT = coerce_to_T(pos, Tnew)
        push!(sysT.cables, posT)
        sysT.num_cables += 1
        assigned = unique(vcat((cp.conn for cp in sysT.cables)...))
        sysT.num_phases = count(x -> x > 0, assigned)
        return sysT
    end
end

"""
$(TYPEDSIGNATURES)

Convenience `add!` that accepts a cable design and coordinates (and optional mapping).
Builds a [`CablePosition`](@ref) and forwards to `add!(system, pos)`.
"""
function add!(
    system::LineCableSystem{T},
    cable::CableDesign,
    horz::Number,
    vert::Number,
    conn::Union{Dict{String,Int},Nothing}=nothing,
) where {T}
    pos = CablePosition(cable, horz, vert, conn)
    return add!(system, pos)  # may mutate or return a promoted system
end

include("linecablesystem/dataframe.jl")
include("linecablesystem/base.jl")
