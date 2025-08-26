import ..LineCableModels: add!

"""
$(TYPEDEF)

Represents a physically defined cable with position and phase mapping within a system.

$(TYPEDFIELDS)
"""
struct CablePosition
    "The [`CableDesign`](@ref) object assigned to this cable position."
    design_data::CableDesign
    "Horizontal coordinate \\[m\\]."
    horz::Number
    "Vertical coordinate \\[m\\]."
    vert::Number
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
    function CablePosition(
        cable::Union{CableDesign,Nothing},
        horz::Number,
        vert::Number,
        conn::Union{Dict{String,Int},Nothing}=nothing,
    )
        # Validate cable design is not empty
        @assert !isnothing(cable) "A valid CableDesign must be provided"
        @assert !isempty(cable.components) "CableDesign must contain at least one component"

        # Find outermost radius by checking both conductor and insulator groups of last component
        last_comp = cable.components[end]
        r_cond = last_comp.conductor_group.radius_ext
        r_ins = last_comp.insulator_group.radius_ext
        r_max = max(r_cond, r_ins)

        # Validate vertical position
        @assert vert != 0 "Vertical position cannot be exactly at the air/earth interface (z=0)"
        @assert abs(vert) >= r_max """
        Vertical position |$(vert)| must be greater than or equal to cable's outer radius $(r_max)
        to prevent crossing the air/earth interface at z=0
        """


        # Create phase mapping vector
        components = [comp.id for comp in cable.components]
        if isnothing(conn)
            conn_vector = [i == 1 ? 1 : 0 for i in 1:length(components)]  # Default: First component gets phase 1
        else
            conn_vector = [get(conn, name, 0) for name in components]  # Ensure correct mapping order
        end
        # Validate there is at least one ungrounded conductor
        !all(iszero, conn_vector) || @warn ("At least one component must be assigned to a non-zero phase.")

        for component_id in keys(conn)
            if !(component_id in [c.id for c in cable.components])
                throw(ArgumentError("Component ID '$component_id' not found in the cable design."))
            end
        end


        return new(cable, horz, vert, conn_vector)
    end
end

"""
$(TYPEDEF)

Represents a cable system configuration, defining the physical structure, cables, and their positions.

$(TYPEDFIELDS)
"""
mutable struct LineCableSystem
    "Unique identifier for the system."
    system_id::String
    "Length of the cable system \\[m\\]."
    line_length::Number
    "Number of cables in the system."
    num_cables::Int
    "Number of actual phases in the system."
    num_phases::Int
    "Cross-section cable positions."
    cables::Vector{CablePosition}

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
    function LineCableSystem(
        system_id::String,
        line_length::Number,
        cable::CablePosition,
    )
        # Initialize with the first cable definition
        num_cables = 1

        # Count unique nonzero phases from the first cable
        assigned_phases = unique(cable.conn)
        num_phases = count(x -> x > 0, assigned_phases)

        return new(system_id, line_length, num_cables, num_phases, [cable])
    end
end

"""
$(TYPEDSIGNATURES)

Adds a new cable position to an existing [`LineCableSystem`](@ref), updating its phase mapping and cable count.

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
function add!(
    system::LineCableSystem,
    cable::CableDesign,
    horz::Number,
    vert::Number,
    conn::Union{Dict{String,Int},Nothing}=nothing,
)
    max_phase =
        isempty(system.cables) ? 0 : maximum(maximum.(getfield.(system.cables, :conn)))

    component_names = [comp.id for comp in cable.components]  # Get component IDs from vector

    new_conn = if isnothing(conn)
        Dict(name => (i == 1 ? max_phase + 1 : 0) for (i, name) in enumerate(component_names))
    else
        Dict(name => get(conn, name, 0) for name in component_names)  # Ensures correct mapping order
    end

    # Validate that the coordinates do not overlap with existing cables
    for cable_pos in system.cables
        if cable_pos.horz == horz && cable_pos.vert == vert
            throw(ArgumentError("Cable position overlaps with existing cable"))
        end
    end

    push!(system.cables, CablePosition(cable, horz, vert, new_conn))

    # Update num_cables
    system.num_cables += 1

    # Update num_phases by counting unique nonzero phases
    assigned_phases = unique(vcat([cable_pos.conn for cable_pos in system.cables]...))
    system.num_phases = count(x -> x > 0, assigned_phases)
    system
end