"""
$(TYPEDEF)

Hold the numerical input of a computation in the coaxial model. Each conductor of the system
is an equivalent tube, and the dielectric layers of its interval surround it. Each cable has
a center and an outer radius. Its plan comes from [`reductions`](@ref) under the `bundle`
and `kron` options of the formulation. `CoaxialInput(problem, options)` builds the input
once from the blueprints of the designs and from the placements and the connections of the
system. The values have the problem's type `T`, and the conductors follow the terminals of
the system, cable by cable and from the inside out.

The members of a cable that has no conductor around them count as cables.

$(TYPEDFIELDS)
"""
struct CoaxialInput{T <: Real, L <: AbstractVector{<:Material{T}}}
    "Cable of each conductor."
    cable::Vector{Int}
    "Conductor index of each conductor in the circuit, `k ≥ 1` or zero."
    connection::Vector{Int}
    "Label of each conductor, `cable:<i>:<terminal>`."
    label::Vector{String}
    "Inner radius of the tube of each conductor \\[m\\]."
    r_in::Vector{T}
    "Outer radius of the tube of each conductor \\[m\\]."
    r_ex::Vector{T}
    "Material of the tube of each conductor."
    material::Vector{Material{T, :conductor}}
    """
    For each conductor, the outer radius of its dielectric interval, or the conductor's own
    outer radius for an empty interval \\[m\\].
    """
    r_ins_ex::Vector{T}
    """
    Relative permeability of the dielectric interval of each conductor, corrected for the
    helical path of its current \\[dimensionless\\].
    """
    mu_r_ins::Vector{T}
    "Range of the dielectric layers in the interval of each conductor."
    dielectric_ranges::Vector{UnitRange{Int}}
    "Inner radius of each dielectric layer \\[m\\]."
    dielectric_r_in::Vector{T}
    "Outer radius of each dielectric layer \\[m\\]."
    dielectric_r_ex::Vector{T}
    "Material of each dielectric layer."
    dielectric_material::L
    "Range of the conductors of each concentric assembly."
    assemblies::Vector{UnitRange{Int}}
    "Assemblies in the interior of each assembly, its members in the blueprints."
    members::Vector{Vector{Int}}
    "Axis of each assembly in the frame of its cable, in \\[m\\]."
    positions::Vector{Tuple{T, T}}
    "Range of the interior layers of each assembly, from its filler around its members."
    interior_ranges::Vector{UnitRange{Int}}
    "Horizontal coordinate of the center of each cable in the frame of the system \\[m\\]."
    x::Vector{T}
    "Vertical coordinate of the center of each cable in the frame of the system \\[m\\]."
    y::Vector{T}
    "Outer radius \\[m\\] of each cable."
    radius::Vector{T}
    "Layer of the earth model at the center of each cable: 1 for the air, `k` for the `(k - 1)`-th earth layer."
    earth_layer::Vector{Int}
    "Conductor pairs `(first, other)` to merge with [`bundle!`](@ref), or an empty list."
    bundles::Vector{Tuple{Int, Int}}
    """
    Indices to keep in result order, with the first conductor of each positive index first.
    """
    kept::Vector{Int}
    "Indices to remove: other conductors with `bundle`, index zero with `kron`."
    eliminated::Vector{Int}
    """
    Label of each kept row: the label of its conductor, or the labels of a bundle joined by
    `+`.
    """
    terminals::Vector{String}
    function CoaxialInput(
            problem::LineParametersProblem{T}, options::FormulationOptions) where {T <:
                                                                                   Real}
        system = problem.system
        blueprints = [flatten(design) for design in system.designs]
        tubes = [tube for blueprint in blueprints for tube in blueprint.conductors]
        dielectrics = [dielectric for blueprint in blueprints
                       for dielectric in blueprint.dielectrics]
        tube_starts = cumsum([0;
                              [length(blueprint.conductors) for blueprint in blueprints]])
        dielectric_starts = cumsum([0;
                                    [length(blueprint.dielectrics)
                                     for blueprint in blueprints]])
        connection = [row.connection for row in system.terminals]
        (; label, bundles, kept, eliminated, terminals) = reductions(system, options)
        r_in = [convert(T, tube.r_in) for tube in tubes]
        r_ex = [convert(T, tube.r_ex) for tube in tubes]
        material = [convert(Material{T}, tube.material) for tube in tubes]
        mu_r_ins = [convert(T, mu) for blueprint in blueprints for mu in blueprint.mu_r_ins]
        dielectric_ranges = [span .+ dielectric_starts[k]
                             for (k, blueprint) in enumerate(blueprints)
                             for span in blueprint.dielectric_ranges]
        dielectric_r_in = [convert(T, dielectric.r_in) for dielectric in dielectrics]
        dielectric_r_ex = [convert(T, dielectric.r_ex) for dielectric in dielectrics]
        r_ins_ex = [isempty(span) ? outer : dielectric_r_ex[last(span)]
                    for (span, outer) in zip(dielectric_ranges, r_ex)]
        materials = [convert(Material{T}, dielectric.material)
                     for dielectric in dielectrics]
        dielectric_material = collect(Union{map(typeof, materials)...}, materials)
        assemblies = [span .+ tube_starts[k]
                      for (k, blueprint) in enumerate(blueprints)
                      for span in blueprint.assembly_ranges]
        assembly_starts = cumsum([0;
                                  [length(blueprint.assembly_ranges)
                                   for blueprint in blueprints]])
        members = [indices .+ assembly_starts[k]
                   for (k, blueprint) in enumerate(blueprints)
                   for indices in blueprint.members]
        positions = [convert.(T, position)
                     for blueprint in blueprints for position in blueprint.positions]
        interior_ranges = [span .+ dielectric_starts[k]
                           for (k, blueprint) in enumerate(blueprints)
                           for span in blueprint.interior_ranges]
        # The cables, each the row of its blueprint that the earth formulas see: the root, or
        # each member of a root without a conductor.
        rows = [(k, row) for (k, blueprint) in enumerate(blueprints)
                for row in (isempty(first(blueprint.assembly_ranges)) &&
                     !isempty(first(blueprint.members)) ? first(blueprint.members) : [1])]
        centers = [system.placements[k] * Placement(blueprints[k].positions[row]...)
                   for (k, row) in rows]
        x = [convert(T, center.x) for center in centers]
        y = [convert(T, center.y) for center in centers]
        radius = [isone(row) ? convert(T, support(system.designs[k].outer)) :
                  r_ins_ex[last(assemblies[assembly_starts[k] + row])]
                  for (k, row) in rows]
        # Each cable's conductors end at its last one, and the cables follow the conductors.
        ends = [isone(row) ? tube_starts[k + 1] :
                last(assemblies[assembly_starts[k] + row])
                for (k, row) in rows]
        cable = [searchsortedfirst(ends, i) for i in eachindex(tubes)]
        depths = cumsum([stratum.thickness for stratum in Base.tail(problem.earth.layers)])
        earth_layer = [center > 0 ? 1 : 1 + searchsortedfirst(depths, -center)
                       for center in y]
        return new{T, typeof(dielectric_material)}(
            cable, connection, label, r_in, r_ex, material,
            r_ins_ex, mu_r_ins, dielectric_ranges, dielectric_r_in,
            dielectric_r_ex, dielectric_material, assemblies,
            members, positions, interior_ranges, x, y, radius, earth_layer, bundles, kept,
            eliminated, terminals)
    end
end

# The tubes with a hole that the coaxial input reduces its conductors to.
function domain(input::CoaxialInput{T}) where {T <: Real}
    N = typeof(nominal(zero(T)))
    none = (typemax(N), typemin(N))
    wall = extrema(
        nominal((input.r_ex[i] - input.r_in[i]) / (input.r_ex[i] + input.r_in[i]))
        for i in eachindex(input.r_in) if !iszero(input.r_in[i]); init = none)
    return (wall = ClosedInterval(wall...),)
end
