"""
$(TYPEDEF)

Hold the numerical input of a computation in the coaxial model. Each conductor of the system
is an equivalent tube, and the dielectric layers of its interval surround it. Each cable has
a center and an outer radius. Its plan comes from [`reductions`](@ref) under the `bundle`
and `kron` options of the formulation. `CoaxialInput(problem, options)` builds the input
once from the blueprints of the designs and from the placements and the connections of the
system. The values have the problem's type `T`, and the conductors follow the terminals of
the system, cable by cable and from the inside out.

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
    For each conductor, the relative permeability of its dielectric interval, with the
    correction for the helical current path \\[dimensionless\\].
    """
    mu_r_ins::Vector{T}
    "Range of the dielectric layers in the interval of each conductor."
    layers::Vector{UnitRange{Int}}
    "Inner radius of each dielectric layer \\[m\\]."
    layer_r_in::Vector{T}
    "Outer radius of each dielectric layer \\[m\\]."
    layer_r_ex::Vector{T}
    "Material of each dielectric layer."
    layer_material::L
    "Range of the conductors of each concentric assembly."
    assemblies::Vector{UnitRange{Int}}
    "Horizontal coordinate of the center of each cable in the frame of the system \\[m\\]."
    x::Vector{T}
    "Vertical coordinate of the center of each cable in the frame of the system \\[m\\]."
    y::Vector{T}
    "Outer radius of each cable \\[m\\]."
    radius::Vector{T}
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
        dielectrics = [layer for blueprint in blueprints for layer in blueprint.dielectrics]
        tube_starts = cumsum([0;
                              [length(blueprint.conductors) for blueprint in blueprints]])
        layer_starts = cumsum([0;
                               [length(blueprint.dielectrics) for blueprint in blueprints]])
        cable = [row.cable for row in system.terminals]
        connection = [row.connection for row in system.terminals]
        (; label, bundles, kept, eliminated, terminals) = reductions(system, options)
        r_in = [convert(T, tube.r_in) for tube in tubes]
        r_ex = [convert(T, tube.r_ex) for tube in tubes]
        material = [convert(Material{T}, tube.material) for tube in tubes]
        mu_r_ins = [convert(T, mu) for blueprint in blueprints for mu in blueprint.mu_r_ins]
        layers = [span .+ layer_starts[k]
                  for (k, blueprint) in enumerate(blueprints)
                  for span in blueprint.dielectric_ranges]
        layer_r_in = [convert(T, layer.r_in) for layer in dielectrics]
        layer_r_ex = [convert(T, layer.r_ex) for layer in dielectrics]
        r_ins_ex = [isempty(span) ? outer : layer_r_ex[last(span)]
                    for (span, outer) in zip(layers, r_ex)]
        materials = [convert(Material{T}, layer.material) for layer in dielectrics]
        layer_material = collect(Union{map(typeof, materials)...}, materials)
        assemblies = [span .+ tube_starts[k]
                      for (k, blueprint) in enumerate(blueprints)
                      for span in blueprint.assembly_ranges]
        centers = [system.placements[k] * Placement(only(blueprint.positions)...)
                   for (k, blueprint) in enumerate(blueprints)]
        x = [convert(T, center.x) for center in centers]
        y = [convert(T, center.y) for center in centers]
        radius = [convert(T, support(design.outer)) for design in system.designs]
        return new{T, typeof(layer_material)}(
            cable, connection, label, r_in, r_ex, material,
            r_ins_ex, mu_r_ins, layers, layer_r_in, layer_r_ex, layer_material, assemblies,
            x, y, radius, bundles, kept,
            eliminated, terminals)
    end
end
