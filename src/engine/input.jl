"""
$(TYPEDEF)

Hold the numerical input of a computation in the coaxial model. Each conductor of the system
is an equivalent tube, and the dielectric layers of its interval surround it. Each cable has
a center and an outer radius. Its plan comes from [`reductions`](@ref) under the `bundle`
and `kron` options of the formulation. `CoaxialInput(problem, options)` builds the input
once from the flattened designs and from the placements and the connections of the system.
The values have the problem's type `T`, and the conductors follow the terminals of the
system, cable by cable and from the inside out.

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
    "Assemblies in the interior of each assembly with a conductor, empty for an assembly without one."
    members::Vector{Vector{Int}}
    "Range of the conductors that the outermost conductor of each assembly encloses, itself included."
    enclosed::Vector{UnitRange{Int}}
    "Axis of each assembly in the frame of its cable \\[m\\]."
    positions::Vector{Tuple{T, T}}
    "Range of the interior layers of each assembly, from its filler around its members."
    interior_ranges::Vector{UnitRange{Int}}
    "Range of the layers around each assembly: the interior layers of the assembly that holds it, empty for an assembly in no interior."
    exterior_ranges::Vector{UnitRange{Int}}
    "Horizontal coordinate of the center of each cable in the frame of the system \\[m\\]."
    x::Vector{T}
    "Vertical coordinate of the center of each cable in the frame of the system \\[m\\]."
    y::Vector{T}
    """
    The outer radius of the tube of each cable: the outer radius of the dielectric interval
    of its last conductor \\[m\\].
    """
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
    @doc """
    $(TYPEDSIGNATURES)

    Read the input from the flattened stack of each cable design, as [`flatten`](@ref) gives
    it. A stack is one assembly of conductors about one axis. The cable lies about its
    origin. A member of an assembly or of a repeat lies about the centroid of its first
    conductor. A member without a conductor lies at its placement, and it gives its members
    and its regions to the assembly around it.

    A terminal over a circle or an annulus is its tube. Any other terminal takes Ametani's
    cylinder of the area ``S`` and the perimeter ``l`` of its region:

    ```math
    r_e = \\frac{l}{2\\pi}, \\qquad r_i = \\sqrt{r_e^2 - \\frac{S}{\\pi}}.
    ```

    An inner radius within `tolerance(r_e)` of zero, decided on `nominal`, is zero. The
    cylinder keeps the resistance of the region. Its relative permeability gives it the
    geometric mean radius ``g`` of the region, with ``S_0`` its own section:

    ```math
    \\mu = \\frac{\\log r_e - \\log g}{\\log r_e - \\log \\operatorname{gmd}(S_0, S_0)}.
    ```

    The dielectric regions after a terminal in its stack are the layers of its interval. A
    region that an offset declares around another outline than a circle is the annulus of
    its thickness outside the layer before it. A tube that such an offset declares around
    a filler other than a circle is the first conductor of its stack. It takes Ametani's
    cylinder of its own cross-section. The interval of a tube reaches from ``r_e``
    to ``r_d``, and its layer ``k`` of relative permeability ``\\mu_k`` reaches from ``a_k``
    to ``b_k``. The current of a laid conductor makes ``N`` turns per unit length, which the
    helix of its lay gives at the outer radius of its tube. The interval has

    ```math
    \\mu_{ins} = \\frac{\\sum_k \\mu_k \\log(b_k / a_k)}{\\log(r_d / r_e)}
        \\left[1 + \\frac{2 \\pi^2 N^2 (r_d^2 - r_e^2)}{\\log(r_d / r_e)}\\right],
    ```

    with the last factor evaluated as ``1 + 4 \\pi^2 N^2 r_d r_e \\sinh(u) / u``, where
    ``u = \\log(r_d / r_e)``, and ``\\sinh(u) / u`` at its limit 1 for ``r_d = r_e``. An
    interval without layers has ``\\mu_{ins} = 1``.

    The interior of an assembly with members begins with its filler. The filler is one
    layer of the fill's material, from the axis to the container of the enclosure. A
    container that is not round ends at the outer radius of Ametani's cylinder,
    ``l / 2\\pi`` for its perimeter ``l``.

    The regions of the stack inside the enclosure join the filler, and so do the regions of
    a member without a conductor. The filler then takes their composite [`Material`](@ref),
    each with its fraction of the area. The annuli between the enclosure and the first
    conductor of the stack are the other layers of the interior. A region inside the first
    conductor of a stack without members lies in no interval.

    # Errors

    - Throws the errors of [`flatten`](@ref).
    """
    function CoaxialInput(
            problem::LineParametersProblem{T}, options::FormulationOptions) where {T <:
                                                                                   Real}
        system = problem.system
        rows = (r_in = T[], r_ex = T[], material = Material{T, :conductor}[], turns = T[],
            dielectric_ranges = UnitRange{Int}[], dielectric_r_in = T[],
            dielectric_r_ex = T[], dielectric_material = Material{T}[],
            assemblies = UnitRange{Int}[], members = Vector{Int}[],
            positions = Tuple{T, T}[], interior_ranges = UnitRange{Int}[],
            fillers = Vector{Pair{Material{T}, T}}[])
        roots = Int[]
        tube_ends = Int[]
        for design in system.designs
            push!(roots, length(rows.assemblies) + 1)
            resolve(rows, CoaxialInput, design)
            push!(tube_ends, length(rows.r_in))
        end
        (; r_in, r_ex, material, dielectric_ranges, dielectric_r_in, dielectric_r_ex,
            assemblies, positions, interior_ranges) = rows
        connection = [row.connection for row in system.terminals]
        (; label, bundles, kept, eliminated, terminals) = reductions(system, options)
        r_ins_ex = [isempty(span) ? outer : dielectric_r_ex[last(span)]
                    for (span, outer) in zip(dielectric_ranges, r_ex)]
        dielectric_material = collect(
            Union{map(typeof, rows.dielectric_material)...}, rows.dielectric_material)
        mu_r_ins = map(eachindex(r_ex)) do i
            span = dielectric_ranges[i]
            isempty(span) && return one(T)
            u = log(r_ins_ex[i] / r_ex[i])
            μ = dielectric_material[first(span)].mu_r
            return (μ +
                    sum(
                l -> (dielectric_material[l].mu_r - μ) *
                     log(dielectric_r_ex[l] / dielectric_r_in[l]) / u,
                span[2:end]; init = zero(T))) *
                   (1 +
                    4π^2 * rows.turns[i]^2 * r_ins_ex[i] * r_ex[i] *
                    (iszero(nominal(u)) ? one(u) : sinh(u) / u))
        end
        # The cables that the earth formulas see: the root of a design, or each member of a
        # root without a conductor.
        cables = [(k, row)
                  for (k, root) in enumerate(roots)
                  for row in (isempty(assemblies[root]) && !isempty(rows.members[root]) ?
                       rows.members[root] : [root])]
        members = [isempty(assemblies[a]) ? Int[] : rows.members[a]
                   for a in eachindex(assemblies)]
        enclosed = map(eachindex(assemblies)) do a
            head = a
            while !isempty(members[head])
                head = first(members[head])
            end
            return first(assemblies[head]):last(assemblies[a])
        end
        exterior_ranges = fill(1:0, length(assemblies))
        for (c, indices) in enumerate(members), m in indices

            exterior_ranges[m] = interior_ranges[c]
        end
        centers = [system.placements[k] * Placement(positions[row]...)
                   for (k, row) in cables]
        x = [convert(T, center.x) for center in centers]
        y = [convert(T, center.y) for center in centers]
        radius = [r_ins_ex[last(assemblies[row])] for (k, row) in cables]
        # Each cable's conductors end at its last one, and the cables follow the conductors.
        ends = [row == roots[k] ? tube_ends[k] : last(assemblies[row])
                for (k, row) in cables]
        cable = [searchsortedfirst(ends, i) for i in eachindex(r_in)]
        depths = cumsum([stratum.thickness for stratum in Base.tail(problem.earth.layers)])
        earth_layer = [center > 0 ? 1 : 1 + searchsortedfirst(depths, -center)
                       for center in y]
        return new{T, typeof(dielectric_material)}(
            cable, connection, label, r_in, r_ex, material,
            r_ins_ex, mu_r_ins, dielectric_ranges, dielectric_r_in,
            dielectric_r_ex, dielectric_material, assemblies,
            members, enclosed, positions, interior_ranges, exterior_ranges,
            x, y, radius, earth_layer, bundles, kept,
            eliminated, terminals)
    end
end

# The cable is the first assembly of its design, about its origin. `rows` holds the rows of
# the input of a system, and each part of a flattened stack adds its own.
function resolve(rows::NamedTuple, ::Type{CoaxialInput}, design::CableDesign)
    T = eltype(rows.r_in)
    l = length(rows.dielectric_r_in)
    push!(rows.assemblies, 1:0)
    push!(rows.members, Int[])
    push!(rows.fillers, Pair{Material{T}, T}[])
    push!(rows.interior_ranges, (l + 1):l)
    push!(rows.positions, (zero(T), zero(T)))
    a = length(rows.assemblies)
    foreach(item -> resolve(rows, CoaxialInput, Placement(0, 0), a, a, item),
        flatten(design).items)
    if isempty(rows.assemblies[a])
        n = length(rows.r_in)
        rows.assemblies[a] = (n + 1):n
    end
    return rows
end

# A member in the frame `frame` of its cable is an assembly about the centroid of its first
# conductor, and a member of the assembly `host`. A member without a conductor lies at its
# placement. Its members and its regions are those of `host`.
function resolve(rows::NamedTuple, ::Type{CoaxialInput}, frame::Placement, host::Int,
        stack::Stack)
    T = eltype(rows.r_in)
    l = length(rows.dielectric_r_in)
    push!(rows.assemblies, 1:0)
    push!(rows.members, Int[])
    push!(rows.fillers, Pair{Material{T}, T}[])
    push!(rows.interior_ranges, (l + 1):l)
    a = length(rows.assemblies)
    k = findfirst(item -> item isa Terminal, stack.items)
    if k === nothing
        push!(rows.positions, (frame.x, frame.y))
        foreach(item -> resolve(rows, CoaxialInput, frame, a, host, item), stack.items)
        n, l = length(rows.r_in), length(rows.dielectric_r_in)
        rows.assemblies[a] = (n + 1):n
        rows.interior_ranges[a] = (l + 1):l
    else
        push!(rows.positions, frame(centroid(boundary(stack.items[k]))))
        foreach(item -> resolve(rows, CoaxialInput, frame, a, a, item), stack.items)
        push!(rows.members[host], a)
    end
    return rows
end

# The tube of a round conductor is its region. It is the next conductor of the assembly
# `a`, whose interior ends where its first conductor begins.
function resolve(rows::NamedTuple, ::Type{CoaxialInput}, ::Placement, a::Int, ::Int,
        terminal::Terminal{<:Region{<:Union{Circle, Annulus}}})
    region = terminal.part
    push!(rows.r_in, r_in(region.primitive))
    push!(rows.r_ex, r_ex(region.primitive))
    push!(rows.material, convert(Material{eltype(rows.r_in)}, region.material))
    push!(rows.turns, 0)
    n, l = length(rows.r_in), length(rows.dielectric_r_in)
    push!(rows.dielectric_ranges, (l + 1):l)
    if isempty(rows.assemblies[a])
        isempty(rows.interior_ranges[a]) && (rows.interior_ranges[a] = (l + 1):l)
        rows.assemblies[a] = n:n
    else
        rows.assemblies[a] = first(rows.assemblies[a]):n
    end
    return rows
end

# Any other conductor takes Ametani's cylinder of the area and the perimeter of its
# cross-section `S`.
function resolve(rows::NamedTuple, ::Type{CoaxialInput}, frame::Placement, a::Int,
        owner::Int, terminal::Terminal{<:Region})
    return resolve(rows, CoaxialInput, frame, a, owner, terminal, terminal.part.primitive)
end
function resolve(rows::NamedTuple, ::Type{CoaxialInput}, frame::Placement, a::Int,
        owner::Int, terminal::Terminal{<:Region}, S::AbstractShape)
    r_e = perimeter(boundary(S)) / 2π
    r_i = nominal(r_e^2 - area(S) / π) <= tolerance(r_e)^2 ? zero(r_e) :
          sqrt(r_e^2 - area(S) / π)
    return resolve(rows, CoaxialInput, frame, a, owner, terminal, S,
        iszero(nominal(r_i)) ? Circle(r_e) : Annulus(r_i, r_e))
end

# A part that an offset declares around another outline than a circle begins where the
# layer before it ends. `before` holds the outer radii of what lies before it in its
# assembly. In an assembly with a conductor, they are the layers of its last conductor, or
# its conductors. Before the first conductor, they are the layers of the interior. A part
# around a filler other than a circle has none.
function resolve(rows::NamedTuple, ::Type{CoaxialInput}, frame::Placement, a::Int,
        owner::Int,
        part::Union{Region{<:Offset},
            Terminal{<:CablePart{<:Region{<:Offset}, <:PerforatedShape}}})
    tubes = rows.assemblies[a]
    span = isempty(tubes) ? rows.interior_ranges[a] : rows.dielectric_ranges[last(tubes)]
    return resolve(rows, CoaxialInput, frame, a, owner, part,
        isempty(span) ? view(rows.r_ex, tubes) : view(rows.dielectric_r_ex, span))
end

# Such a conductor is the tube of its thickness outside the last of the radii `before`.
# With none before it, it takes Ametani's cylinder of its own cross-section.
function resolve(rows::NamedTuple, ::Type{CoaxialInput}, frame::Placement, a::Int,
        owner::Int, terminal::Terminal{<:CablePart{<:Region{<:Offset}, <:PerforatedShape}},
        before::AbstractVector{<:Real})
    placed = terminal.part
    conductor = Terminal(terminal.name, placed.source)
    isempty(before) &&
        return resolve(rows, CoaxialInput, frame, a, owner, conductor, placed.shape)
    r = last(before)
    return resolve(rows, CoaxialInput, frame, a, owner, conductor, placed.shape,
        Annulus(r, r + placed.source.primitive.t))
end

# Such a dielectric layer is the annulus of its thickness outside the last of the radii
# `before`. With none before it, it has no radius to begin at.
function resolve(rows::NamedTuple, ::Type{CoaxialInput}, frame::Placement, a::Int,
        owner::Int, region::Region{<:Offset}, before::AbstractVector{<:Real})
    isempty(before) && throw(ArgumentError(
        "the region :$(region.tag) declares an Offset from the outline of a filler that is not round; " *
        "the coaxial input begins such a layer at the outer radius of a conductor or of a layer around one"))
    r = last(before)
    return resolve(rows, CoaxialInput, frame, a, owner,
        Region(region.tag, Annulus(r, r + region.primitive.t), region.material))
end

# The tube `tube` stands for a conductor of the cross-section `S`. Its material gives it the
# resistance and the geometric mean radius of `S`.
function resolve(rows::NamedTuple, ::Type{CoaxialInput}, frame::Placement, a::Int,
        owner::Int, terminal::Terminal{<:Region}, S::AbstractShape,
        tube::Union{Circle, Annulus})
    medium = terminal.part.material
    return resolve(rows, CoaxialInput, frame, a, owner,
        Terminal(terminal.name,
            Region(terminal.part.tag, tube, medium.rho / area(S), gmr(S, medium.mu_r),
                medium.T0, medium.alpha)))
end

# The helix of a laid conductor gives its turns at the outer radius of its tube.
function resolve(rows::NamedTuple, ::Type{CoaxialInput}, frame::Placement, a::Int,
        owner::Int, terminal::Terminal{<:Lay})
    (; helix, part) = terminal.part
    resolve(rows, CoaxialInput, frame, a, owner, Terminal(terminal.name, part))
    rows.turns[end] = helix.dir / pitch(helix, rows.r_ex[end])
    return rows
end

# A dielectric annulus after a conductor of its stack is a layer of the interval of that
# conductor. Before the first conductor and outside the members of the stack, it is a layer
# of the interior. Before the members, or in a stack without a conductor, it joins the
# filler of the assembly `owner`.
function resolve(rows::NamedTuple, ::Type{CoaxialInput}, ::Placement, a::Int, owner::Int,
        region::Region{<:Annulus})
    T = eltype(rows.r_in)
    if owner != a || isempty(rows.assemblies[a]) && isempty(rows.members[a])
        push!(rows.fillers[owner],
            convert(Material{T}, region.material) => area(region.primitive))
        return rows
    end
    push!(rows.dielectric_r_in, r_in(region.primitive))
    push!(rows.dielectric_r_ex, r_ex(region.primitive))
    push!(rows.dielectric_material, convert(Material{T}, region.material))
    l = length(rows.dielectric_r_in)
    if isempty(rows.assemblies[a])
        rows.interior_ranges[a] = (isempty(rows.interior_ranges[a]) ? l :
                                   first(rows.interior_ranges[a])):l
    else
        i = last(rows.assemblies[a])
        rows.dielectric_ranges[i] = first(rows.dielectric_ranges[i]):l
    end
    return rows
end

# Any other dielectric region declared by a primitive lies inside the first conductor or the
# members of its stack, and joins the filler.
function resolve(rows::NamedTuple, ::Type{CoaxialInput}, ::Placement, ::Int, owner::Int,
        region::Region{<:AbstractPrimitive})
    push!(rows.fillers[owner],
        convert(Material{eltype(rows.r_in)}, region.material) => area(region.primitive))
    return rows
end

# The filler of an enclosure ends at the outer radius of the circle that stands for its
# container. A circle or an annulus is its own circle. Any other container takes the outer
# radius of Ametani's cylinder, from its perimeter.
function resolve(rows::NamedTuple, ::Type{CoaxialInput}, frame::Placement, a::Int,
        owner::Int, enclosure::Enclosure{<:CablePart})
    return resolve(rows, CoaxialInput, frame, a, owner, enclosure,
        enclosure.fill.shape.outer)
end
function resolve(rows::NamedTuple, ::Type{CoaxialInput}, frame::Placement, a::Int,
        owner::Int, enclosure::Enclosure{<:CablePart}, container::AbstractShape)
    return resolve(rows, CoaxialInput, frame, a, owner, enclosure,
        Circle(perimeter(boundary(container)) / 2π))
end

# An enclosure gives its members, then its fill. In the stack of an assembly the fill is the
# filler of the interior, from the axis to the circle `container`. The regions that joined
# the filler share its material by their areas. The cross-section of the resolved fill is its
# container without the outlines of its members.
function resolve(rows::NamedTuple, ::Type{CoaxialInput}, frame::Placement, a::Int,
        owner::Int, enclosure::Enclosure{<:CablePart}, container::Union{Circle, Annulus})
    T = eltype(rows.r_in)
    resolve(rows, CoaxialInput, frame, a, owner, enclosure.contents)
    fill, media = enclosure.fill, rows.fillers[owner]
    (owner != a || !isempty(media)) &&
        push!(media, convert(Material{T}, fill.source.material) => area(fill.shape))
    owner == a || return rows
    push!(rows.dielectric_r_in, zero(T))
    push!(rows.dielectric_r_ex, r_ex(container))
    push!(rows.dielectric_material,
        isempty(media) ? convert(Material{T}, fill.source.material) :
        Material(first.(media), last.(media) ./ sum(last, media)))
    l = length(rows.dielectric_r_in)
    rows.interior_ranges[a] = l:l
    return rows
end

# Each member of an assembly and each copy of a repeat is a member at its placement.
function resolve(rows::NamedTuple, ::Type{CoaxialInput}, frame::Placement, ::Int,
        owner::Int, assembly::Assembly)
    foreach(
        member -> resolve(rows, CoaxialInput, frame * member.first, owner, member.second),
        assembly.members)
    return rows
end
function resolve(rows::NamedTuple, ::Type{CoaxialInput}, frame::Placement, ::Int,
        owner::Int, repeat::Repeat)
    foreach(at -> resolve(rows, CoaxialInput, frame * at, owner, repeat.part),
        repeat.pattern)
    return rows
end

# The walls of tubes with a hole, and m·D from the centers and tube radii of buried cables.
function domain(input::CoaxialInput{T}, problem::LineParametersProblem) where {T <: Real}
    N = typeof(nominal(zero(T)))
    none = (typemax(N), typemin(N))
    wall = extrema(
        nominal((input.r_ex[i] - input.r_in[i]) / (input.r_ex[i] + input.r_in[i]))
        for i in eachindex(input.r_in) if !iszero(input.r_in[i]); init = none)
    (; earth, frequencies) = problem
    deepest = last(earth.layers)
    buried = findall(i -> nominal(input.y[i]) < 0, eachindex(input.y))
    reach = extrema(
        (i == j ? nominal(input.radius[i]) :
         nominal(hypot(input.x[i] - input.x[j], input.y[i] - input.y[j]))
        for i in buried, j in buried); init = none)
    m = map(f -> nominal(sqrt(2π * f * μ₀ * deepest.mu_r / deepest.rho)),
        (first(frequencies), last(frequencies)))
    return (wall = ClosedInterval(wall...),
        mD = isempty(buried) ? ClosedInterval(none...) :
             ClosedInterval(first(m) * first(reach), last(m) * last(reach)))
end
