function _identity_pose(pose::Pose2)
    return iszero(pose.x) && iszero(pose.y) && iszero(pose.φ)
end

function _display_pose(pose::Pose2)
    _identity_pose(pose) && return nothing
    return sprint(show, pose; context = :compact => true)
end

function _bounded_collection(value)
    value isa AbstractArray && return "$(length(value)) values"
    value isa Tuple && return "$(length(value)) values"
    return sprint(show, value; context = :compact => true)
end

TextDisplay.@showfields Pose2 "Pose" pose -> (
    x = TextDisplay.engineering(pose.x, :meter),
    y = TextDisplay.engineering(pose.y, :meter),
    φ = TextDisplay.angle(pose.φ)
)

TextDisplay.@showfields Disk "Disk" primitive -> (
    r = TextDisplay.engineering(primitive.r, :meter),
    at = _display_pose(primitive.at)
)

TextDisplay.@showfields Rectangle "Rectangle" primitive -> (
    w = TextDisplay.engineering(primitive.w, :meter),
    h = TextDisplay.engineering(primitive.h, :meter),
    at = _display_pose(primitive.at)
)

TextDisplay.@showfields Ellipse "Ellipse" primitive -> (
    a = TextDisplay.engineering(primitive.a, :meter),
    b = TextDisplay.engineering(primitive.b, :meter),
    at = _display_pose(primitive.at)
)

TextDisplay.@showfields Sector "Sector" primitive -> (
    rᵢ = TextDisplay.engineering(primitive.ri, :meter),
    rₒ = TextDisplay.engineering(primitive.ro, :meter),
    φ₀ = iszero(primitive.φ0) ? nothing : TextDisplay.angle(primitive.φ0),
    Δφ = TextDisplay.angle(primitive.span),
    at = _display_pose(primitive.at)
)

TextDisplay.@showfields Annulus "Annulus" primitive -> (
    rᵢ = TextDisplay.engineering(primitive.ri, :meter),
    rₒ = TextDisplay.engineering(primitive.ro, :meter),
    at = _display_pose(primitive.at)
)

TextDisplay.@showfields Polygon "Polygon" primitive -> (
    vertices = length(primitive.points),
    at = _display_pose(primitive.at)
)

TextDisplay.@showfields RoundedSector "RoundedSector" primitive -> (
    Δφ = TextDisplay.angle(primitive.span),
    rᵢ = TextDisplay.engineering(primitive.r_base, :meter),
    rₒ = TextDisplay.engineering(primitive.r_back, :meter),
    fillet = iszero(primitive.fillet) ? nothing :
             TextDisplay.engineering(primitive.fillet, :meter)
)

TextDisplay.@showfields Shell "Shell" layer -> (
    t = TextDisplay.engineering(layer.t, :meter),
)

TextDisplay.@showfields DifferenceShape "DifferenceShape" shape -> (
    outer = string(nameof(typeof(shape.outer))),
    holes = length(shape.holes)
)

TextDisplay.@showfields AssemblyShape "AssemblyShape" shape -> (
    members = length(shape.members),
)

TextDisplay.@showfields AssemblyMember "AssemblyMember" member -> (
    item = string(nameof(typeof(member.item))),
    at = _display_pose(member.at)
)

TextDisplay.name(::Type{<:EmptyBoundary}) = "Empty boundary"
Base.summary(io::IO, ::EmptyBoundary) = print(io, "Empty boundary")
Base.show(io::IO, ::EmptyBoundary) = print(io, "EmptyBoundary()")
Base.show(io::IO, ::MIME"text/plain", value::EmptyBoundary) = show(io, value)

function _ring_count(pattern::Ring)
    pattern.n isa _DeferredCardinality && return "capacity()"
    return pattern.n
end

TextDisplay.@showfields Ring "Ring" pattern -> (
    n = _ring_count(pattern),
    r = pattern.r === nothing ? nothing : TextDisplay.engineering(pattern.r, :meter),
    φ₀ = iszero(pattern.φ0) ? nothing : TextDisplay.angle(pattern.φ0),
    Δφ = isapprox(pattern.span, oftype(pattern.span, 2π)) ? nothing :
         TextDisplay.angle(pattern.span),
    gap = iszero(pattern.gap_frac) ? nothing : TextDisplay.value(pattern.gap_frac)
)

TextDisplay.@showfields Hexa "Hexagonal course" pattern -> (
    course = pattern.course,
    n = 6pattern.course,
    φ₀ = iszero(pattern.φ0) ? nothing : TextDisplay.angle(pattern.φ0),
    gap = iszero(pattern.gap_frac) ? nothing : TextDisplay.value(pattern.gap_frac)
)

TextDisplay.@showfields Polar "Polar" pattern -> (
    nᵣ = pattern.nr,
    nφ = pattern.nφ,
    r₀ = TextDisplay.engineering(pattern.r0, :meter),
    Δr = TextDisplay.engineering(pattern.dr, :meter),
    φ₀ = iszero(pattern.φ0) ? nothing : TextDisplay.angle(pattern.φ0),
    Δφ = isapprox(pattern.span, oftype(pattern.span, 2π)) ? nothing :
         TextDisplay.angle(pattern.span)
)

TextDisplay.@showfields Fill "Fill" pattern -> (
    r = TextDisplay.engineering(pattern.r, :meter),
    φ = iszero(pattern.φ) ? nothing : TextDisplay.angle(pattern.φ),
    φ₀ = iszero(pattern.φ0) ? nothing : TextDisplay.angle(pattern.φ0),
    Δφ = isapprox(pattern.span, oftype(pattern.span, 2π)) ? nothing :
         TextDisplay.angle(pattern.span)
)

TextDisplay.@showfields Lattice "Lattice" pattern -> (
    nₓ = pattern.nx,
    nᵧ = pattern.ny,
    Δx = TextDisplay.engineering(pattern.dx, :meter),
    Δy = TextDisplay.engineering(pattern.dy, :meter)
)

TextDisplay.@showfields DiameterFactor "DiameterFactor" factor -> (
    k = TextDisplay.value(factor.k),
)

TextDisplay.@showfields FillFactor "FillFactor" factor -> (
    η = TextDisplay.value(factor.η),
)

TextDisplay.@showfields TabulatedCompaction "TabulatedCompaction" compaction -> (
    data = _bounded_collection(compaction.data),
)

TextDisplay.@showfields AffineCompaction "AffineCompaction" compaction -> (
    map = _bounded_collection(compaction.map),
)

TextDisplay.@showfields LayRatio "LayRatio" lay -> (
    q = TextDisplay.value(lay.q),
)

TextDisplay.@showfields Pitch "Pitch" lay -> (
    p = TextDisplay.engineering(lay.p, :meter),
)

TextDisplay.@showfields LayAngle "LayAngle" lay -> (
    α = TextDisplay.angle(lay.α),
)

TextDisplay.@showfields Helix{<:LayRatio} "Helix" path -> (
    q = TextDisplay.value(path.lay.q),
    dir = path.dir > 0 ? "+1" : "−1",
    φ₀ = iszero(path.φ0) ? nothing : TextDisplay.angle(path.φ0)
)

TextDisplay.@showfields Helix{<:Pitch} "Helix" path -> (
    p = TextDisplay.engineering(path.lay.p, :meter),
    dir = path.dir > 0 ? "+1" : "−1",
    φ₀ = iszero(path.φ0) ? nothing : TextDisplay.angle(path.φ0)
)

TextDisplay.@showfields Helix{<:LayAngle} "Helix" path -> (
    α = TextDisplay.angle(path.lay.α),
    dir = path.dir > 0 ? "+1" : "−1",
    φ₀ = iszero(path.φ0) ? nothing : TextDisplay.angle(path.φ0)
)

TextDisplay.name(::Type{<:Helix}) = "Helix"
Base.summary(io::IO, ::Helix) = print(io, "Helix")
function Base.show(io::IO, path::Helix)
    return TextDisplay.fields(io,
        "Helix",
        (
            lay = sprint(show, path.lay; context = :compact => true),
            dir = path.dir > 0 ? "+1" : "−1",
            φ₀ = iszero(path.φ0) ? nothing : TextDisplay.angle(path.φ0)
        ))
end
function Base.show(io::IO, ::MIME"text/plain", path::Helix)
    get(io, :compact, false) && return show(io, path)
    return TextDisplay.fields(io,
        "Helix",
        (
            lay = sprint(show, path.lay; context = :compact => true),
            dir = path.dir > 0 ? "+1" : "−1",
            φ₀ = iszero(path.φ0) ? nothing : TextDisplay.angle(path.φ0)
        );
        multiline = true)
end

_datasheet_unit(::Val{:U0}) = "kV"
_datasheet_unit(::Val{:U}) = "kV"
_datasheet_unit(::Val{:conductor_cross_section}) = "mm²"
_datasheet_unit(::Val{:screen_cross_section}) = "mm²"
_datasheet_unit(::Val{:armor_cross_section}) = "mm²"
_datasheet_unit(::Val{:resistance}) = "Ω/km"
_datasheet_unit(::Val{:capacitance}) = "μF/km"
_datasheet_unit(::Val{:inductance}) = "mH/km"
_datasheet_unit(::Val) = nothing

function _datasheet_value(name::Symbol, item)
    item === nothing && return nothing
    text = item isa Real ? TextDisplay.value(item) :
           sprint(show, item; context = :compact => true)
    unit = _datasheet_unit(Val(name))
    return unit === nothing ? text : "$text $unit"
end

function _datasheet_fields(info::DatasheetInfo)
    names = keys(info)
    values = Tuple(_datasheet_value(name, info[name]) for name in names)
    return NamedTuple{names}(values)
end

TextDisplay.@showfields DatasheetInfo "DatasheetInfo" info -> _datasheet_fields(info)

_compact_text(value) = sprint(show, value; context = :compact => true)

function _part_label(region::Region)
    primitive = _compact_text(region.primitive)
    return "Region :$(region.tag) · $primitive · $(region.material.kind)"
end

function _part_label(stack::Stack)
    count = length(stack.items)
    return "Stack · $count $(count == 1 ? "part" : "parts")"
end

function _part_label(group::Group)
    attributes = String[]
    !_identity_pose(group.at) && push!(attributes, _compact_text(group.at))
    one_member = group.pattern === nothing ||
                 group.pattern isa Ring && group.pattern.n == 1 &&
                 (group.pattern.r === nothing || iszero(group.pattern.r))
    push!(attributes, one_member ? "one member" : _compact_text(group.pattern))
    group.path === nothing || push!(attributes, _compact_text(group.path))
    group.compact === nothing || push!(attributes, _compact_text(group.compact))
    return "Group :$(group.name) · $(join(attributes, " · "))"
end

function _part_label(assembly::Assembly)
    attributes = String[]
    !_identity_pose(assembly.at) && push!(attributes, _compact_text(assembly.at))
    if assembly.item isa Tuple
        count = length(assembly.item)
        push!(attributes, "$count explicit $(count == 1 ? "member" : "members")")
    else
        push!(attributes, _compact_text(assembly.pattern))
        assembly.path === nothing || push!(attributes, _compact_text(assembly.path))
        assembly.compact === nothing || push!(attributes, _compact_text(assembly.compact))
    end
    return "Assembly · $(join(attributes, " · "))"
end

function _part_label(enclosure::Enclosure)
    attributes = String[_compact_text(enclosure.primitive)]
    !_identity_pose(enclosure.at) && push!(attributes, _compact_text(enclosure.at))
    return "Enclosure :$(enclosure.tag) · $(join(attributes, " · "))"
end

_part_children(::Region) = ()
_part_children(stack::Stack) = Tuple(_part_node(item) for item in stack.items)
_part_children(group::Group) = (_part_node(group.item),)
function _part_children(assembly::Assembly)
    if assembly.item isa Tuple
        return Tuple((
                         label = _identity_pose(member.at) ? _part_label(member.item) :
                                 string(_part_label(member.item), " · ", _compact_text(member.at)),
                         children = _part_children(member.item),
                         noun = "parts"
                     ) for member in assembly.item)
    end
    return (_part_node(assembly.item),)
end
function _part_children(enclosure::Enclosure)
    fill = enclosure.fill isa Material ?
           "fill · Material · $(enclosure.fill.kind)" :
           string("fill · ", _part_label(enclosure.fill))
    children = Any[_part_node(enclosure.item), (
        label = fill, noun = "parts")]
    enclosure.wall === nothing || push!(children,
        (
            label = string("wall · ", _part_label(enclosure.wall)),
            children = _part_children(enclosure.wall),
            noun = "parts"
        ))
    return Tuple(children)
end

function _part_node(part::AbstractCablePart)
    (
        label = _part_label(part),
        children = _part_children(part),
        noun = "parts"
    )
end

TextDisplay.name(::Type{<:Region}) = "Region"
Base.summary(io::IO, region::Region) = print(io, "Region :", region.tag)
Base.show(io::IO, region::Region) = print(io, _part_label(region))
Base.show(io::IO, ::MIME"text/plain", region::Region) = show(io, region)

TextDisplay.name(::Type{<:Stack}) = "Stack"
function Base.summary(io::IO, stack::Stack)
    count = length(stack.items)
    print(io, "Stack with $count $(count == 1 ? "part" : "parts")")
end
function Base.show(io::IO, stack::Stack)
    count = length(stack.items)
    print(io, "Stack($count $(count == 1 ? "part" : "parts"))")
end
function Base.show(io::IO, ::MIME"text/plain", stack::Stack)
    get(io, :compact, false) && return show(io, stack)
    return TextDisplay.tree(io, _part_label(stack), _part_children(stack); noun = "parts")
end

TextDisplay.name(::Type{<:Group}) = "Group"
Base.summary(io::IO, group::Group) = print(io, "Group :", group.name)
function Base.show(io::IO, group::Group)
    print(io, "Group(:", group.name, "; ")
    if group.pattern === nothing
        print(io, "one member")
    else
        show(IOContext(io, :compact => true), group.pattern)
    end
    group.path === nothing ||
        (print(io, ", path="); show(IOContext(io, :compact => true), group.path))
    group.compact === nothing ||
        (print(io, ", compact="); show(IOContext(io, :compact => true), group.compact))
    print(io, ")")
end
function Base.show(io::IO, ::MIME"text/plain", group::Group)
    get(io, :compact, false) && return show(io, group)
    return TextDisplay.tree(io, _part_label(group), _part_children(group); noun = "parts")
end

TextDisplay.name(::Type{<:Assembly}) = "Assembly"
function Base.summary(io::IO, assembly::Assembly)
    assembly.item isa Tuple ?
    print(io, "Assembly with $(length(assembly.item)) members") :
    print(io, "Assembly")
end
function Base.show(io::IO, assembly::Assembly)
    if assembly.item isa Tuple
        print(io, "Assembly($(length(assembly.item)) members)")
    else
        print(io, "Assembly(")
        show(IOContext(io, :compact => true), assembly.pattern)
        print(io, ")")
    end
end
function Base.show(io::IO, ::MIME"text/plain", assembly::Assembly)
    get(io, :compact, false) && return show(io, assembly)
    return TextDisplay.tree(io, _part_label(assembly), _part_children(assembly); noun = "members")
end

TextDisplay.name(::Type{<:Enclosure}) = "Enclosure"
Base.summary(io::IO, enclosure::Enclosure) = print(io, "Enclosure :", enclosure.tag)
function Base.show(io::IO, enclosure::Enclosure)
    print(io, "Enclosure(:", enclosure.tag, "; ")
    show(IOContext(io, :compact => true), enclosure.primitive)
    print(io, ")")
end
function Base.show(io::IO, ::MIME"text/plain", enclosure::Enclosure)
    get(io, :compact, false) && return show(io, enclosure)
    return TextDisplay.tree(
        io, _part_label(enclosure), _part_children(enclosure); noun = "parts"
    )
end

TextDisplay.name(::Type{<:PlacedRegion}) = "PlacedRegion"
Base.summary(io::IO, region::PlacedRegion) = print(io, "Placed region :", region.source.tag)
function Base.show(io::IO, region::PlacedRegion)
    terminal = region.terminal === nothing ? "" : "; terminal=:$(region.terminal)"
    print(io, "PlacedRegion(:", region.source.tag, terminal, "; ")
    show(IOContext(io, :compact => true), region.primitive)
    print(io, ")")
end
function Base.show(io::IO, ::MIME"text/plain", region::PlacedRegion)
    get(io, :compact, false) && return show(io, region)
    attributes = String[
    "source    $(_part_label(region.source))",
    "primitive $(_compact_text(region.primitive))"
]
    region.terminal === nothing || push!(attributes, "terminal  :$(region.terminal)")
    isempty(region.placement.patterns) || push!(attributes,
        "patterns  $(length(region.placement.patterns))")
    isempty(region.paths) || push!(attributes, "paths     $(length(region.paths))")
    return TextDisplay.tree(io, "PlacedRegion :$(region.source.tag)", attributes)
end

TextDisplay.name(::Type{<:CableGeometry}) = "CableGeometry"
function Base.summary(io::IO, geometry::CableGeometry)
    print(io, "Cable geometry with $(length(geometry.regions)) regions")
end
function Base.show(io::IO, geometry::CableGeometry)
    print(io, "CableGeometry(regions=$(length(geometry.regions)); outer=")
    show(IOContext(io, :compact => true), geometry.outer)
    print(io, ")")
end
function Base.show(io::IO, ::MIME"text/plain", geometry::CableGeometry)
    get(io, :compact, false) && return show(io, geometry)
    children = Any[(
        label = "outer    $(_compact_text(geometry.outer))", noun = "regions")]
    append!(children,
        (
            label = string(
                "region :", region.source.tag,
                region.terminal === nothing ? "" : " · terminal :$(region.terminal)",
                " · ", _compact_text(region.primitive)
            ),
            noun = "regions"
        ) for region in geometry.regions)
    count = length(geometry.regions)
    return TextDisplay.tree(
        io,
        "CableGeometry · $count $(count == 1 ? "region" : "regions")",
        Tuple(children);
        noun = "regions"
    )
end

TextDisplay.name(::Type{<:CableDesign}) = "CableDesign"
function Base.summary(io::IO, design::CableDesign)
    print(io, "CableDesign \"", design.cable_id, "\"")
end
function Base.show(io::IO, design::CableDesign)
    print(
        io,
        "CableDesign(\"", design.cable_id, "\"; terminals=",
        length(design.terminal_order), ", regions=", length(design.geometry.regions), ")"
    )
end
function Base.show(io::IO, ::MIME"text/plain", design::CableDesign)
    get(io, :compact, false) && return show(io, design)
    terminals = isempty(design.terminal_order) ? "none" : join(design.terminal_order, ", ")
    root = _part_node(design.root)
    children = (
        (label = "terminals  $terminals", noun = "parts"),
        (label = "regions    $(length(design.geometry.regions))", noun = "parts"),
        (label = "diameter   $(TextDisplay.engineering(2outer_radius(design), :meter))",
            noun = "parts"),
        (label = "root       $(root.label)", children = root.children, noun = "parts")
    )
    return TextDisplay.tree(io, "CableDesign \"$(design.cable_id)\"", children)
end

TextDisplay.name(::Type{<:LineCableSystem}) = "LineCableSystem"
function Base.summary(io::IO, system::LineCableSystem)
    print(io, "LineCableSystem \"", system.system_id, "\"")
end
function Base.show(io::IO, system::LineCableSystem)
    print(
        io,
        "LineCableSystem(\"", system.system_id, "\"; cables=", ncables(system),
        ", terminals=", length(system.terminal_order), ")"
    )
end
function Base.show(io::IO, ::MIME"text/plain", system::LineCableSystem)
    get(io, :compact, false) && return show(io, system)
    cables = Tuple((
                       label = string(
                           design.cable_id,
                           " · ", _compact_text(position),
                           " · ", join(
                               ("$terminal→$phase"
                               for (terminal, phase) in
                                   zip(design.terminal_order, connections)),
                               ", "
                           )
                       ),
                       noun = "cables"
                   )
    for (design, position, connections) in zip(
        system.designs, system.positions, system.connections
    ))
    children = (
        (label = "length     $(TextDisplay.engineering(system.line_length, :meter))",
            noun = "cables"),
        (label = "terminals  $(length(system.terminal_order))", noun = "cables"),
        (label = "cables", children = cables, noun = "cables")
    )
    return TextDisplay.tree(io, "LineCableSystem \"$(system.system_id)\"", children)
end

TextDisplay.name(::Type{<:CablesLibrary}) = "CablesLibrary"
function Base.summary(io::IO, library::CablesLibrary)
    count = length(library)
    print(io, "CablesLibrary with $count $(count == 1 ? "design" : "designs")")
end
function Base.show(io::IO, library::CablesLibrary)
    count = length(library)
    print(io, "CablesLibrary($count $(count == 1 ? "design" : "designs"))")
end
function Base.show(io::IO, ::MIME"text/plain", library::CablesLibrary)
    get(io, :compact, false) && return show(io, library)
    ids = sort!(collect(keys(library)))
    children = Tuple(begin
                         design = library[cable_id]
                         (
                             label = string(
                                 cable_id, " · ", length(design.terminal_order), " terminals · ",
                                 length(design.geometry.regions), " regions · D=",
                                 TextDisplay.engineering(2outer_radius(design), :meter)
                             ),
                             noun = "designs"
                         )
                     end
    for cable_id in ids)
    count = length(library)
    return TextDisplay.tree(
        io,
        "CablesLibrary · $count $(count == 1 ? "design" : "designs")",
        children;
        noun = "designs"
    )
end

TextDisplay.name(::Type{<:RoundedSectorShape}) = "Rounded sector shape"
Base.summary(io::IO, ::RoundedSectorShape) = print(io, "Resolved rounded sector")
function Base.show(io::IO, shape::RoundedSectorShape)
    print(io, "ResolvedRoundedSector(")
    show(IOContext(io, :compact => true), shape.primitive)
    _identity_pose(shape.at) ||
        (print(io, "; at="); show(IOContext(io, :compact => true), shape.at))
    print(io, ")")
end
Base.show(io::IO, ::MIME"text/plain", shape::RoundedSectorShape) = show(io, shape)

TextDisplay.name(::Type{<:ShellShape}) = "Shell shape"
Base.summary(io::IO, ::ShellShape) = print(io, "Resolved shell")
function Base.show(io::IO, shape::ShellShape)
    print(io, "ResolvedShell(inner=")
    show(IOContext(io, :compact => true), shape.inner)
    print(io, ", outer=")
    show(IOContext(io, :compact => true), shape.outer)
    print(io, ")")
end
function Base.show(io::IO, ::MIME"text/plain", shape::ShellShape)
    get(io, :compact, false) && return show(io, shape)
    return TextDisplay.tree(io,
        "Resolved shell",
        (
            (label = "inner  $(_compact_text(shape.inner))", noun = "boundaries"),
            (label = "outer  $(_compact_text(shape.outer))", noun = "boundaries")
        ))
end

function _preview_geometry_count(geometry)
    points = geometry isa GeometryBasics.Polygon ? geometry.exterior : geometry
    applicable(length, points) && return length(points)
    return count(_ -> true, points)
end

TextDisplay.name(::Type{<:PreviewShape}) = "Preview shape"
function Base.summary(io::IO, shape::PreviewShape)
    print(io, "Preview shape :", shape.tag)
end
function Base.show(io::IO, shape::PreviewShape)
    print(io, "PreviewShape(:", shape.tag, "; vertices=",
        _preview_geometry_count(shape.geometry), ")")
end
Base.show(io::IO, ::MIME"text/plain", shape::PreviewShape) = show(io, shape)
