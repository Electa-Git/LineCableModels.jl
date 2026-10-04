function Region(tag, primitive, material; combine::Symbol = :product)
    return parameterize(
        DataModel.Region, DataModel.Region, (tag, primitive, material); combine
    )
end

# `DataModel.Region(::Symbol, primitive, material)` is deliberately permissive
# for scalar declaration types, so explicit finite inputs need narrower methods
# to use the collection constructor instead of that scalar constructor.
const _FiniteRegionInput = Union{AbstractGrid, Gridspace}

function Region(tag::_FiniteRegionInput, primitive, material; combine::Symbol = :product)
    parameterize(DataModel.Region, DataModel.Region, (tag, primitive, material); combine)
end
function Region(tag, primitive::_FiniteRegionInput, material; combine::Symbol = :product)
    parameterize(DataModel.Region, DataModel.Region, (tag, primitive, material); combine)
end
function Region(tag, primitive, material::_FiniteRegionInput; combine::Symbol = :product)
    parameterize(DataModel.Region, DataModel.Region, (tag, primitive, material); combine)
end
function Region(
        tag::Symbol,
        primitive::_FiniteRegionInput,
        material;
        combine::Symbol = :product
)
    parameterize(DataModel.Region, DataModel.Region, (tag, primitive, material); combine)
end
function Region(
        tag::Symbol,
        primitive,
        material::_FiniteRegionInput;
        combine::Symbol = :product
)
    parameterize(DataModel.Region, DataModel.Region, (tag, primitive, material); combine)
end
function Region(
        tag::Symbol,
        primitive::_FiniteRegionInput,
        material::_FiniteRegionInput;
        combine::Symbol = :product
)
    parameterize(DataModel.Region, DataModel.Region, (tag, primitive, material); combine)
end
function Region(
        tag::_FiniteRegionInput,
        primitive::_FiniteRegionInput,
        material;
        combine::Symbol = :product
)
    parameterize(DataModel.Region, DataModel.Region, (tag, primitive, material); combine)
end
function Region(
        tag::_FiniteRegionInput,
        primitive,
        material::_FiniteRegionInput;
        combine::Symbol = :product
)
    parameterize(DataModel.Region, DataModel.Region, (tag, primitive, material); combine)
end
function Region(
        tag,
        primitive::_FiniteRegionInput,
        material::_FiniteRegionInput;
        combine::Symbol = :product
)
    parameterize(DataModel.Region, DataModel.Region, (tag, primitive, material); combine)
end
function Region(
        tag::_FiniteRegionInput,
        primitive::_FiniteRegionInput,
        material::_FiniteRegionInput;
        combine::Symbol = :product
)
    parameterize(DataModel.Region, DataModel.Region, (tag, primitive, material); combine)
end

function Stack(items...; combine::Symbol = :product)
    isempty(items) && throw(ArgumentError("layers require at least one part"))
    return parameterize(DataModel.Stack, DataModel.Stack, items; combine)
end

"""
$(TYPEDSIGNATURES)

Compose physical parts in outward order.

# Arguments

- `parts`: physical declarations ordered from the center outward.

# Keywords

- `combine=:product`: gridspace composition rule.

# Returns

- A `Stack`, or a `Gridspace{Stack}` when a direct argument varies.
"""
layers(parts...; combine::Symbol = :product) = Stack(parts...; combine)

function Group(
        name,
        item;
        at = DataModel.Pose2(0, 0, 0),
        pattern = nothing,
        path = nothing,
        compact = nothing,
        boundary = nothing,
        combine::Symbol = :product
)
    values = (name, at, item, pattern, path, compact, boundary)
    return parameterize(DataModel.Group, DataModel.Group, values; combine)
end

function Assembly(
        item;
        at = DataModel.Pose2(0, 0, 0),
        pattern,
        names = nothing,
        path = nothing,
        compact = nothing,
        combine::Symbol = :product
)
    values = (at, item, pattern, path, compact, names)
    return parameterize(DataModel.Assembly, DataModel.Assembly, values; combine)
end

function _explicit_assembly(members...)
    placed = map(members) do member
        member isa DataModel.AssemblyMember ? member :
        member isa DataModel.AbstractCablePart ? DataModel.AssemblyMember(member) :
        throw(ArgumentError("assembly members must be physical cable parts"))
    end
    return DataModel.Assembly(
        DataModel.Pose2(0, 0, 0), Tuple(placed), nothing, nothing, nothing, nothing
    )
end

"""
$(TYPEDSIGNATURES)

Preserve explicit physical members and their independent terminal identities.

With `pattern`, retain one prototype and one placement pattern. Without a
pattern, retain one or more heterogeneous members and their local poses.

# Arguments

- `members`: physical members, optionally placed with [`at`](@ref).

# Keywords

- `pattern=nothing`: placement pattern for one repeated prototype, or explicit members.
- `names=nothing`: exact terminal names for repeated terminal-bearing members.
- `path=nothing`: shared longitudinal path declaration.
- `compact=nothing`: explicit compaction law.
- `combine=:product`: gridspace composition rule.

# Returns

- An `Assembly`, or a `Gridspace{Assembly}` when a direct argument varies.
"""
function assembly(
        members...;
        pattern = nothing,
        names = nothing,
        path = nothing,
        compact = nothing,
        combine::Symbol = :product
)
    isempty(members) && throw(ArgumentError("assembly requires at least one member"))
    if pattern === nothing
        names === nothing && path === nothing && compact === nothing ||
            throw(ArgumentError(
                "explicit assembly members own their names, paths and compaction; " *
                "provide pattern to repeat one prototype"))
        return parameterize(DataModel.Assembly, _explicit_assembly, members; combine)
    end
    length(members) == 1 || throw(ArgumentError(
        "a repeated assembly requires exactly one prototype with pattern"))
    return Assembly(
        first(members);
        pattern,
        names,
        path,
        compact,
        combine
    )
end

function Enclosure(
        tag,
        item;
        at = DataModel.Pose2(0, 0, 0),
        primitive,
        fill,
        wall = nothing,
        combine::Symbol = :product
)
    values = (tag, at, primitive, item, fill, wall)
    return parameterize(DataModel.Enclosure, DataModel.Enclosure, values; combine)
end

_terminal_eligible(region::DataModel.Region) = region.material.kind === :conductor ? 1 : 0
_terminal_eligible(stack::DataModel.Stack) = sum(_terminal_eligible, stack.items; init = 0)
_terminal_eligible(group::DataModel.Group) = _terminal_eligible(group.item)
function _terminal_eligible(
        assembly::DataModel.Assembly{<:Any, <:DataModel.AbstractCablePart}
)
    _terminal_eligible(assembly.item)
end
function _terminal_eligible(assembly::DataModel.Assembly{<:Any, <:Tuple})
    return sum(member -> _terminal_eligible(member.item), assembly.item; init = 0)
end
function _terminal_eligible(enclosure::DataModel.Enclosure)
    count = _terminal_eligible(enclosure.item)
    enclosure.fill isa DataModel.Region && (count += _terminal_eligible(enclosure.fill))
    enclosure.wall === nothing || (count += _terminal_eligible(enclosure.wall))
    return count
end

function _terminal(name, parts...)
    root = DataModel.Stack(parts...)
    _terminal_eligible(root) > 0 || throw(ArgumentError(
        "terminal :$name requires a conductive descendant"
    ))
    return DataModel.Group(
        name, DataModel.Pose2(0, 0, 0), root, nothing, nothing, nothing
    )
end

"""
$(TYPEDSIGNATURES)

Coalesce every conductive descendant of an ordered physical subtree into one
retained terminal.

# Arguments

- `name`: retained electrical terminal name.
- `parts`: physical declarations ordered from the center outward.

# Keywords

- `combine=:product`: gridspace composition rule.

# Returns

- A terminal-owning `Group`, or a `Gridspace{Group}` when a direct argument
  varies.

# Errors

- Throws `ArgumentError` when the realized subtree contains no conductive
  descendant.
"""
function terminal(name, parts...; combine::Symbol = :product)
    isempty(parts) && throw(ArgumentError("terminal requires at least one part"))
    return parameterize(DataModel.Group, _terminal, (name, parts...); combine)
end

"""
$(TYPEDSIGNATURES)

Bind an intrinsic primitive definition to one material and local physical tag.

# Arguments

- `material`: constitutive material record.
- `primitive`: intrinsic cross-sectional primitive definition.

# Keywords

- `tag=:solid`: local physical identity.
- `combine=:product`: gridspace composition rule.

# Returns

- A `Region`, or a `Gridspace{Region}` when a direct argument varies.
"""
function solid(material, primitive; tag = :solid, combine::Symbol = :product)
    Region(tag, primitive, material; combine)
end

"""
$(TYPEDSIGNATURES)

Declare one outward material layer of thickness `t` \\[m\\].

# Arguments

- `material`: constitutive material record.

# Keywords

- `t`: normal layer thickness \\[m\\].
- `tag=:shell`: local physical identity.
- `combine=:product`: gridspace composition rule.

# Returns

- A contextual `Region`, or a `Gridspace{Region}` when a direct argument
  varies.
"""
function shell(material; t, tag = :shell, combine::Symbol = :product)
    Region(tag, Shell(t), material; combine)
end

"""
$(TYPEDSIGNATURES)

Declare a conductive core region.

# Arguments

- `material`: material with `kind == :conductor`.
- `primitive`: intrinsic core geometry. The keyword form constructs a disk of
  radius `r` \\[m\\].

# Keywords

- `r`: disk radius \\[m\\].
- `tag=:core`: local physical identity.
- `combine=:product`: gridspace composition rule.

# Returns

- A conductive `Region`, or a `Gridspace{Region}` when a direct argument
  varies.
"""
function core(material, primitive; tag = :core, combine::Symbol = :product)
    caller = (resolved_material,
        resolved_primitive,
        resolved_tag) -> DataModel.Region(
        resolved_tag,
        resolved_primitive,
        validate(resolved_material, DataModel.Region, :core, (:conductor,))
    )
    return parameterize(
        DataModel.Region, caller, (material, primitive, tag); combine
    )
end

function core(material; r, tag = :core, combine::Symbol = :product)
    core(material, Disk(r); tag, combine)
end

function _role_shell(role::Symbol, allowed, material, t, tag)
    return DataModel.Region(
        tag,
        DataModel.Shell(t),
        validate(material, DataModel.Region, role, allowed)
    )
end

function _shell_role(
        role::Symbol, allowed::Tuple, material, t, tag;
        combine::Symbol
)
    caller = (resolved_material,
        resolved_t,
        resolved_tag) -> _role_shell(
        role, allowed, resolved_material, resolved_t, resolved_tag)
    return parameterize(DataModel.Region, caller, (material, t, tag); combine)
end

"""Declare an insulating layer of thickness `t` \\[m\\]."""
function insulation(material; t, tag = :insulation, combine::Symbol = :product)
    _shell_role(:insulation, (:insulator,), material, t, tag; combine)
end
"""Declare a semiconductive or conductive screen layer of thickness `t` \\[m\\]."""
function screen(material; t, tag = :screen, combine::Symbol = :product)
    _shell_role(:screen, (:semicon, :conductor), material, t, tag; combine)
end
"""Declare a conductive sheath layer of thickness `t` \\[m\\]."""
function sheath(material; t, tag = :sheath, combine::Symbol = :product)
    _shell_role(:sheath, (:conductor,), material, t, tag; combine)
end
"""Declare a nonconducting bedding layer of thickness `t` \\[m\\]."""
function bedding(material; t, tag = :bedding, combine::Symbol = :product)
    _shell_role(:bedding, (:insulator,), material, t, tag; combine)
end
"""Declare an insulating jacket layer of thickness `t` \\[m\\]."""
function jacket(material; t, tag = :jacket, combine::Symbol = :product)
    _shell_role(:jacket, (:insulator,), material, t, tag; combine)
end

"""
$(TYPEDSIGNATURES)

Bind a nonconducting filler material to an intrinsic primitive definition.

# Arguments

- `material`: material with `kind == :insulator`.
- `primitive`: intrinsic filler geometry.

# Keywords

- `tag=:filler`: local physical identity.
- `combine=:product`: gridspace composition rule.

# Returns

- A filler `Region`, or a `Gridspace{Region}` when a direct argument varies.
"""
function filler(material, primitive; tag = :filler, combine::Symbol = :product)
    caller = (resolved_material,
        resolved_primitive,
        resolved_tag) -> DataModel.Region(
        resolved_tag,
        resolved_primitive,
        validate(resolved_material, DataModel.Region, :filler, (:insulator,))
    )
    return parameterize(
        DataModel.Region, caller, (material, primitive, tag); combine
    )
end

_path(::Nothing, dir, φ0) = nothing
_path(path::DataModel.Helix, dir, φ0) = path
_path(lay, dir, φ0) = DataModel.Helix(lay; dir, φ0)

function _wires(
        material, shape, pattern, path, compact, tag
)
    source = DataModel.Region(
        tag,
        shape,
        validate(material, DataModel.Region, :wires, (:conductor,))
    )
    return DataModel.Group(
        tag, DataModel.Pose2(0, 0, 0), source, pattern, path, compact
    )
end

"""
$(TYPEDSIGNATURES)

Declare one repeated course of conductive members.

Supply `pattern` and `path` directly for general placement, or supply `n`, `r`,
and `lay` for the practical ring-course form.

# Arguments

- `material`: material with `kind == :conductor`.

# Keywords

- `shape`: intrinsic member primitive.
- `pattern=nothing`: explicit member placement pattern.
- `path=nothing`: explicit longitudinal path.
- `n=nothing`: exact ring cardinality or `capacity()`.
- `r=nothing`: member-center ring radius \\[m\\].
- `gap_frac=0`: fractional adjacent clearance \\[dimensionless\\].
- `lay=nothing`: one lay law used to construct a `Helix`.
- `dir=1`: helix handedness, `1` or `-1` \\[dimensionless\\].
- `φ0=0`: initial angular position \\[rad\\].
- `compact=nothing`: explicit compaction law.
- `tag=:wire`: local member and group identity.
- `combine=:product`: gridspace composition rule.

# Returns

- A repeated-member `Group`, or a `Gridspace{Group}` when a direct argument
  varies.
"""
function wires(
        material;
        shape,
        pattern = nothing,
        path = nothing,
        n = nothing,
        r = nothing,
        gap_frac = 0,
        lay = nothing,
        dir = 1,
        φ0 = 0,
        compact = nothing,
        tag = :wire,
        combine::Symbol = :product
)
    caller = function (
            resolved_material, resolved_shape, resolved_pattern, resolved_path,
            resolved_n, resolved_r, resolved_gap, resolved_lay, resolved_dir,
            resolved_φ0, resolved_compact, resolved_tag
    )
        if resolved_pattern !== nothing
            resolved_n === nothing && resolved_r === nothing || throw(ArgumentError(
                "pattern cannot be combined with n or r"
            ))
            resolved_lay === nothing || throw(ArgumentError(
                "pattern-oriented wires use path rather than lay"
            ))
            return _wires(
                resolved_material,
                resolved_shape,
                resolved_pattern,
                resolved_path,
                resolved_compact,
                resolved_tag
            )
        end
        resolved_n === nothing && throw(ArgumentError(
            "ring-course wires require n"
        ))
        resolved_r === nothing && throw(ArgumentError(
            "ring-course wires require r"
        ))
        resolved_path === nothing || throw(ArgumentError(
            "ring-course wires use lay rather than path"
        ))
        return _wires(
            resolved_material,
            resolved_shape,
            DataModel.Ring(
                resolved_n;
                r = resolved_r,
                φ0 = resolved_φ0,
                gap_frac = resolved_gap
            ),
            _path(resolved_lay, resolved_dir, resolved_φ0),
            resolved_compact,
            resolved_tag
        )
    end
    values = (
        material, shape, pattern, path, n, r, gap_frac, lay, dir, φ0,
        compact, tag
    )
    return parameterize(DataModel.Group, caller, values; combine)
end

function _course_schedule(value, count::Int, name::Symbol)
    if value isa Union{Tuple, AbstractVector}
        length(value) == count || throw(DimensionMismatch(
            "$name requires exactly $count course values"
        ))
        return Tuple(value)
    end
    return ntuple(_ -> value, count)
end

function _count_schedule(value, count::Int)
    if value isa Union{Tuple, AbstractVector}
        length(value) == count || throw(DimensionMismatch(
            "n requires exactly $count course values"
        ))
        return Tuple(value)
    end
    value isa Integer && !(value isa Bool) && value > 0 &&
        return ntuple(index -> index * Int(value), count)
    value === DataModel.capacity() && return ntuple(_ -> value, count)
    throw(ArgumentError(
        "n must be a positive base count, exact course schedule, or capacity()"
    ))
end

function _stranding_paths(lay, dir, φ0)
    lay === nothing && return nothing
    if lay isa Union{Tuple, AbstractVector}
        count = length(lay)
        count > 0 || throw(ArgumentError("lay schedule cannot be empty"))
        directions = _course_schedule(dir, count, :dir)
        angles = _course_schedule(φ0, count, :φ0)
        return ntuple(index -> _path(lay[index], directions[index], angles[index]), count)
    end
    dir isa Union{Tuple, AbstractVector} && throw(ArgumentError(
        "a dir schedule requires a lay schedule of the same length"
    ))
    φ0 isa Union{Tuple, AbstractVector} && throw(ArgumentError(
        "a φ0 schedule requires a lay schedule of the same length"
    ))
    return _path(lay, dir, φ0)
end

function _stranded(
        material,
        center,
        shape,
        lay,
        dir,
        φ0,
        compact,
        prescribed_boundary,
        fill
)
    material = validate(material, DataModel.Region, :stranded, (:conductor,))
    fill = validate(fill, DataModel.Region, :stranded_fill, (:insulator, :semicon))
    prescribed_boundary isa Union{DataModel.Disk, DataModel.Sector} ||
        throw(ArgumentError(
            "stranded requires a nonhollow Disk or Sector boundary"
        ))
    shape isa Union{DataModel.Disk, DataModel.Rectangle} || throw(ArgumentError(
        "stranded members must be Disk wires or Rectangle strands"
    ))
    if prescribed_boundary isa DataModel.Sector
        shape isa DataModel.Disk || throw(ArgumentError(
            "a sector stranded core requires circular Disk wires"
        ))
        center === nothing || throw(ArgumentError(
            "a sector bundle infers its center strand from shape; omit center"
        ))
        compact === nothing || throw(ArgumentError(
            "sector stranded cores are intrinsically compacted; omit compact"
        ))
        compaction = nothing
    else
        if shape isa DataModel.Disk
            center === nothing && (center = shape)
            compaction = compact === nothing ? false : compact
            compaction isa Bool || throw(ArgumentError(
                "circular stranded compaction is selected with compact=true or false"
            ))
        else
            center isa DataModel.Disk || throw(ArgumentError(
                "a rectangular stranded core requires center=Disk(...)"
            ))
            compact === nothing || throw(ArgumentError(
                "rectangular strands are intrinsically bent; omit compact"
            ))
            compaction = true
        end
        center isa DataModel.Disk || throw(ArgumentError(
            "a disk-bounded stranded core requires one circular center wire"
        ))
    end

    paths = _stranding_paths(lay, dir, φ0)
    parts = DataModel.AbstractCablePart[]
    center === nothing || push!(parts,
        DataModel.Group(
            :strand,
            DataModel.Pose2(0, 0, 0),
            DataModel.Region(:wire, center, material),
            nothing,
            nothing,
            nothing
        ))
    push!(parts,
        DataModel.Group(
            :strand,
            DataModel.Pose2(0, 0, 0),
            DataModel.Region(:wire, shape, material),
            nothing,
            paths,
            nothing
        ))
    bounded = DataModel.Group(
        :core,
        DataModel.Pose2(0, 0, 0),
        DataModel.Stack(parts),
        nothing,
        nothing,
        compaction,
        prescribed_boundary
    )
    shape isa DataModel.Rectangle && return bounded
    return DataModel.Enclosure(
        :stranded,
        DataModel.Pose2(0, 0, 0),
        prescribed_boundary,
        bounded,
        fill,
        nothing
    )
end

"""
$(TYPEDSIGNATURES)

Fill one authoritative core boundary with the maximum admissible inventory of
equal source strands.

Circular source bundles have one center strand and `6k` wire strands in course `k`.
A sector admits the largest complete inventory satisfying
``[1 + 3L(L+1)]\\,\\pi a^2 \\leq A_{\\mathrm{sector}}``. Its mapped sites define
a prescribed-area power diagram. Disks clipped to those cells retain each
source area. This geometric reconstruction preserves the declared copper
fill and does not model mechanical forming or plasticity.

Rectangular strands occupy complete, area-preserving annular courses, including
a full annulus when a course contains one strand. Their requested geometric boundary is
a packing limit. The resolved geometric boundary is the occupied metal disk. Subsequent
layers start there, without an automatically generated outer filler film.

# Arguments

- `material`: material with `kind == :conductor`.

# Keywords

- `center=nothing`: circular center member for a disk-bounded core. Circular
  strands default to one center wire equal to `shape`. Rectangular strands
  require an explicit `Disk`. A sector bundle infers its center strand from
  `shape` and does not admit a separate center declaration.
- `shape`: circular wire or rectangular strand primitive.
- `boundary`: nonhollow `Disk` or `Sector` core boundary. For rectangular
  strands, a `Disk` packing limit rather than an imposed finished radius.
- `fill=air`: interstitial insulating or semiconducting material, unused by
  contiguous rectangular courses. The default
  is lossless air (``\\rho=\\infty``, ``\\epsilon_r=\\mu_r=1``).
- `lay=nothing`: one lay law or a schedule matching the inferred radial courses.
- `dir=1`: helix handedness or a schedule matching `lay`.
- `φ0=0`: helix initial angle or a schedule matching `lay` \\[rad\\].
- `compact=nothing`: circular disk-bounded strands remain circular by default.
  `true` requests area-preserving deformation. Sector and rectangular
  wire stranding are intrinsically deformed and do not admit this keyword.
- `combine=:product`: gridspace composition rule.

# Returns

- A bounded `Group` for rectangular strands, or a material-complete
  `Enclosure` for circular strands. Gridded arguments return a `Gridspace`
  of the corresponding part type, or `AbstractCablePart` for mixed shapes.
"""
function stranded(
        material;
        center = nothing,
        shape,
        lay = nothing,
        dir = 1,
        φ0 = 0,
        compact = nothing,
        boundary,
        fill = Materials.Material(:insulator, Inf),
        combine::Symbol = :product
)
    caller = function (
            resolved_material, resolved_center, resolved_shape, resolved_lay,
            resolved_dir, resolved_φ0, resolved_compact, resolved_boundary,
            resolved_fill
    )
        return _stranded(
            resolved_material,
            resolved_center,
            resolved_shape,
            resolved_lay,
            resolved_dir,
            resolved_φ0,
            resolved_compact,
            resolved_boundary,
            resolved_fill
        )
    end
    values = (
        material, center, shape, lay, dir, φ0, compact, boundary, fill
    )
    shape_type = shape isa AbstractGrid ? eltype(shape) : typeof(shape)
    target = if shape_type <: DataModel.Rectangle || shape isa Gridspace{<:DataModel.Rectangle}
        DataModel.Group
    elseif shape_type <: DataModel.Disk || shape isa Gridspace{<:DataModel.Disk}
        DataModel.Enclosure
    else
        DataModel.AbstractCablePart
    end
    return parameterize(target, caller, values; combine)
end

function _inner_radius(primitive::DataModel.Disk)
    return hypot(primitive.at.x, primitive.at.y) - primitive.r
end

function _inner_radius(primitive::DataModel.Polygon)
    cosine = cos(primitive.at.φ)
    sine = sin(primitive.at.φ)
    points = map(primitive.points) do point
        return (
            primitive.at.x + cosine * point[1] - sine * point[2],
            primitive.at.y + sine * point[1] + cosine * point[2]
        )
    end
    return minimum(eachindex(points)) do index
        next = mod1(index + 1, length(points))
        first_point = points[index]
        last_point = points[next]
        edge = (
            last_point[1] - first_point[1],
            last_point[2] - first_point[2]
        )
        length_squared = edge[1]^2 + edge[2]^2
        fraction = iszero(length_squared) ? zero(length_squared) : clamp(
            -(first_point[1] * edge[1] + first_point[2] * edge[2]) /
            length_squared,
            zero(length_squared),
            one(length_squared)
        )
        hypot(
            first_point[1] + fraction * edge[1],
            first_point[2] + fraction * edge[2]
        )
    end
end

function _milliken(
        material,
        shape,
        segment,
        segments,
        lay,
        dir,
        φ0,
        fill
)
    material = validate(material, DataModel.Region, :milliken, (:conductor,))
    shape isa DataModel.Disk || throw(ArgumentError(
        "milliken segment strands must be circular Disk wires"
    ))
    segment isa DataModel.Sector || throw(ArgumentError(
        "milliken segment must be one Sector boundary"
    ))
    segments isa Integer && !(segments isa Bool) && segments >= 2 ||
        throw(ArgumentError("milliken segments must be an integer of at least two"))
    pitch = 2pi / segments
    isapprox(segment.span, pitch) || throw(DomainError(
        segment.span,
        "milliken segment span must equal its angular pitch $pitch"
    ))
    fill = validate(fill, DataModel.Region, :milliken_fill, (:insulator, :semicon))

    segment_part = _stranded(
        material,
        nothing,
        shape,
        lay,
        dir,
        φ0,
        nothing,
        segment,
        fill
    )
    _, _, segment_primitives = DataModel.bounded_members(segment_part.item)
    center_radius = minimum(_inner_radius, segment_primitives)
    center_radius > zero(center_radius) || throw(DomainError(
        center_radius,
        "milliken segment strands leave no positive radius for the center wire"
    ))
    center = DataModel.Disk(center_radius)
    members = DataModel.AssemblyMember[]
    center_terminal = DataModel.Group(
        :milliken_center,
        DataModel.Pose2(0, 0, 0),
        DataModel.Region(:wire, center, material),
        nothing,
        nothing,
        nothing
    )
    push!(members, DataModel.AssemblyMember(center_terminal))
    for index in 1:Int(segments)
        name = Symbol(:milliken_segment_, index)
        terminal = DataModel.Group(
            name,
            DataModel.Pose2(0, 0, 0),
            segment_part.item,
            nothing,
            nothing,
            nothing
        )
        push!(members, DataModel.AssemblyMember(
            terminal,
            DataModel.Pose2(0, 0, (index - 1) * pitch)
        ))
    end
    assembly = DataModel.Assembly(
        DataModel.Pose2(0, 0, 0),
        Tuple(members),
        nothing,
        nothing,
        nothing,
        nothing
    )
    core = DataModel.Group(
        :core,
        DataModel.Pose2(0, 0, 0),
        DataModel.Stack(assembly),
        nothing,
        nothing,
        nothing
    )
    return DataModel.Enclosure(
        :milliken,
        DataModel.Pose2(0, 0, 0),
        DataModel.Disk(segment.r_back),
        core,
        fill,
        nothing
    )
end

"""
$(TYPEDSIGNATURES)

Declare one Milliken conductor as a cable-center wire surrounded by equal
stranded sector segments. Each segment contains its own bundle-center strand
and complete `6k` courses. Each conductive descendant resolves to one `:core`
terminal.

# Arguments

- `material`: conductor material shared by the center and segment strands.

# Keywords

- `shape`: circular source wire repeated inside every segment.
- `segment`: authoritative geometric boundary of one sector segment. Its span must equal
  ``2\\pi/N`` for `segments=N`.
- `segments=6`: number of equal sector segments \\[dimensionless\\].
- `lay=nothing`: common strand lay law.
- `dir=1`: helix handedness \\[dimensionless\\].
- `φ0=0`: helix initial angle \\[rad\\].
- `fill=air`: interstitial material around the center and every segment strand.
- `combine=:product`: gridspace composition rule.

The center-wire radius is inferred from the resolved segment packing so that
the center wire is tangent to the innermost strand of every equal segment.

# Returns

- A material-complete `Enclosure`, or a `Gridspace{Enclosure}` when a direct
  argument varies.
"""
function milliken(
        material;
        shape,
        segment,
        segments = 6,
        lay = nothing,
        dir = 1,
        φ0 = 0,
        fill = Materials.Material(:insulator, Inf),
        combine::Symbol = :product
)
    values = (
        material, shape, segment, segments, lay, dir, φ0, fill
    )
    return parameterize(DataModel.Enclosure, _milliken, values; combine)
end

function _rope(
        item,
        course_count,
        counts,
        lays,
        directions,
        angles,
        compactions,
        gaps
)
    item isa DataModel.AbstractCablePart || throw(ArgumentError(
        "rope item must be a physical cable part"
    ))
    central = DataModel.Group(
        :rope, DataModel.Pose2(0, 0, 0), item, nothing, nothing, nothing
    )
    parts = DataModel.AbstractCablePart[central]
    for course in 1:course_count
        push!(parts,
            DataModel.Group(
                :rope,
                DataModel.Pose2(0, 0, 0),
                item,
                DataModel.Ring(
                    counts[course];
                    r = nothing,
                    φ0 = angles[course],
                    gap_frac = gaps[course]
                ),
                _path(lays[course], directions[course], angles[course]),
                compactions[course]
            ))
    end
    return DataModel.Stack(parts)
end

"""
$(TYPEDSIGNATURES)

Repeat one physical item as a central child and concentric outer courses.

# Arguments

- `item`: physical cable part repeated by the rope.

# Keywords

- `layers`: number of outer courses \\[dimensionless\\].
- `n=6`: base count, exact course schedule or deferred maximum count `capacity()`.
- `lay=nothing`: one lay law or one law per outer course. Homogeneous schedules
  may be declared as `LayRatio(q...)`, `Pitch(p...)` or `LayAngle(α...)`.
- `dir=1`: one handedness or one value per outer course.
- `φ0=0`: one initial angle or one value per outer course \\[rad\\].
- `compact=nothing`: one compaction law or one law per outer course.
  Homogeneous scalar schedules may be declared as `FillFactor(η...)`.
- `gap_frac=0`: one clearance fraction or one value per outer course
  \\[dimensionless\\].
- `combine=:product`: gridspace composition rule.

# Returns

- A nested `Stack`, or a `Gridspace{Stack}` when a direct argument varies.
"""
function rope(
        item;
        layers,
        n = 6,
        lay = nothing,
        dir = 1,
        φ0 = 0,
        compact = nothing,
        gap_frac = 0,
        combine::Symbol = :product
)
    caller = function (
            resolved_item, resolved_layers, resolved_n, resolved_lay,
            resolved_dir, resolved_φ0, resolved_compact, resolved_gap
    )
        resolved_layers isa Integer && !(resolved_layers isa Bool) &&
        resolved_layers >= 0 || throw(ArgumentError(
            "layers must be a nonnegative integer"
        ))
        count = Int(resolved_layers)
        return _rope(
            resolved_item,
            count,
            _count_schedule(resolved_n, count),
            _course_schedule(resolved_lay, count, :lay),
            _course_schedule(resolved_dir, count, :dir),
            _course_schedule(resolved_φ0, count, :φ0),
            _course_schedule(resolved_compact, count, :compact),
            _course_schedule(resolved_gap, count, :gap_frac)
        )
    end
    values = (item, layers, n, lay, dir, φ0, compact, gap_frac)
    return parameterize(DataModel.Stack, caller, values; combine)
end

"""
$(TYPEDSIGNATURES)

Declare a repeated armor-wire course whose radius is resolved from the current
outer boundary.

# Arguments

- `material`: material with `kind == :conductor`.

# Keywords

- `shape`: intrinsic armor-member primitive.
- `n`: exact cardinality or deferred maximum count `capacity()`.
- `lay=nothing`: one helical lay law.
- `dir=1`: helix handedness, `1` or `-1` \\[dimensionless\\].
- `φ0=0`: initial angular position \\[rad\\].
- `compact=nothing`: explicit compaction law.
- `gap_frac=0`: fractional adjacent clearance \\[dimensionless\\].
- `tag=:armor`: local member and group identity.
- `combine=:product`: gridspace composition rule.

# Returns

- An armor `Group`, or a `Gridspace{Group}` when a direct argument varies.
"""
function armor(
        material;
        shape,
        n,
        lay = nothing,
        dir = 1,
        φ0 = 0,
        compact = nothing,
        gap_frac = 0,
        tag = :armor,
        combine::Symbol = :product
)
    caller = function (
            resolved_material, resolved_shape, resolved_n, resolved_lay,
            resolved_dir, resolved_φ0, resolved_compact, resolved_gap,
            resolved_tag
    )
        source = DataModel.Region(
            resolved_tag,
            resolved_shape,
            validate(resolved_material, DataModel.Region, :armor, (:conductor,))
        )
        return DataModel.Group(
            resolved_tag,
            DataModel.Pose2(0, 0, 0),
            source,
            DataModel.Ring(
                resolved_n;
                r = nothing,
                φ0 = resolved_φ0,
                gap_frac = resolved_gap
            ),
            _path(resolved_lay, resolved_dir, resolved_φ0),
            resolved_compact
        )
    end
    values = (material, shape, n, lay, dir, φ0, compact, gap_frac, tag)
    return parameterize(DataModel.Group, caller, values; combine)
end

function _tape(material, section, n, lay, gap, compact, tag)
    section isa DataModel.Rectangle || throw(ArgumentError(
        "tape section must be a Rectangle(width, thickness)"
    ))
    source = DataModel.Region(tag, section, material)
    return DataModel.Group(
        :tapes,
        DataModel.Pose2(0, 0, 0),
        source,
        DataModel.Ring(n; r = nothing, gap_frac = gap),
        _path(lay, 1, 0),
        compact
    )
end

"""
$(TYPEDSIGNATURES)

Declare a repeated conductive, semiconductive, or insulating tape system.

# Arguments

- `material`: tape material.

# Keywords

- `section`: intrinsic tape cross-section.
- `n`: exact angular cardinality or deferred maximum count `capacity()`.
- `lay=nothing`: one helical lay law.
- `gap_frac=0`: fractional angular clearance \\[dimensionless\\].
- `compact=nothing`: explicit compaction law.
- `tag=:tape`: local tape identity.
- `combine=:product`: gridspace composition rule.

# Returns

- A tape `Group`, or a `Gridspace{Group}` when a direct argument varies.
"""
function tape(
        material;
        section,
        n,
        lay = nothing,
        gap_frac = 0,
        compact = nothing,
        tag = :tape,
        combine::Symbol = :product
)
    values = (material, section, n, lay, gap_frac, compact, tag)
    return parameterize(DataModel.Group, _tape, values; combine)
end

"""
$(TYPEDSIGNATURES)

Arrange independent cable parts as repeated or explicit cores.

The pattern-backed form retains one prototype. The variadic form preserves
heterogeneous members and their local poses.

# Arguments

- `members`: explicit core members.

# Keywords

- `n=nothing`: repeated cardinality. Omit for explicit members.
- `r=nothing`: member-center ring radius \\[m\\], required when `n` is supplied.
- `names=nothing`: exact terminal names required for repeated terminal-bearing members.
- `φ0=0`: starting angle \\[rad\\].
- `span=2π`: angular span \\[rad\\].
- `path=nothing`: shared longitudinal path.
- `compact=nothing`: explicit compaction law.
- `combine=:product`: gridspace composition rule.

# Returns

- An `Assembly`, or a `Gridspace{Assembly}` when a direct argument varies.

# Notes

Origin-centered repeated sectors require the sector span to equal their angular
pitch. Their resolved sides must not overlap. An outer insulating layer may
close the clearance to zero. Bare sectors require positive side clearance.
"""
function cores(
        members...;
        n = nothing,
        r = nothing,
        names = nothing,
        φ0 = 0,
        span = 2π,
        path = nothing,
        compact = nothing,
        combine::Symbol = :product
)
    if n === nothing
        r === nothing && φ0 == 0 && span == 2π || throw(ArgumentError(
            "cores requires n when specifying a repeated radius or angular placement"))
        return assembly(members...; names, path, compact, combine)
    end
    r === nothing && throw(ArgumentError("repeated cores require an explicit ring radius r"))
    return assembly(members...;
        pattern = DataModel.Ring(n; r, φ0, span, combine),
        names, path, compact, combine)
end

function _enclosure_item(items::Tuple, formation)
    if formation === nothing
        return length(items) == 1 && only(items) isa DataModel.AbstractCablePart ?
               only(items) : _explicit_assembly(items...)
    end
    length(items) == 1 || throw(ArgumentError(
        "a pattern-backed enclosure requires one repeated prototype"
    ))
    only(items) isa DataModel.AbstractCablePart || throw(ArgumentError(
        "a pattern-backed enclosure requires an unplaced physical prototype"
    ))
    return DataModel.Assembly(
        DataModel.Pose2(0, 0, 0),
        only(items),
        formation,
        nothing,
        nothing,
        nothing
    )
end

function _enclose(tag, items, shape, fill, wall, pose, formation)
    item = _enclosure_item(Tuple(items), formation)
    resolved_pose = pose === nothing ? DataModel.Pose2(0, 0, 0) : pose
    return DataModel.Enclosure(tag, resolved_pose, shape, item, fill, wall)
end

"""
$(TYPEDSIGNATURES)

Contain one or more physical members inside a pipe cross-section.

# Arguments

- `items`: enclosed physical members. Several members form an explicit
  assembly.

# Keywords

- `shape`: intrinsic containing primitive.
- `fill`: filling material or explicit filling region.
- `wall=nothing`: optional outward wall declaration.
- `at=nothing`: pipe pose relative to its parent frame.
- `combine=:product`: gridspace composition rule.

# Returns

- An `Enclosure`, or a `Gridspace{Enclosure}` when a direct argument varies.
"""
function pipe(
        items...;
        shape,
        fill,
        wall = nothing,
        at = nothing,
        combine::Symbol = :product
)
    isempty(items) && throw(ArgumentError("pipe requires enclosed content"))
    caller = (selected...) -> begin
        count = length(items)
        physical = selected[1:count]
        _enclose(:pipe, physical, selected[(count + 1):end]..., nothing)
    end
    values = (items..., shape, fill, wall, at)
    return parameterize(DataModel.Enclosure, caller, values; combine)
end

"""
$(TYPEDSIGNATURES)

Contain one or more physical members inside a duct cross-section.

A `formation` repeats one prototype without expanding it. Explicitly placed
members may differ in geometry and terminal structure.

# Arguments

- `items`: enclosed physical members.

# Keywords

- `shape`: intrinsic containing primitive.
- `fill`: filling material or explicit filling region.
- `wall=nothing`: optional outward wall declaration.
- `formation=nothing`: placement pattern for one repeated prototype.
- `at=nothing`: duct pose relative to its parent frame.
- `combine=:product`: gridspace composition rule.

# Returns

- An `Enclosure`, or a `Gridspace{Enclosure}` when a direct argument varies.
"""
function duct(
        items...;
        shape,
        fill,
        wall = nothing,
        formation = nothing,
        at = nothing,
        combine::Symbol = :product
)
    isempty(items) && throw(ArgumentError("duct requires enclosed content"))
    caller = (selected...) -> begin
        count = length(items)
        physical = selected[1:count]
        _enclose(:duct, physical, selected[(count + 1):end]...)
    end
    values = (items..., shape, fill, wall, at, formation)
    return parameterize(DataModel.Enclosure, caller, values; combine)
end

function build(
        ::Type{DataModel.CableDesign},
        cable_id,
        parts::Tuple,
        nominal_data;
        combine::Symbol = :product
)
    isempty(parts) && throw(ArgumentError("a cable design requires one physical part"))
    caller = (selected...) -> begin
        id = first(selected)
        data = last(selected)
        physical = selected[2:(end - 1)]
        build(DataModel.CableDesign, id, Tuple(physical), data)
    end
    values = (cable_id, parts..., nominal_data)
    return parameterize(DataModel.CableDesign, caller, values; combine)
end

function build(
        ::Type{DataModel.CableDesign},
        cable_id,
        parts...;
        nominal_data = nothing,
        combine::Symbol = :product
)
    values = (cable_id, parts..., nominal_data)
    any(value -> value isa Union{AbstractGrid, Gridspace}, values) || throw(MethodError(
        build, (DataModel.CableDesign, cable_id, parts...)
    ))
    return build(
        DataModel.CableDesign,
        cable_id,
        Tuple(parts),
        nominal_data;
        combine
    )
end
