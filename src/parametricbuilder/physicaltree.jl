function Region(tag, primitive, material; combine::Symbol = :product)
    return parameterize(
        DataModel.Region, DataModel.Region, (tag, primitive, material); combine
    )
end

# `DataModel.Region(::Symbol, primitive, material)` is deliberately permissive
# for scalar declaration types, so explicit finite inputs need narrower methods
# to reach the common construction boundary instead of that scalar constructor.
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

- `parts`: Physical declarations ordered from the centre outward.

# Keywords

- `combine=:product`: Gridspace composition rule.

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

The pattern-oriented method retains one prototype and one placement pattern.
The variadic method retains heterogeneous members and their local poses.

# Arguments

- `members`: Physical members, optionally placed with [`at`](@ref).
- `item`: Repeated physical prototype.

# Keywords

- `pattern`: Placement pattern for a repeated prototype.
- `names=nothing`: Exact terminal names for repeated terminal-bearing members.
- `path=nothing`: Shared longitudinal path declaration.
- `compact=nothing`: Explicit compaction law.
- `combine=:product`: Gridspace composition rule.

# Returns

- An `Assembly`, or a `Gridspace{Assembly}` when a direct argument varies.
"""
function assembly(members...; combine::Symbol = :product)
    isempty(members) && throw(ArgumentError("assembly requires at least one member"))
    return parameterize(DataModel.Assembly, _explicit_assembly, members; combine)
end

function assembly(
        item;
        pattern,
        names = nothing,
        path = nothing,
        compact = nothing,
        combine::Symbol = :product
)
    return Assembly(
        item;
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

- `name`: Retained electrical terminal name.
- `parts`: Physical declarations ordered from the centre outward.

# Keywords

- `combine=:product`: Gridspace composition rule.

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

function _require_material(material, role::Symbol, allowed::Tuple)
    material isa Material || throw(ArgumentError(
        "$role material must resolve to Material"
    ))
    material.kind in allowed || throw(ArgumentError(
        "$role material must have kind $(join(string.(allowed), " or ")); " *
        "got :$(material.kind)"
    ))
    return material
end

"""
$(TYPEDSIGNATURES)

Bind an intrinsic primitive definition to one material and local physical tag.

# Arguments

- `material`: Constitutive material record.
- `primitive`: Intrinsic cross-sectional primitive definition.

# Keywords

- `tag=:solid`: Local physical identity.
- `combine=:product`: Gridspace composition rule.

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

- `material`: Constitutive material record.

# Keywords

- `t`: Normal layer thickness \\[m\\].
- `tag=:shell`: Local physical identity.
- `combine=:product`: Gridspace composition rule.

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

- `material`: Material with `kind == :conductor`.
- `primitive`: Intrinsic core geometry. The keyword form constructs a disk of
  radius `r` \\[m\\].

# Keywords

- `r`: Disk radius \\[m\\].
- `tag=:core`: Local physical identity.
- `combine=:product`: Gridspace composition rule.

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
        _require_material(resolved_material, :core, (:conductor,))
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
        _require_material(material, role, allowed)
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

- `material`: Material with `kind == :insulator`.
- `primitive`: Intrinsic filler geometry.

# Keywords

- `tag=:filler`: Local physical identity.
- `combine=:product`: Gridspace composition rule.

# Returns

- A filler `Region`, or a `Gridspace{Region}` when a direct argument varies.
"""
function filler(material, primitive; tag = :filler, combine::Symbol = :product)
    caller = (resolved_material,
        resolved_primitive,
        resolved_tag) -> DataModel.Region(
        resolved_tag,
        resolved_primitive,
        _require_material(
            resolved_material, :filler, (:insulator,)
        )
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
        _require_material(material, :wires, (:conductor,))
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

- `material`: Material with `kind == :conductor`.

# Keywords

- `shape`: Intrinsic member primitive.
- `pattern=nothing`: Explicit member placement pattern.
- `path=nothing`: Explicit longitudinal path.
- `n=nothing`: Exact ring cardinality or `capacity()`.
- `r=nothing`: Member-centre ring radius \\[m\\].
- `gap_frac=0`: Fractional adjacent clearance \\[dimensionless\\].
- `lay=nothing`: One lay law used to construct a `Helix`.
- `dir=1`: Helix handedness, `1` or `-1` \\[dimensionless\\].
- `φ0=0`: Initial angular position \\[rad\\].
- `compact=nothing`: Explicit compaction law.
- `tag=:wire`: Local member and group identity.
- `combine=:product`: Gridspace composition rule.

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

function _stranded(
        material,
        center,
        shapes,
        course_count,
        counts,
        lays,
        directions,
        angles,
        compaction,
        prescribed_boundary
)
    course_count isa Integer && !(course_count isa Bool) && course_count >= 0 ||
        throw(ArgumentError("layers must be a nonnegative integer"))
    material = _require_material(material, :stranded, (:conductor,))
    prescribed_boundary isa Union{DataModel.Disk, DataModel.Sector} ||
        throw(ArgumentError(
            "stranded requires a nonhollow Disk or Sector boundary"
        ))
    all(shape -> shape isa DataModel.Disk, shapes) || throw(ArgumentError(
        "stranded members must be circular Disk primitives; use tape for flat strips"
    ))
    if prescribed_boundary isa DataModel.Sector
        center === nothing || throw(ArgumentError(
            "a sector stranded formation has no central member"
        ))
        course_count > 0 || throw(ArgumentError(
            "a sector stranded formation requires at least one strand course"
        ))
    elseif center === nothing
        isempty(shapes) && throw(ArgumentError(
            "a zero-course stranded formation requires an explicit centre member"
        ))
        center = first(shapes)
    end
    center === nothing || center isa DataModel.Disk || throw(ArgumentError(
        "a disk stranded formation requires one circular centre member"
    ))

    natural_sector = prescribed_boundary isa DataModel.Sector && compaction === nothing
    all(counts) do count
        count isa Integer && !(count isa Bool) && count > 0 ||
            natural_sector && count === DataModel.capacity()
    end || throw(ArgumentError(
        "stranded requires positive member counts; only an uncompacted sector accepts capacity()"
    ))
    resolved_angles = Any[angles...]
    angle_zero = prescribed_boundary isa DataModel.Sector ?
                 zero(float(prescribed_boundary.span)) :
                 zero(float(prescribed_boundary.r))
    for course in eachindex(resolved_angles)
        resolved_angles[course] === nothing || continue
        resolved_angles[course] = course == 1 ?
                                  angle_zero :
                                  resolved_angles[course - 1] +
                                  π / counts[course - 1]
    end
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
    for course in 1:course_count
        pattern = DataModel.Ring(
            counts[course];
            r = nothing,
            φ0 = resolved_angles[course]
        )
        push!(parts,
            DataModel.Group(
                :strand,
                DataModel.Pose2(0, 0, 0),
                DataModel.Region(:wire, shapes[course], material),
                pattern,
                _path(lays[course], directions[course], resolved_angles[course]),
                nothing
            ))
    end
    stack = DataModel.Stack(parts)
    bounded = DataModel.Group(
        :core,
        DataModel.Pose2(0, 0, 0),
        stack,
        nothing,
        nothing,
        compaction,
        prescribed_boundary
    )
    return DataModel.Stack(bounded)
end

"""
$(TYPEDSIGNATURES)

Declare one boundary-constrained conductive core made from discrete members.

# Arguments

- `material`: Material with `kind == :conductor`.

# Keywords

- `center=nothing`: Circular centre member for a disk-bounded core. It defaults
  to `shape`. Sector-bounded cores have no centre member.
- `shape`: One circular wire primitive, or one circular primitive per course.
- `boundary`: Authoritative nonhollow `Disk` or `Sector` core boundary.
- `layers=nothing`: Number of outer courses \\[dimensionless\\]. Disk boundaries
  and compacted sector boundaries require it. An uncompacted sector infers one
  uniformly packed wire inventory.
- `n=nothing`: Positive base count or exact positive course schedule. It
  defaults to six per circular course and to `capacity()` for an uncompacted
  sector.
- `lay=nothing`: One lay law or one law per outer course. Homogeneous schedules
  may be declared as `LayRatio(q...)`, `Pitch(p...)`, or `LayAngle(α...)`.
- `dir=1`: One handedness or one value per outer course.
- `φ0=nothing`: Optional initial angle or schedule. Omission deterministically
  staggers consecutive courses.
- `compact=nothing`: Preserve natural circular wires. A `FillFactor` explicitly
  requests area-preserving deformation inside the boundary.
- `combine=:product`: Gridspace composition rule.

# Returns

- A boundary-owning core `Stack`, or a `Gridspace{Stack}` when a direct
  argument varies.
"""
function stranded(
        material;
        center = nothing,
        shape,
        layers = nothing,
        n = nothing,
        lay = nothing,
        dir = 1,
        φ0 = nothing,
        compact = nothing,
        boundary,
        combine::Symbol = :product
)
    caller = function (
            resolved_material, resolved_center, resolved_shape, resolved_layers, resolved_n,
            resolved_lay, resolved_dir, resolved_φ0, resolved_compact,
            resolved_boundary
    )
        natural_sector = resolved_boundary isa DataModel.Sector &&
                         resolved_compact === nothing
        count = if natural_sector
            resolved_layers in (nothing, 1) || throw(ArgumentError(
                "an uncompacted sector infers one uniform circular-wire inventory"
            ))
            1
        else
            resolved_layers isa Integer && !(resolved_layers isa Bool) &&
            resolved_layers >= 0 || throw(ArgumentError(
                "layers must be a nonnegative integer"
            ))
            Int(resolved_layers)
        end
        counts = _count_schedule(
            resolved_n === nothing ?
            (natural_sector ? DataModel.capacity() : 6) : resolved_n,
            count
        )
        return _stranded(
            resolved_material,
            resolved_center,
            _course_schedule(resolved_shape, count, :shape),
            count,
            counts,
            _course_schedule(resolved_lay, count, :lay),
            _course_schedule(resolved_dir, count, :dir),
            _course_schedule(resolved_φ0, count, :φ0),
            resolved_compact,
            resolved_boundary
        )
    end
    values = (
        material, center, shape, layers, n, lay, dir, φ0, compact, boundary
    )
    return parameterize(DataModel.Stack, caller, values; combine)
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

- `item`: Physical cable part repeated by the rope.

# Keywords

- `layers`: Number of outer courses \\[dimensionless\\].
- `n=6`: Base count, exact course schedule, or deferred `capacity()` policy.
- `lay=nothing`: One lay law or one law per outer course. Homogeneous schedules
  may be declared as `LayRatio(q...)`, `Pitch(p...)`, or `LayAngle(α...)`.
- `dir=1`: One handedness or one value per outer course.
- `φ0=0`: One initial angle or one value per outer course \\[rad\\].
- `compact=nothing`: One compaction law or one law per outer course.
  Homogeneous scalar schedules may be declared as `FillFactor(η...)`.
- `gap_frac=0`: One clearance fraction or one value per outer course
  \\[dimensionless\\].
- `combine=:product`: Gridspace composition rule.

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

- `material`: Material with `kind == :conductor`.

# Keywords

- `shape`: Intrinsic armor-member primitive.
- `n`: Exact cardinality or deferred `capacity()` policy.
- `lay=nothing`: One helical lay law.
- `dir=1`: Helix handedness, `1` or `-1` \\[dimensionless\\].
- `φ0=0`: Initial angular position \\[rad\\].
- `compact=nothing`: Explicit compaction law.
- `gap_frac=0`: Fractional adjacent clearance \\[dimensionless\\].
- `tag=:armor`: Local member and group identity.
- `combine=:product`: Gridspace composition rule.

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
            _require_material(resolved_material, :armor, (:conductor,))
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

- `material`: Tape material.

# Keywords

- `section`: Intrinsic tape cross-section.
- `n`: Exact angular cardinality or deferred `capacity()` policy.
- `lay=nothing`: One helical lay law.
- `gap_frac=0`: Fractional angular clearance \\[dimensionless\\].
- `compact=nothing`: Explicit compaction law.
- `tag=:tape`: Local tape identity.
- `combine=:product`: Gridspace composition rule.

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

- `item`: Repeated core prototype.
- `members`: Explicit core members.

# Keywords

- `n`: Repeated cardinality.
- `r`: Member-centre ring radius \\[m\\].
- `names`: Exact terminal names for repeated members.
- `φ0=0`: Starting angle \\[rad\\].
- `span=2π`: Angular span \\[rad\\].
- `path=nothing`: Shared longitudinal path.
- `compact=nothing`: Explicit compaction law.
- `combine=:product`: Gridspace composition rule.

# Returns

- An `Assembly`, or a `Gridspace{Assembly}` when a direct argument varies.
"""
function cores(
        item;
        n,
        r,
        names,
        φ0 = 0,
        span = 2π,
        path = nothing,
        compact = nothing,
        combine::Symbol = :product
)
    caller = function (
            resolved_item, resolved_n, resolved_r, resolved_names,
            resolved_φ0, resolved_span, resolved_path, resolved_compact
    )
        return DataModel.Assembly(
            DataModel.Pose2(0, 0, 0),
            resolved_item,
            DataModel.Ring(
                resolved_n;
                r = resolved_r,
                φ0 = resolved_φ0,
                span = resolved_span
            ),
            resolved_path,
            resolved_compact,
            resolved_names
        )
    end
    values = (item, n, r, names, φ0, span, path, compact)
    return parameterize(DataModel.Assembly, caller, values; combine)
end

cores(members...; combine::Symbol = :product) = assembly(members...; combine)

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

- `items`: Enclosed physical members. Several members form an explicit
  assembly.

# Keywords

- `shape`: Intrinsic containing primitive.
- `fill`: Filling material or explicit filling region.
- `wall=nothing`: Optional outward wall declaration.
- `at=nothing`: Pipe pose relative to its parent frame.
- `combine=:product`: Gridspace composition rule.

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

- `items`: Enclosed physical members.

# Keywords

- `shape`: Intrinsic containing primitive.
- `fill`: Filling material or explicit filling region.
- `wall=nothing`: Optional outward wall declaration.
- `formation=nothing`: Placement pattern for one repeated prototype.
- `at=nothing`: Duct pose relative to its parent frame.
- `combine=:product`: Gridspace composition rule.

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
