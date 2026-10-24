struct LayersBuilder <: Builder{Stack} end

"""
```julia
layers(parts...; combine = :product)
```

Compose physical parts in outward order.

# Arguments

- `parts`: physical declarations ordered from the center outward.

# Keywords

- `combine=:product`: composition of the sources when an argument or a keyword is a `Grid`
  or a `Gridspace`, `:product` or `:zip`.

# Returns

- A `Stack`, or a `Gridspace{Stack}` when an argument varies.
"""
const layers = LayersBuilder()

build(::LayersBuilder, parts...) = Stack(parts...)

function _explicit_assembly(members...)
    placed = map(_assembly_member, members)
    return DataModel.Assembly(
        DataModel.Pose2(0, 0, 0), Tuple(placed), nothing, nothing, nothing, nothing
    )
end

# An explicit assembly takes placed members, and it places a bare cable part at its origin.
_assembly_member(member::DataModel.AssemblyMember) = member
_assembly_member(member::DataModel.AbstractCablePart) = DataModel.AssemblyMember(member)
function _assembly_member(member)
    throw(ArgumentError(
        "assembly members must be physical cable parts; got $(nameof(typeof(member)))"
    ))
end

struct AssemblyBuilder <: Builder{Assembly} end

"""
```julia
assembly(members...; pattern = nothing, names = nothing, path = nothing,
    compact = nothing, combine = :product)
```

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
- `combine=:product`: composition of the sources when an argument or a keyword is a `Grid`
  or a `Gridspace`, `:product` or `:zip`.

# Returns

- An `Assembly`, or a `Gridspace{Assembly}` when an argument or a keyword varies.
"""
const assembly = AssemblyBuilder()

function build(::AssemblyBuilder, members...; pattern = nothing, names = nothing,
        path = nothing, compact = nothing)
    isempty(members) && throw(ArgumentError("assembly requires at least one member"))
    if pattern === nothing
        names === nothing && path === nothing && compact === nothing ||
            throw(ArgumentError(
                "explicit assembly members own their names, paths and compaction; " *
                "provide pattern to repeat one prototype"))
        return _explicit_assembly(members...)
    end
    length(members) == 1 || throw(ArgumentError(
        "a repeated assembly requires exactly one prototype with pattern"))
    return Assembly(first(members); pattern, names, path, compact)
end

_terminal_eligible(::DataModel.Region{<:Any, <:Materials.Material{<:Real, :conductor}}) = 1
_terminal_eligible(::DataModel.Region) = 0
_terminal_eligible(stack::DataModel.Stack) = sum(_terminal_eligible, stack.items; init = 0)
_terminal_eligible(group::DataModel.Group) = _terminal_eligible(group.member)
function _terminal_eligible(assembly::DataModel.Assembly)
    return sum(member -> _terminal_eligible(member.part), assembly.members; init = 0)
end
# A material fill has no terminal.
_terminal_eligible(::Materials.Material) = 0
function _terminal_eligible(enclosure::DataModel.Enclosure)
    count = _terminal_eligible(enclosure.contents)
    count += _terminal_eligible(enclosure.fill)
    enclosure.wall === nothing || (count += _terminal_eligible(enclosure.wall))
    return count
end

struct TerminalBuilder <: Builder{Group} end

"""
```julia
terminal(name, parts...; combine = :product)
```

Coalesce every conductive descendant of an ordered physical subtree into one
retained terminal.

# Arguments

- `name`: retained electrical terminal name.
- `parts`: physical declarations ordered from the center outward.

# Keywords

- `combine=:product`: composition of the sources when an argument or a keyword is a `Grid`
  or a `Gridspace`, `:product` or `:zip`.

# Returns

- A terminal-owning `Group`, or a `Gridspace{Group}` when an argument varies.

# Errors

- Throws `ArgumentError` when the realized subtree contains no conductive
  descendant.
"""
const terminal = TerminalBuilder()

function build(::TerminalBuilder, name, parts...)
    isempty(parts) && throw(ArgumentError("terminal requires at least one part"))
    root = DataModel.Stack(parts...)
    _terminal_eligible(root) > 0 || throw(ArgumentError(
        "terminal :$name requires a conductive descendant"
    ))
    return DataModel.Group(
        name, DataModel.Pose2(0, 0, 0), root, nothing, nothing, nothing
    )
end

struct SolidBuilder <: Builder{Region} end

"""
```julia
solid(material, primitive = nothing; tag = :solid, combine = :product)
```

Bind an intrinsic primitive definition to one material and local physical tag.

# Arguments

- `material`: constitutive material record.
- `primitive=nothing`: intrinsic cross-sectional primitive definition. Omitted, the region
  takes the outline of its context, which only the center of a `milliken` supplies.

# Keywords

- `tag=:solid`: local physical identity.
- `combine=:product`: composition of the sources when an argument or a keyword is a `Grid`
  or a `Gridspace`, `:product` or `:zip`.

# Returns

- A `Region`, or a `Gridspace{Region}` when an argument or a keyword varies.
"""
const solid = SolidBuilder()

build(::SolidBuilder, material, primitive = nothing; tag = :solid) =
    Region(tag, primitive, material)

struct ShellBuilder <: Builder{Region} end

"""
```julia
shell(material; t, tag = :shell, combine = :product)
```

Declare one outward material layer of thickness `t` \\[m\\].

# Arguments

- `material`: constitutive material record.

# Keywords

- `t`: normal layer thickness \\[m\\].
- `tag=:shell`: local physical identity.
- `combine=:product`: composition of the sources when an argument or a keyword is a `Grid`
  or a `Gridspace`, `:product` or `:zip`.

# Returns

- A contextual `Region`, or a `Gridspace{Region}` when an argument or a keyword varies.
"""
const shell = ShellBuilder()

build(::ShellBuilder, material; t, tag = :shell) = Region(tag, Shell(t), material)

struct CoreBuilder <: Builder{Region} end

"""
```julia
core(material, primitive; tag = :core, combine = :product)
core(material; r, tag = :core, combine = :product)
```

Declare a conductive core region.

# Arguments

- `material`: material with `kind == :conductor`.
- `primitive`: intrinsic core geometry. The keyword form constructs a disk of
  radius `r` \\[m\\].

# Keywords

- `r`: disk radius \\[m\\].
- `tag=:core`: local physical identity.
- `combine=:product`: composition of the sources when an argument or a keyword is a `Grid`
  or a `Gridspace`, `:product` or `:zip`.

# Returns

- A conductive `Region`, or a `Gridspace{Region}` when an argument or a keyword varies.
"""
const core = CoreBuilder()

function build(::CoreBuilder, material, primitive; tag = :core)
    return DataModel.Region(
        tag, primitive, validate(material, DataModel.Region, :core, (:conductor,)))
end
build(builder::CoreBuilder, material; r, tag = :core) =
    build(builder, material, Disk(r); tag)

function _role_shell(role::Symbol, allowed, material, t, tag)
    return DataModel.Region(
        tag,
        DataModel.Shell(t),
        validate(material, DataModel.Region, role, allowed)
    )
end

struct InsulationBuilder <: Builder{Region} end

"""
```julia
insulation(material; t, tag = :insulation, combine = :product)
```

Declare an insulating layer of thickness `t` \\[m\\] around the preceding part.

# Arguments

- `material`: material with `kind == :insulator`.

# Keywords

- `t`: normal layer thickness \\[m\\].
- `tag=:insulation`: local physical identity.
- `combine=:product`: composition of the sources when an argument or a keyword is a `Grid`
  or a `Gridspace`, `:product` or `:zip`.

# Returns

- A `Region`, or a `Gridspace{Region}` when an argument or a keyword varies.

# Examples

```julia
xlpe = Material(:insulator, 1.97e14, 2.3)
pe = Material(:insulator, 1e15, 2.25)
insulation(xlpe; t = 5.5e-3)                                    # one Region
insulation(xlpe; t = Grid((5.5e-3, 6e-3)))                      # two Regions
insulation(Grid((pe, xlpe)); t = Grid((5.5e-3, 6e-3)), combine = :zip)   # pe with 5.5 mm
```

# Extension

`insulation` lifts every call in one grammar method, then builds each point from scalars
through `build`. A material type of your own needs that scalar method only:

```julia
LineCableModels.build(::typeof(insulation), m::MyMaterial; t, tag = :insulation) =
    Region(tag, Shell(t), Material(:insulator, m.rho, m.eps_r))
```
"""
const insulation = InsulationBuilder()

build(::InsulationBuilder, material; t, tag = :insulation) =
    _role_shell(:insulation, (:insulator,), material, t, tag)

struct ScreenBuilder <: Builder{Region} end

"""
```julia
screen(material; t, tag = :screen, combine = :product)
```

Declare a semiconductive or conductive screen layer of thickness `t` \\[m\\]. It takes the
keywords of [`insulation`](@ref) and lifts the same way.
"""
const screen = ScreenBuilder()

build(::ScreenBuilder, material; t, tag = :screen) =
    _role_shell(:screen, (:semicon, :conductor), material, t, tag)

struct SheathBuilder <: Builder{Region} end

"""
```julia
sheath(material; t, tag = :sheath, combine = :product)
```

Declare a conductive sheath layer of thickness `t` \\[m\\]. It takes the keywords of
[`insulation`](@ref) and lifts the same way.
"""
const sheath = SheathBuilder()

build(::SheathBuilder, material; t, tag = :sheath) =
    _role_shell(:sheath, (:conductor,), material, t, tag)

struct BeddingBuilder <: Builder{Region} end

"""
```julia
bedding(material; t, tag = :bedding, combine = :product)
```

Declare a nonconducting bedding layer of thickness `t` \\[m\\]. It takes the keywords of
[`insulation`](@ref) and lifts the same way.
"""
const bedding = BeddingBuilder()

build(::BeddingBuilder, material; t, tag = :bedding) =
    _role_shell(:bedding, (:insulator,), material, t, tag)

struct JacketBuilder <: Builder{Region} end

"""
```julia
jacket(material; t, tag = :jacket, combine = :product)
```

Declare an insulating jacket layer of thickness `t` \\[m\\]. It takes the keywords of
[`insulation`](@ref) and lifts the same way.
"""
const jacket = JacketBuilder()

build(::JacketBuilder, material; t, tag = :jacket) =
    _role_shell(:jacket, (:insulator,), material, t, tag)

struct FillerBuilder <: Builder{Region} end

"""
```julia
filler(material, primitive; tag = :filler, combine = :product)
```

Bind a nonconducting filler material to an intrinsic primitive definition.

# Arguments

- `material`: material with `kind == :insulator`.
- `primitive`: intrinsic filler geometry.

# Keywords

- `tag=:filler`: local physical identity.
- `combine=:product`: composition of the sources when an argument or a keyword is a `Grid`
  or a `Gridspace`, `:product` or `:zip`.

# Returns

- A filler `Region`, or a `Gridspace{Region}` when an argument or a keyword varies.
"""
const filler = FillerBuilder()

function build(::FillerBuilder, material, primitive; tag = :filler)
    return DataModel.Region(
        tag, primitive, validate(material, DataModel.Region, :filler, (:insulator,)))
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

struct WiresBuilder <: Builder{Group} end

"""
```julia
wires(material; shape, pattern = nothing, path = nothing, n = nothing, r = nothing,
    gap_frac = 0, lay = nothing, dir = 1, φ0 = 0, compact = nothing, tag = :wire,
    combine = :product)
```

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
- `combine=:product`: composition of the sources when an argument or a keyword is a `Grid`
  or a `Gridspace`, `:product` or `:zip`.

# Returns

- A repeated-member `Group`, or a `Gridspace{Group}` when an argument or a keyword varies.
"""
const wires = WiresBuilder()

function build(::WiresBuilder, material; shape, pattern = nothing, path = nothing,
        n = nothing, r = nothing, gap_frac = 0, lay = nothing, dir = 1, φ0 = 0,
        compact = nothing, tag = :wire)
    if pattern !== nothing
        n === nothing && r === nothing || throw(ArgumentError(
            "pattern cannot be combined with n or r"
        ))
        lay === nothing || throw(ArgumentError(
            "pattern-oriented wires use path rather than lay"
        ))
        return _wires(material, shape, pattern, path, compact, tag)
    end
    n === nothing && throw(ArgumentError("ring-course wires require n"))
    r === nothing && throw(ArgumentError("ring-course wires require r"))
    path === nothing || throw(ArgumentError("ring-course wires use lay rather than path"))
    return _wires(material, shape, DataModel.Ring(n; r, φ0, gap_frac), _path(lay, dir, φ0),
        compact, tag)
end

# A per-course keyword takes one value for every course, or one value per course.
function _course_schedule(values::Union{Tuple, AbstractVector}, count::Int, name::Symbol)
    length(values) == count || throw(DimensionMismatch(
        "$name requires exactly $count course values"
    ))
    return Tuple(values)
end
_course_schedule(value, count::Int, name::Symbol) = ntuple(_ -> value, count)

# `n` declares the member counts of the courses: `capacity()` infers them, a positive base
# count puts k times it in course k, and a schedule gives one count per course from the
# inside out. A Bool is not a count. `stranded` and `rope` share this meaning.
_course_declaration(values::Union{Tuple, AbstractVector}) = values
function _course_declaration(value::Integer)
    value > 0 && return Int(value)
    throw(ArgumentError(
        "n must be a positive base count, exact course schedule, or capacity()"
    ))
end
function _course_declaration(value::Bool)
    throw(ArgumentError(
        "n must be a positive base count, exact course schedule, or capacity()"
    ))
end
function _course_declaration(value)
    value === DataModel.capacity() && return nothing
    throw(ArgumentError(
        "n must be a positive base count, exact course schedule, or capacity()"
    ))
end

# A rope expands its course declaration over its layers.
function _count_schedule(values::Union{Tuple, AbstractVector}, count::Int)
    length(values) == count || throw(DimensionMismatch(
        "n requires exactly $count course values"
    ))
    return Tuple(values)
end
_count_schedule(base::Int, count::Int) = ntuple(index -> index * base, count)
_count_schedule(::Nothing, count::Int) = ntuple(_ -> DataModel.capacity(), count)

# A stranded core stores its declaration, and the formation expands it when it resolves. A
# schedule is stored as a vector, as it is serialized.
_stranding_counts(values::Union{Tuple, AbstractVector}) = collect(values)
_stranding_counts(declaration) = declaration

# Without a lay, a stranded core has no paths. A lay schedule gives one path per inferred
# course, and one lay law gives one path for every course.
_stranding_paths(::Nothing, dir, φ0) = nothing
function _stranding_paths(lay::Union{Tuple, AbstractVector}, dir, φ0)
    count = length(lay)
    count > 0 || throw(ArgumentError("lay schedule cannot be empty"))
    directions = _course_schedule(dir, count, :dir)
    angles = _course_schedule(φ0, count, :φ0)
    return ntuple(index -> _path(lay[index], directions[index], angles[index]), count)
end
function _stranding_paths(lay, dir, φ0)
    return _path(lay, validate(dir, stranded, :dir), validate(φ0, stranded, :φ0))
end

# In a sector, a rectangle, an ellipse or a star-shaped polygon, the outline deforms the
# Disk wires around a center strand that `shape` gives.
function _stranded(
        material,
        center,
        shape::DataModel.Disk,
        counts,
        lay,
        dir,
        φ0,
        compact,
        prescribed_boundary::Union{
            DataModel.Sector, DataModel.Rectangle, DataModel.Ellipse, DataModel.Polygon},
        fill
)
    material = validate(material, DataModel.Region, :stranded, (:conductor,))
    fill = validate(fill, DataModel.Region, :stranded_fill, (:insulator, :semicon))
    center === nothing || throw(ArgumentError(
        "a non-circular stranded core infers its center strand from shape; omit center"
    ))
    compact === nothing || throw(ArgumentError(
        "a non-circular boundary deforms its wires intrinsically; omit compact"
    ))
    bounded = _stranded_core(
        material, nothing, shape, counts, _stranding_paths(lay, dir, φ0), nothing,
        prescribed_boundary)
    return DataModel.Enclosure(
        :stranded, DataModel.Pose2(0, 0, 0), prescribed_boundary, bounded, fill, nothing)
end

# In a circular core, Disk wires surround one circular center wire. The core stores a
# declared center only, and resolution places the default, a copy of `shape`, at the cable
# center. Only `compact = true` deforms the wires. Without a geometric boundary, the core
# takes the outline of its context.
function _stranded(
        material,
        center::Union{Nothing, DataModel.Disk},
        shape::DataModel.Disk,
        counts,
        lay,
        dir,
        φ0,
        compact::Union{Nothing, Bool},
        prescribed_boundary::Union{Nothing, DataModel.Disk},
        fill
)
    material = validate(material, DataModel.Region, :stranded, (:conductor,))
    fill = validate(fill, DataModel.Region, :stranded_fill, (:insulator, :semicon))
    compaction = compact === nothing ? false : compact
    bounded = _stranded_core(
        material, center, shape, counts, _stranding_paths(lay, dir, φ0), compaction,
        prescribed_boundary)
    return DataModel.Enclosure(
        :stranded, DataModel.Pose2(0, 0, 0), prescribed_boundary, bounded, fill, nothing)
end

# In a rectangular core, bent courses of a Rectangle strand surround a declared Disk center
# wire, or the outline of a preceding part. The bounded Group is the part, because its
# contiguous courses fill their whole disk.
function _stranded(
        material,
        center::Union{Nothing, DataModel.Disk},
        shape::DataModel.Rectangle,
        counts,
        lay,
        dir,
        φ0,
        compact,
        prescribed_boundary::DataModel.Disk,
        fill
)
    material = validate(material, DataModel.Region, :stranded, (:conductor,))
    fill = validate(fill, DataModel.Region, :stranded_fill, (:insulator, :semicon))
    compact === nothing || throw(ArgumentError(
        "rectangular strands are intrinsically bent; omit compact"
    ))
    return _stranded_core(
        material, center, shape, counts, _stranding_paths(lay, dir, φ0), true,
        prescribed_boundary)
end

function _stranded(
        material,
        center::Union{Nothing, DataModel.Disk},
        shape::DataModel.Disk,
        counts,
        lay,
        dir,
        φ0,
        compact,
        prescribed_boundary::Union{Nothing, DataModel.Disk},
        fill
)
    throw(ArgumentError(
        "circular stranded compaction is selected with compact=true or false"
    ))
end

function _stranded(
        material,
        center,
        shape::DataModel.Disk,
        counts,
        lay,
        dir,
        φ0,
        compact,
        prescribed_boundary::Union{Nothing, DataModel.Disk},
        fill
)
    throw(ArgumentError(
        "a disk-bounded stranded core requires one circular center wire"
    ))
end

function _stranded(
        material,
        center,
        shape::DataModel.Rectangle,
        counts,
        lay,
        dir,
        φ0,
        compact,
        prescribed_boundary::DataModel.Disk,
        fill
)
    throw(ArgumentError("a rectangular stranded core requires center=Disk(...)"))
end

function _stranded(
        material,
        center,
        shape::DataModel.Rectangle,
        counts,
        lay,
        dir,
        φ0,
        compact,
        prescribed_boundary::Union{
            DataModel.Sector, DataModel.Rectangle, DataModel.Ellipse, DataModel.Polygon},
        fill
)
    throw(ArgumentError(
        "rectangular strands take a Disk boundary; a Sector, Rectangle, Ellipse or " *
        "Polygon boundary takes Disk wires"
    ))
end

function _stranded(
        material,
        center,
        shape::Union{DataModel.Disk, DataModel.Rectangle},
        counts,
        lay,
        dir,
        φ0,
        compact,
        prescribed_boundary,
        fill
)
    throw(ArgumentError(
        "stranded requires a nonhollow Disk, Sector, Rectangle, Ellipse or Polygon boundary"
    ))
end

function _stranded(material, center, shape, counts, lay, dir, φ0, compact,
        prescribed_boundary, fill)
    throw(ArgumentError("stranded members must be Disk wires or Rectangle strands"))
end

# The bounded Group of a stranded core contains one Stack: the optional center strand,
# then the course wires with their declared counts and their paths.
function _stranded_core(material, center, shape, counts, paths, compaction,
        prescribed_boundary)
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
            counts,
            paths,
            nothing
        ))
    return DataModel.Group(
        :core,
        DataModel.Pose2(0, 0, 0),
        DataModel.Stack(parts),
        nothing,
        nothing,
        compaction,
        prescribed_boundary
    )
end

struct StrandedBuilder <: Builder{DataModel.AbstractCablePart} end

"""
```julia
stranded(material; center = nothing, shape, n = capacity(), lay = nothing, dir = 1,
    φ0 = 0, compact = nothing, boundary = nothing, fill = air, combine = :product)
```

Fill one authoritative core boundary with equal source strands, in inferred or declared
courses.

Circular source bundles have one center strand. Unless `n` declares the counts, natural
wires fill the first course with as many wires as fit around it, and each next course with
six more. Compacted wires take `6k` in course `k`. Declared natural courses must fit on
their radii inside the disk boundary. Declared compacted courses take their radii from
their areas and must fit the area of the disk boundary. Declared rectangular courses must
fill the disk boundary within one strand. The core stores the declaration only, never the inferred counts.

After a part, the wires fill the ring between the circular outline of that part and the
disk boundary, and the formation has no center of its own. Natural wires and bent
rectangular courses start at that outline, by the same rules.

A `Sector`, `Rectangle`, `Ellipse` or `Polygon` core boundary holds a center wire and
concentric courses, mapped into the core boundary by the area swept from its centroid.
Unless `n` declares the counts, it admits the largest complete inventory that meets
``[1 + 3L(L+1)]\\,\\pi a^2 \\leq A``, with `6k` wires in course `k`. The declared
counts and the center wire must fit the area of the core boundary. A `Polygon` must be
star-shaped about its centroid. The mapped sites define a prescribed-area power diagram,
and disks clipped to those cells keep each source area. This geometric reconstruction
preserves the declared copper fill and does not model mechanical forming or plasticity.

Rectangular strands occupy complete, contiguous annular courses, including a full annulus
when a course contains one strand. The disk boundary is their finished outline, as for
circular strands, and the course heights scale so that the outermost course ends on it.
The resistance is that of the declared rectangular strands. Subsequent layers start at the
disk boundary, without an outer filler film.

# Arguments

- `material`: material with `kind == :conductor`.

# Keywords

- `center=nothing`: circular center member for a disk-bounded core, stored only when
  declared. At the cable center, circular strands default to one center wire equal to
  `shape`, and rectangular strands require an explicit `Disk`. A formation after a part
  does not admit a center. A non-circular boundary infers its center strand from `shape`
  and does not admit a separate center declaration.
- `shape`: circular wire or rectangular strand primitive.
- `n=capacity()`: the course counts. `capacity()` infers them, a positive base count puts
  `k n` members in course `k`, and a schedule gives one count per course from the inside
  out.
- `boundary`: nonhollow `Disk`, `Sector`, `Rectangle`, `Ellipse` or star-shaped `Polygon`
  core boundary, or a `Disk` around a preceding part. It is the finished outline of the
  formation. Omitted, a core of Disk wires takes the outline of its context, which only
  the center of a `milliken` supplies.
- `fill=air`: interstitial insulating or semiconducting material, unused by
  contiguous rectangular courses. The default
  is lossless air (``\\rho=\\infty``, ``\\epsilon_r=\\mu_r=1``).
- `lay=nothing`: one lay law or a schedule matching the inferred radial courses.
- `dir=1`: helix handedness or a schedule matching `lay`.
- `φ0=0`: helix initial angle or a schedule matching `lay` \\[rad\\].
- `compact=nothing`: circular disk-bounded strands remain circular by default.
  `true` requests area-preserving deformation at the cable center. Non-circular
  boundaries and rectangular strands deform intrinsically and do not admit this keyword.
- `combine=:product`: composition of the sources when an argument or a keyword is a `Grid`
  or a `Gridspace`, `:product` or `:zip`.

# Returns

- A bounded `Group` for rectangular strands, or a material-complete
  `Enclosure` for circular strands. When an argument or a keyword varies, a
  `Gridspace` of the part that the `shape` keyword selects, or of `AbstractCablePart` when
  its source mixes primitive types.
"""
const stranded = StrandedBuilder()

# The `shape` keyword selects what `stranded` builds: a bounded group for rectangular
# strands, and an enclosure for Disk wires. A source of `shape` values selects by the type of
# its points, and a mixed or empty one keeps `AbstractCablePart`. The material is typed so
# that a `shape` value passed on as the subject never returns to this method.
function target(builder::StrandedBuilder, material::ScalarOrGrid{Materials.AbstractMaterial};
        shape, keywords...)
    return target(builder, shape)
end
target(::StrandedBuilder, ::Type{<:DataModel.Rectangle}) = DataModel.Group
target(::StrandedBuilder, ::Type{<:DataModel.Disk}) = DataModel.Enclosure

function build(::StrandedBuilder, material; center = nothing, shape,
        n = DataModel.capacity(), lay = nothing, dir = 1, φ0 = 0, compact = nothing,
        boundary = nothing, fill = Materials.Material(:insulator, Inf))
    return _stranded(material, center, shape, _stranding_counts(_course_declaration(n)),
        lay, dir, φ0, compact, boundary, fill)
end

# One lay law takes one direction and one initial angle. A schedule of either requires a
# lay schedule of the same length.
validate(value, ::typeof(stranded), role::Symbol) = value
function validate(value::Union{Tuple, AbstractVector}, ::typeof(stranded), role::Symbol)
    throw(ArgumentError(
        "a $role schedule requires a lay schedule of the same length; " *
        "got $(length(value)) values"
    ))
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

# The center of a Milliken fills the disk that its segments leave. A solid or a stranded
# core declared without its outline takes that disk. As in a segment, the fill of the
# Milliken takes the interstices of a stranded core at the center.
function _milliken_center(::Nothing, material, disk)
    return _milliken_center(build(solid, material), material, disk)
end
function _milliken_center(center::DataModel.Region{Nothing}, material, disk)
    return DataModel.Region(center.tag, disk, center.material)
end
function _milliken_center(
        center::DataModel.Enclosure{<:DataModel.Pose2, Nothing, <:DataModel.Group},
        material,
        disk
)
    core = center.contents
    return DataModel.Group(
        core.name, core.at, core.member, core.pattern, core.path, core.compact, disk)
end
function _milliken_center(center, material, disk)
    throw(ArgumentError(
        "a milliken center is declared without its outline, as solid(material) or " *
        "stranded(material; shape, …)"
    ))
end

function _milliken(
        material,
        shape::DataModel.Disk,
        segment::DataModel.Sector,
        segments::Integer,
        lay,
        dir,
        φ0,
        fill,
        center
)
    material = validate(material, DataModel.Region, :milliken, (:conductor,))
    segments >= 2 ||
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
        nothing,
        lay,
        dir,
        φ0,
        nothing,
        segment,
        fill
    )
    _, _, segment_primitives = DataModel.bounded_members(
        DataModel.EmptyBoundary(), segment_part.contents)
    center_radius = minimum(_inner_radius, segment_primitives)
    center_radius > zero(center_radius) || throw(DomainError(
        center_radius,
        "milliken segment strands leave no positive radius for the center wire"
    ))
    members = DataModel.AssemblyMember[]
    center_terminal = DataModel.Group(
        :milliken_center,
        DataModel.Pose2(0, 0, 0),
        _milliken_center(center, material, DataModel.Disk(center_radius)),
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
            segment_part.contents,
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

function _milliken(
        material,
        shape::DataModel.Disk,
        segment::DataModel.Sector,
        segments,
        lay,
        dir,
        φ0,
        fill,
        center
)
    throw(ArgumentError("milliken segments must be an integer of at least two"))
end

function _milliken(material, shape::DataModel.Disk, segment, segments, lay, dir, φ0, fill,
        center)
    throw(ArgumentError("milliken segment must be one Sector boundary"))
end

function _milliken(material, shape, segment, segments, lay, dir, φ0, fill, center)
    throw(ArgumentError("milliken segment strands must be circular Disk wires"))
end

struct MillikenBuilder <: Builder{DataModel.Enclosure} end

"""
```julia
milliken(material; shape, segment, segments = 6, lay = nothing, dir = 1, φ0 = 0,
    fill = air, center = nothing, combine = :product)
```

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
- `center=nothing`: `solid(material)` by default. A part declared without its outline
  fills the center disk, such as a `stranded` core without a geometric boundary.
- `combine=:product`: composition of the sources when an argument or a keyword is a `Grid`
  or a `Gridspace`, `:product` or `:zip`.

The center-wire radius is inferred from the resolved segment packing so that
the center wire is tangent to the innermost strand of every equal segment.

# Returns

- A material-complete `Enclosure`, or a `Gridspace{Enclosure}` when an argument or a
  keyword varies.
"""
const milliken = MillikenBuilder()

build(::MillikenBuilder, material; shape, segment, segments = 6, lay = nothing,
    dir = 1, φ0 = 0, fill = Materials.Material(:insulator, Inf), center = nothing) =
    _milliken(material, shape, segment, segments, lay, dir, φ0, fill, center)

function _rope(item, course_count, counts, lays, directions, angles, compactions, gaps)
    throw(ArgumentError("rope item must be a physical cable part"))
end

function _rope(
        item::DataModel.AbstractCablePart,
        course_count,
        counts,
        lays,
        directions,
        angles,
        compactions,
        gaps
)
    throw(ArgumentError(
        "rope repeats a stranded member; declare plain wires as " *
        "stranded(material; shape, n = …, boundary)"
    ))
end

# As in a stranded core, every second course of a rope starts half its spacing later, and
# `φ0` is the initial angle of its helix. The count of a `capacity()` course is known only
# when the geometry resolves, so that course starts at zero.
_course_phase(course::Int, count::Int) = isodd(course) ? 0.0 : π / count
_course_phase(course::Int, count) = 0.0

# A rope repeats one stranded core, the output of `stranded`, as its central member and
# its outer courses.
function _rope(
        item::Union{
            DataModel.Group{<:DataModel.Pose2, <:DataModel.Stack, <:Any, <:Any, <:Any,
                <:Union{DataModel.Disk, DataModel.Sector}},
            DataModel.Enclosure{<:DataModel.Pose2, <:Union{DataModel.Disk, DataModel.Sector},
                <:DataModel.Group{<:DataModel.Pose2, <:DataModel.Stack, <:Any, <:Any,
                    <:Any, <:Union{DataModel.Disk, DataModel.Sector}}}
        },
        course_count,
        counts,
        lays,
        directions,
        angles,
        compactions,
        gaps
)
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
                    φ0 = _course_phase(course, counts[course]),
                    gap_frac = gaps[course]
                ),
                _path(lays[course], directions[course], angles[course]),
                compactions[course]
            ))
    end
    return DataModel.Stack(parts)
end

struct RopeBuilder <: Builder{Stack} end

"""
```julia
rope(item; layers, n = 6, lay = nothing, dir = 1, φ0 = 0, compact = nothing,
    gap_frac = 0, combine = :product)
```

Repeat one stranded core as a central member and concentric outer courses.

# Arguments

- `item`: stranded core repeated by the rope, the result of [`stranded`](@ref). Plain
  wires are a `stranded` formation with `n`.

# Keywords

- `layers`: number of outer courses \\[dimensionless\\].
- `n=6`: the course counts, in the forms of `stranded`. `capacity()` takes the maximum
  count of each course.
- `lay=nothing`: one lay law or one law per outer course. Homogeneous schedules
  may be declared as `LayRatio(q...)`, `Pitch(p...)` or `LayAngle(α...)`.
- `dir=1`: one handedness or one value per outer course.
- `φ0=0`: one helix initial angle or one value per outer course \\[rad\\]. As in
  `stranded`, every second course starts half its member spacing later. A course of
  `capacity()` members starts unshifted.
- `compact=nothing`: one compaction law or one law per outer course.
  Homogeneous scalar schedules may be declared as `FillFactor(η...)`.
- `gap_frac=0`: one clearance fraction or one value per outer course
  \\[dimensionless\\].
- `combine=:product`: composition of the sources when an argument or a keyword is a `Grid`
  or a `Gridspace`, `:product` or `:zip`.

# Returns

- A nested `Stack`, or a `Gridspace{Stack}` when an argument or a keyword varies.
"""
const rope = RopeBuilder()

function build(::RopeBuilder, item; layers, n = 6, lay = nothing, dir = 1, φ0 = 0,
        compact = nothing, gap_frac = 0)
    count = Int(validate(layers, rope))
    return _rope(item, count, _count_schedule(_course_declaration(n), count),
        _course_schedule(lay, count, :lay), _course_schedule(dir, count, :dir),
        _course_schedule(φ0, count, :φ0), _course_schedule(compact, count, :compact),
        _course_schedule(gap_frac, count, :gap_frac))
end

# The course count of a rope: a nonnegative integer. A Bool is not a count.
function validate(layers::Integer, ::typeof(rope))
    layers >= 0 || throw(ArgumentError("layers must be a nonnegative integer; got $layers"))
    return layers
end
function validate(layers::Bool, ::typeof(rope))
    throw(ArgumentError("layers must be a nonnegative integer; got $layers"))
end
function validate(layers, ::typeof(rope))
    throw(ArgumentError(
        "layers must be a nonnegative integer; got $(nameof(typeof(layers)))"
    ))
end

struct ArmorBuilder <: Builder{Group} end

"""
```julia
armor(material; shape, n, lay = nothing, dir = 1, φ0 = 0, compact = nothing,
    gap_frac = 0, tag = :armor, combine = :product)
```

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
- `combine=:product`: composition of the sources when an argument or a keyword is a `Grid`
  or a `Gridspace`, `:product` or `:zip`.

# Returns

- An armor `Group`, or a `Gridspace{Group}` when an argument or a keyword varies.
"""
const armor = ArmorBuilder()

function build(::ArmorBuilder, material; shape, n, lay = nothing, dir = 1, φ0 = 0,
        compact = nothing, gap_frac = 0, tag = :armor)
    source = DataModel.Region(
        tag, shape, validate(material, DataModel.Region, :armor, (:conductor,)))
    return DataModel.Group(tag, DataModel.Pose2(0, 0, 0), source,
        DataModel.Ring(n; r = nothing, φ0, gap_frac), _path(lay, dir, φ0), compact)
end

function _tape(material, section, n, lay, gap, compact, tag)
    throw(ArgumentError("tape section must be a Rectangle(width, thickness)"))
end

function _tape(material, section::DataModel.Rectangle, n, lay, gap, compact, tag)
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

struct TapeBuilder <: Builder{Group} end

"""
```julia
tape(material; section, n, lay = nothing, gap_frac = 0, compact = nothing, tag = :tape,
    combine = :product)
```

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
- `combine=:product`: composition of the sources when an argument or a keyword is a `Grid`
  or a `Gridspace`, `:product` or `:zip`.

# Returns

- A tape `Group`, or a `Gridspace{Group}` when an argument or a keyword varies.
"""
const tape = TapeBuilder()

build(::TapeBuilder, material; section, n, lay = nothing, gap_frac = 0,
    compact = nothing, tag = :tape) =
    _tape(material, section, n, lay, gap_frac, compact, tag)

struct CoresBuilder <: Builder{Assembly} end

"""
```julia
cores(members...; n = nothing, r = nothing, names = nothing, φ0 = 0, span = 2π,
    path = nothing, compact = nothing, combine = :product)
```

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
- `combine=:product`: composition of the sources when an argument or a keyword is a `Grid`
  or a `Gridspace`, `:product` or `:zip`.

# Returns

- An `Assembly`, or a `Gridspace{Assembly}` when an argument or a keyword varies.

# Notes

Origin-centered repeated sectors require the sector span to equal their angular
pitch. Their resolved sides must not overlap. An outer insulating layer may
close the clearance to zero. Bare sectors require positive side clearance.
"""
const cores = CoresBuilder()

function build(::CoresBuilder, members...; n = nothing, r = nothing, names = nothing,
        φ0 = 0, span = 2π, path = nothing, compact = nothing)
    if n === nothing
        r === nothing && φ0 == 0 && span == 2π || throw(ArgumentError(
            "cores requires n when specifying a repeated radius or angular placement"))
        return build(assembly, members...; names, path, compact)
    end
    r === nothing && throw(ArgumentError("repeated cores require an explicit ring radius r"))
    return build(assembly, members...; pattern = DataModel.Ring(n; r, φ0, span), names, path,
        compact)
end

function _enclose(tag, items, shape, fill, wall, pose)
    item = _enclosed_item(items)
    resolved_pose = pose === nothing ? DataModel.Pose2(0, 0, 0) : pose
    return DataModel.Enclosure(tag, resolved_pose, shape, item, fill, wall)
end

# Several members, or one placed member, form an explicit assembly. One cable part is the
# enclosed item.
_enclosed_item(items::Tuple{DataModel.AbstractCablePart}) = only(items)
_enclosed_item(items) = _explicit_assembly(items...)

struct PipeBuilder <: Builder{DataModel.Enclosure} end

"""
```julia
pipe(items...; shape, fill, wall = nothing, at = nothing, combine = :product)
```

Contain one or more physical members inside a pipe cross-section.

# Arguments

- `items`: enclosed physical members. Several members form an explicit
  assembly.

# Keywords

- `shape`: intrinsic containing primitive.
- `fill`: filling material or explicit filling region.
- `wall=nothing`: optional outward wall declaration.
- `at=nothing`: pipe pose relative to its parent frame.
- `combine=:product`: composition of the sources when an argument or a keyword is a `Grid`
  or a `Gridspace`, `:product` or `:zip`.

# Returns

- An `Enclosure`, or a `Gridspace{Enclosure}` when an argument or a keyword varies.
"""
const pipe = PipeBuilder()

function build(::PipeBuilder, items...; shape, fill, wall = nothing,
        at = nothing)
    isempty(items) && throw(ArgumentError(
        "pipe requires enclosed content; build a filled tube as " *
        "filler(material, Disk(ri)) followed by sheath(material; t)"
    ))
    return _enclose(:pipe, items, shape, fill, wall, at)
end

struct DuctBuilder <: Builder{DataModel.Enclosure} end

"""
```julia
duct(items...; shape, fill, wall = nothing, at = nothing, combine = :product)
```

Contain one or more physical members inside a duct cross-section.

# Arguments

- `items`: enclosed physical members. Several members form an explicit
  assembly.

# Keywords

- `shape`: intrinsic containing primitive.
- `fill`: filling material or explicit filling region.
- `wall=nothing`: optional outward wall declaration.
- `at=nothing`: duct pose relative to its parent frame.
- `combine=:product`: composition of the sources when an argument or a keyword is a `Grid`
  or a `Gridspace`, `:product` or `:zip`.

# Returns

- An `Enclosure`, or a `Gridspace{Enclosure}` when an argument or a keyword varies.
"""
const duct = DuctBuilder()

function build(::DuctBuilder, items...; shape, fill, wall = nothing,
        at = nothing)
    isempty(items) && throw(ArgumentError(
        "duct requires enclosed content; build a filled tube as " *
        "filler(material, Disk(ri)) followed by sheath(material; t)"
    ))
    return _enclose(:duct, items, shape, fill, wall, at)
end
