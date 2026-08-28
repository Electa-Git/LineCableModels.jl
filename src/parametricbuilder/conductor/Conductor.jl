"""
    LineCableModels.ParametricBuilder.Conductor

Construct conductor Regions while enforcing `Material.kind == :conductor`.
"""
module Conductor

using ..ParametricBuilder: AbstractGrid, Grid, Gridspace
import ...DataModel
import ...Materials

function _region(tag, primitive, material)
    material isa Materials.Material ||
        throw(ArgumentError("conductor material must resolve to Material"))
    material.kind === :conductor || throw(ArgumentError(
        "conductor material must have kind :conductor; got :$(material.kind)"
    ))
    return DataModel.Region(tag, primitive, material)
end

function Solid(tag::Symbol, material; r, combine::Symbol = :product)
    values = (tag, r, material)
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        build = (resolved_tag, radius, resolved_material) ->
                _region(resolved_tag, DataModel.Disk(radius), resolved_material)
        return Gridspace{DataModel.Region}(build, grids; combine)
    end
    return _region(tag, DataModel.Disk(r), material)
end

function Shell(tag::Symbol, material; t, combine::Symbol = :product)
    values = (tag, t, material)
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        build = (resolved_tag, thickness, resolved_material) ->
                _region(resolved_tag, DataModel.Shell(thickness), resolved_material)
        return Gridspace{DataModel.Region}(build, grids; combine)
    end
    return _region(tag, DataModel.Shell(t), material)
end

_strand_path(::Nothing) = nothing
_strand_path(path::DataModel.Helix) = path
_strand_path(lay::DataModel.LayRatio) = DataModel.Helix(lay)
_strand_path(lay::Real) = DataModel.Helix(DataModel.LayRatio(lay))

function _stranded(
        tag::Symbol,
        strand::DataModel.Disk,
        material::Materials.Material,
        layers::Integer,
        n::Integer,
        lay,
        compact
)
    layers > 0 || throw(DomainError(layers, "layers must be positive"))
    n > 0 || throw(DomainError(n, "the strands-per-layer multiplier must be positive"))
    material.kind === :conductor || throw(ArgumentError(
        "conductor material must have kind :conductor; got :$(material.kind)"
    ))
    path = _strand_path(lay)
    parts = DataModel.AbstractCablePart[]
    for layer in 1:Int(layers)
        count = layer == 1 ? 1 : (layer - 1) * Int(n)
        centre_radius = layer == 1 ? zero(strand.r) : 2 * (layer - 1) * strand.r
        source = DataModel.Region(
            Symbol(tag, :_strands_, layer),
            strand,
            material
        )
        push!(parts, DataModel.Group(
            tag,
            DataModel.Pose2(0, 0, 0),
            source,
            DataModel.Ring(count; r = centre_radius),
            layer == 1 ? nothing : path,
            compact
        ))
    end
    return DataModel.Stack(parts)
end

function _stranded(
        tag::Symbol,
        strand::DataModel.Rectangle,
        material::Materials.Material,
        layers::Integer,
        n::Integer,
        lay,
        compact
)
    layers > 0 || throw(DomainError(layers, "layers must be positive"))
    n > 0 || throw(DomainError(n, "the strands-per-layer multiplier must be positive"))
    material.kind === :conductor || throw(ArgumentError(
        "conductor material must have kind :conductor; got :$(material.kind)"
    ))
    path = _strand_path(lay)
    parts = DataModel.AbstractCablePart[]

    central = DataModel.Region(Symbol(tag, :_strands_, 1), strand, material)
    push!(parts, DataModel.Group(
        tag,
        DataModel.Pose2(0, 0, 0),
        central,
        nothing,
        nothing,
        compact
    ))
    inner_radius = hypot(strand.w, strand.h) / 2

    for layer in 2:Int(layers)
        count = (layer - 1) * Int(n)
        outer_radius = sqrt(
            inner_radius^2 + count * strand.w * strand.h / (one(strand.w) * pi)
        )
        centre_radius = max(
            (inner_radius + outer_radius) / 2,
            inner_radius + strand.h / 2
        )
        chord = 2 * centre_radius * sin((one(centre_radius) * pi) / count)
        tolerance = sqrt(eps(typeof(float(chord)))) * max(one(chord), chord)
        strand.w <= chord + tolerance || throw(DomainError(
            count,
            "rectangular strand layer exceeds its tangential packing limit"
        ))
        poses = map(0:(count - 1)) do member
            theta = 2 * (one(centre_radius) * pi) * member / count
            DataModel.Pose2(
                centre_radius * cos(theta),
                centre_radius * sin(theta),
                theta + (one(theta) * pi) / 2
            )
        end
        source = DataModel.Region(
            Symbol(tag, :_strands_, layer),
            strand,
            material
        )
        push!(parts, DataModel.Group(
            tag,
            DataModel.Pose2(0, 0, 0),
            source,
            poses,
            path,
            compact
        ))
        inner_radius = centre_radius + hypot(strand.w, strand.h) / 2
    end
    return DataModel.Stack(parts)
end

"""
    Conductor.Stranded(tag, strand, material; layers, n=6, lay=LayRatio(11), compact=nothing)

Construct a stranded conductor as an ordinary `Stack` containing one `Group`
per physical layer.

`strand` may be a `Disk` or `Rectangle`. Circular layers use one central strand
followed by `(layer - 1) * n` members on circular loci. Rectangular members keep
their declared area and dimensions, are oriented tangentially, and advance by
the area-preserving annular-radius relation. All layers retain terminal `tag`.

# Arguments

- `tag`: Retained terminal identity.
- `strand`: Intrinsic strand primitive.
- `material`: Conductive material.

# Keywords

- `layers`: Number of physical strand layers, including the centre.
- `n`: Member-count multiplier for noncentral layers.
- `lay`: `LayRatio`, `Helix`, scalar lay ratio, or `nothing`.
- `compact`: Optional ordinary group compaction declaration.
- `combine`: Gridspace composition mode.

# Returns

A `Stack`, or a `Gridspace{Stack}` when a direct argument is a finite source.
"""
function Stranded(
        tag::Symbol,
        strand::Union{DataModel.Disk, DataModel.Rectangle},
        material::Materials.Material;
        layers,
        n = 6,
        lay = DataModel.LayRatio(11),
        compact = nothing,
        combine::Symbol = :product
)
    values = (tag, strand, material, layers, n, lay, compact)
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        return Gridspace{DataModel.Stack}(_stranded, grids; combine)
    end
    return _stranded(tag, strand, material, layers, n, lay, compact)
end

function Stranded(
        tag,
        strand,
        material;
        layers,
        n = 6,
        lay = DataModel.LayRatio(11),
        compact = nothing,
        combine::Symbol = :product
)
    values = (tag, strand, material, layers, n, lay, compact)
    any(value -> value isa Union{AbstractGrid, Gridspace}, values) || throw(
        MethodError(Stranded, (tag, strand, material))
    )
    grids = map(values) do value
        value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
    end
    return Gridspace{DataModel.Stack}(_stranded, grids; combine)
end

end
