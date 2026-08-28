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
                _region(resolved_tag, DataModel.DiskDefinition(radius), resolved_material)
        return Gridspace{DataModel.Region}(build, grids; combine)
    end
    return _region(tag, DataModel.DiskDefinition(r), material)
end

function Shell(tag::Symbol, material; t, combine::Symbol = :product)
    values = (tag, t, material)
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        build = (resolved_tag, thickness, resolved_material) ->
                _region(resolved_tag, DataModel.ShellDefinition(thickness), resolved_material)
        return Gridspace{DataModel.Region}(build, grids; combine)
    end
    return _region(tag, DataModel.ShellDefinition(t), material)
end

_strand_path(::Nothing) = nothing
_strand_path(path::DataModel.Helix) = path
_strand_path(lay::DataModel.LayRatio) = DataModel.Helix(lay)
_strand_path(lay::Real) = DataModel.Helix(DataModel.LayRatio(lay))

function _stranded(
        tag::Symbol,
        strand::DataModel.DiskDefinition,
        material::Materials.Material,
        counts::Tuple,
        lays::Tuple,
        compact
)
    material.kind === :conductor || throw(ArgumentError(
        "conductor material must have kind :conductor; got :$(material.kind)"
    ))
    parts = DataModel.AbstractCablePart[]
    for (layer, count) in enumerate(counts)
        centre_radius = layer == 1 ? zero(strand.r) : 2 * (layer - 1) * strand.r
        layer_path = _strand_path(lays[layer])
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
            layer == 1 ? nothing : layer_path,
            compact
        ))
    end
    return DataModel.Stack(parts)
end

"""
    Conductor.Wires(tag, strand, material; n, r, lay=LayRatio(11), compact=nothing)

Construct one circular wire layer as an ordinary `Group`. `r` is the radius of
the strand-centre locus \\[m\\]. The returned group retains `tag` as its terminal
identity and preserves every strand as a resolved disk.
"""
function Wires(
        tag,
        strand,
        material;
        n,
        r,
        lay = DataModel.LayRatio(11),
        compact = nothing,
        combine::Symbol = :product
)
    values = (tag, strand, material, n, r, lay, compact)
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        caller = (resolved_tag, resolved_strand, resolved_material,
                  resolved_n, resolved_r, resolved_lay, resolved_compact) -> Wires(
            resolved_tag,
            resolved_strand,
            resolved_material;
            n = resolved_n,
            r = resolved_r,
            lay = resolved_lay,
            compact = resolved_compact
        )
        return Gridspace{DataModel.Group}(caller, grids; combine)
    end
    tag isa Symbol && strand isa DataModel.DiskDefinition && material isa Materials.Material ||
        throw(MethodError(Wires, (tag, strand, material)))
    n isa Integer && n > 0 || throw(DomainError(n, "wire count must be positive"))
    isfinite(r) && r >= zero(r) || throw(DomainError(
        r,
        "wire-layer radius must be nonnegative and finite"
    ))
    source = _region(Symbol(tag, :_wires), strand, material)
    return DataModel.Group(
        tag,
        DataModel.Pose2(0, 0, 0),
        source,
        DataModel.Ring(n; r),
        _strand_path(lay),
        compact
    )
end

"""
    Conductor.Strip(tag, sector, material; lay=LayRatio(11), compact=nothing)

Construct one helically laid conductive sector as an ordinary `Group`.
`sector` supplies its intrinsic inner and outer radii \\[m\\] and angular span
\\[rad\\].
"""
function Strip(
        tag,
        sector,
        material;
        lay = DataModel.LayRatio(11),
        compact = nothing,
        combine::Symbol = :product
)
    values = (tag, sector, material, lay, compact)
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        caller = (resolved_tag, resolved_sector, resolved_material,
                  resolved_lay, resolved_compact) -> Strip(
            resolved_tag,
            resolved_sector,
            resolved_material;
            lay = resolved_lay,
            compact = resolved_compact
        )
        return Gridspace{DataModel.Group}(caller, grids; combine)
    end
    tag isa Symbol && sector isa DataModel.SectorDefinition && material isa Materials.Material ||
        throw(MethodError(Strip, (tag, sector, material)))
    source = _region(Symbol(tag, :_strip), sector, material)
    return DataModel.Group(
        tag,
        DataModel.Pose2(0, 0, 0),
        source,
        nothing,
        _strand_path(lay),
        compact
    )
end

"""
    Conductor.Tubular(tag, annulus, material; path=nothing)

Construct one annular conductor as an ordinary `Group`. `annulus` supplies its
intrinsic inner and outer radii \\[m\\].
"""
function Tubular(
        tag,
        annulus,
        material;
        path = nothing,
        combine::Symbol = :product
)
    values = (tag, annulus, material, path)
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        caller = (resolved_tag, resolved_annulus, resolved_material,
                  resolved_path) -> Tubular(
            resolved_tag,
            resolved_annulus,
            resolved_material;
            path = resolved_path
        )
        return Gridspace{DataModel.Group}(caller, grids; combine)
    end
    tag isa Symbol && annulus isa DataModel.AnnulusDefinition && material isa Materials.Material ||
        throw(MethodError(Tubular, (tag, annulus, material)))
    source = _region(Symbol(tag, :_tubular), annulus, material)
    return DataModel.Group(
        tag,
        DataModel.Pose2(0, 0, 0),
        source,
        nothing,
        path,
        nothing
    )
end

function _stranded(
        tag::Symbol,
        strand::DataModel.RectangleDefinition,
        material::Materials.Material,
        counts::Tuple,
        lays::Tuple,
        compact
)
    material.kind === :conductor || throw(ArgumentError(
        "conductor material must have kind :conductor; got :$(material.kind)"
    ))
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

    for layer in 2:length(counts)
        count = counts[layer]
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
            _strand_path(lays[layer]),
            compact
        ))
        inner_radius = centre_radius + hypot(strand.w, strand.h) / 2
    end
    return DataModel.Stack(parts)
end

"""
    Conductor.Stranded(tag, strand, material; counts, lay=LayRatio(11), compact=nothing)

Construct a stranded conductor as an ordinary `Stack` containing one `Group`
per physical layer.

`strand` may be a `DiskDefinition` or `RectangleDefinition`. `counts` states
the member count of every physical layer, beginning with the single central
strand. Rectangular members keep their declared area and dimensions, are
oriented tangentially, and advance by the area-preserving annular-radius
relation. All layers retain terminal `tag`.

# Arguments

- `tag`: Retained terminal identity.
- `strand`: Intrinsic strand definition.
- `material`: Conductive material.

# Keywords

- `counts`: Member count of each physical layer. The first value must be `1`.
- `lay`: One `LayRatio`, `Helix`, scalar lay ratio, or `nothing` for every
  noncentral layer; alternatively, one declaration per noncentral layer. A
  full tuple including the central `nothing` is also accepted.
- `compact`: Optional ordinary group compaction declaration.
- `combine`: Gridspace composition mode.

# Returns

A `Stack`, or a `Gridspace{Stack}` when a direct argument is a finite source.
"""
function Stranded(
        tag::Symbol,
        strand::Union{DataModel.DiskDefinition, DataModel.RectangleDefinition},
        material::Materials.Material;
        counts,
        lay = DataModel.LayRatio(11),
        compact = nothing,
        combine::Symbol = :product
)
    values = (tag, strand, material, counts, lay, compact)
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        caller = (resolved_tag, resolved_strand, resolved_material,
                  resolved_counts, resolved_lay, resolved_compact) -> Stranded(
            resolved_tag,
            resolved_strand,
            resolved_material;
            counts = resolved_counts,
            lay = resolved_lay,
            compact = resolved_compact
        )
        return Gridspace{DataModel.Stack}(caller, grids; combine)
    end
    counts isa Union{Tuple, AbstractVector} || throw(ArgumentError(
        "strand counts must be a tuple or vector"
    ))
    layer_counts = Tuple(counts)
    isempty(layer_counts) && throw(ArgumentError(
        "strand counts cannot be empty"
    ))
    all(count -> count isa Integer && count > 0, layer_counts) || throw(
        DomainError(counts, "every strand-layer count must be a positive integer")
    )
    first(layer_counts) == 1 || throw(ArgumentError(
        "the first strand layer must contain the single central strand"
    ))
    layer_lays = if lay isa Union{Tuple, AbstractVector}
        declared = Tuple(lay)
        if length(declared) == length(layer_counts) - 1
            (nothing, declared...)
        elseif length(declared) == length(layer_counts)
            declared
        else
            throw(DimensionMismatch(
                "lay declarations must align with the noncentral strand layers"
            ))
        end
    else
        (nothing, ntuple(_ -> lay, length(layer_counts) - 1)...)
    end
    first(layer_lays) === nothing || throw(ArgumentError(
        "the central strand cannot have a helical lay"
    ))
    return _stranded(tag, strand, material, layer_counts, layer_lays, compact)
end

function Stranded(
        tag,
        strand,
        material;
        counts,
        lay = DataModel.LayRatio(11),
        compact = nothing,
        combine::Symbol = :product
)
    values = (tag, strand, material, counts, lay, compact)
    any(value -> value isa Union{AbstractGrid, Gridspace}, values) || throw(
        MethodError(Stranded, (tag, strand, material))
    )
    grids = map(values) do value
        value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
    end
    caller = (resolved_tag, resolved_strand, resolved_material,
              resolved_counts, resolved_lay, resolved_compact) -> Stranded(
        resolved_tag,
        resolved_strand,
        resolved_material;
        counts = resolved_counts,
        lay = resolved_lay,
        compact = resolved_compact
    )
    return Gridspace{DataModel.Stack}(caller, grids; combine)
end

end
