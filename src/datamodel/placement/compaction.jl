"""
$(TYPEDEF)

Scale nominal member-centre offsets by a measured diameter ratio.

$(TYPEDFIELDS)
"""
struct DiameterFactor{T <: Real}
    "Compacted diameter divided by nominal diameter \\[dimensionless\\]."
    k::T
    function DiameterFactor{T}(k::T) where {T <: Real}
        isfinite(k) && zero(k) < k <= one(k) ||
            throw(DomainError(k, "diameter factor must lie in (0, 1]"))
        return new{T}(k)
    end
end

_diameter_factor(k) = DiameterFactor{typeof(float(k))}(float(k))
DiameterFactor(k; combine::Symbol = :product) =
    _construction(DiameterFactor, _diameter_factor, (k,); combine)

"""
$(TYPEDEF)

Set the ratio of member material area to the resolved course-envelope area.

$(TYPEDFIELDS)
"""
struct FillFactor{T <: Real}
    "Material-area fraction \\[dimensionless\\]."
    η::T
    function FillFactor{T}(η::T) where {T <: Real}
        isfinite(η) && zero(η) < η <= one(η) || throw(DomainError(
            η, "fill factor must lie in (0, 1]"
        ))
        return new{T}(η)
    end
end

_fill_factor(η) = FillFactor{typeof(float(η))}(float(η))
FillFactor(η; combine::Symbol = :product) =
    _construction(FillFactor, _fill_factor, (η,); combine)

"""
$(TYPEDEF)

Retain manufacturer- or measurement-supplied compaction geometry.

`data` is interpreted only by primitive-specific placement methods.
"""
struct TabulatedCompaction{D}
    data::D
    TabulatedCompaction{D}(data::D) where {D} = new{D}(data)
end

_tabulated_compaction(data) = TabulatedCompaction{typeof(data)}(data)
TabulatedCompaction(data; combine::Symbol = :product) = _construction(
    TabulatedCompaction, _tabulated_compaction, (data,); combine
)

Base.:(==)(left::TabulatedCompaction, right::TabulatedCompaction) =
    left.data == right.data

"""
$(TYPEDEF)

Apply an explicit affine map during primitive-specific compaction.
"""
struct AffineCompaction{M}
    map::M
    AffineCompaction{M}(map::M) where {M} = new{M}(map)
end

_affine_compaction(map) = AffineCompaction{typeof(map)}(map)
AffineCompaction(map; combine::Symbol = :product) = _construction(
    AffineCompaction, _affine_compaction, (map,); combine
)

Base.:(==)(left::AffineCompaction, right::AffineCompaction) =
    left.map == right.map

function _fillfactor_placements(
        pattern::Ring,
        member_area::Real,
        factor::FillFactor,
        inner_radius::Real,
        count::Int
)
    outer_radius = sqrt(
        inner_radius^2 + 2count * member_area / (pattern.span * factor.η)
    )
    span = factor.η * pattern.span / count
    section = Sector(
        inner_radius,
        outer_radius,
        -span / 2,
        span
    )
    step = count == 1 ? zero(pattern.span) : pattern.span / count
    return _ResolvedPlacement[
        _ResolvedPlacement(
            Pose2(0, 0, pattern.φ0 + index * step),
            section
        )
        for index in 0:(count - 1)
    ]
end

function placements(
        pattern::Ring,
        item::Disk,
        factor::FillFactor
)
    pattern.n isa Int || throw(ArgumentError(
        "capacity() requires contextual placement"
    ))
    pattern.r === nothing && throw(ArgumentError(
        "FillFactor requires contextual placement or a resolved inner radius"
    ))
    pattern.n <= capacity(pattern, item, factor) || throw(DomainError(
        pattern.n, "compacted disk course exceeds its capacity"
    ))
    inner = max(zero(pattern.r), pattern.r - item.r)
    return _fillfactor_placements(
        pattern,
        pi * item.r^2,
        factor,
        inner,
        pattern.n
    )
end

function placements(
        pattern::Ring,
        item::Rectangle,
        factor::FillFactor
)
    pattern.n isa Int || throw(ArgumentError(
        "capacity() requires contextual placement"
    ))
    pattern.r === nothing && throw(ArgumentError(
        "FillFactor requires contextual placement or a resolved inner radius"
    ))
    pattern.n <= capacity(pattern, item, factor) || throw(DomainError(
        pattern.n, "compacted rectangle course exceeds its capacity"
    ))
    inner = max(zero(pattern.r), pattern.r - item.h / 2)
    return _fillfactor_placements(
        pattern,
        item.w * item.h,
        factor,
        inner,
        pattern.n
    )
end

function _fillfactor_capacity(
        pattern::Ring,
        member_area::Real,
        radial_half_extent::Real,
        factor::FillFactor
)
    radius = something(pattern.r, radial_half_extent)
    radius > zero(radius) || return 0
    factor.η <= inv(one(pattern.gap_frac) + pattern.gap_frac) || throw(DomainError(
        (factor.η, pattern.gap_frac),
        "fill factor cannot provide the requested member gap"
    ))
    inner = max(zero(radius), radius - radial_half_extent)
    outer = radius + radial_half_extent
    envelope_area = pattern.span * (outer^2 - inner^2) / 2
    count = factor.η * envelope_area / member_area
    return max(0, floor(Int, count + 8eps(float(count))))
end

capacity(pattern::Ring, item::Disk, factor::FillFactor) =
    _fillfactor_capacity(pattern, pi * item.r^2, item.r, factor)

capacity(pattern::Ring, item::Rectangle, factor::FillFactor) =
    _fillfactor_capacity(pattern, item.w * item.h, item.h / 2, factor)
