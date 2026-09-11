"""
$(TYPEDEF)

Set the ratio of member material area to the resolved course-envelope area.

For member areas ``A_i`` inside an envelope of area ``A_B``, the declared
dimensionless factor is

```math
\\eta = \\frac{\\sum_i A_i}{A_B}.
```

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

"""
$(TYPEDSIGNATURES)

Declare one material-area fill factor or a homogeneous schedule of fill
factors.

# Arguments

- `η`: One material-area fraction \\[dimensionless\\].
- `values`: Two or more factors, or one tuple or vector of factors, defining
  one factor per repeated course.

# Keywords

- `combine=:product`: Gridspace composition rule.

# Returns

- A `FillFactor`, a tuple of `FillFactor` values, or the corresponding
  `Gridspace` when an explicit finite source is supplied.
"""
function FillFactor(η; combine::Symbol = :product)
    parameterize(FillFactor, _fill_factor, (η,); combine)
end
function FillFactor(values::Union{Tuple, AbstractVector}; combine::Symbol = :product)
    _normalize_schedule(FillFactor, values; combine)
end
function FillFactor(first, second, remaining...; combine::Symbol = :product)
    _normalize_schedule(FillFactor, (first, second, remaining...); combine)
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
    inner = max(zero(pattern.r), pattern.r - item.h / 2)
    outer = sqrt(
        inner^2 + 2pattern.n * area(item) / (pattern.span * factor.η)
    )
    angular_width = factor.η * pattern.span / pattern.n
    definition = BentStrip(inner, outer, angular_width)
    step = pattern.n == 1 ? zero(pattern.span) : pattern.span / pattern.n
    return _ResolvedPlacement[
        _ResolvedPlacement(
            Pose2(0, 0, pattern.φ0 + index * step),
            definition
        ) for index in 0:(pattern.n - 1)
    ]
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

function capacity(pattern::Ring, item::Rectangle, factor::FillFactor)
    _fillfactor_capacity(pattern, item.w * item.h, item.h / 2, factor)
end
