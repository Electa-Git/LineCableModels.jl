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
    FillFactor{T}(η::T) where {T <: Real} = validate(new{T}(η))
end

function validate(factor::FillFactor)
    isfinite(factor.η) && zero(factor.η) < factor.η <= one(factor.η) || throw(DomainError(
        factor.η, "fill factor must lie in (0, 1]"
    ))
    return factor
end

FillFactor(η::Real) = FillFactor{typeof(float(η))}(float(η))

"""
$(TYPEDSIGNATURES)

Declare one material-area fill factor or a homogeneous schedule of fill
factors.

# Arguments

- `η`: one material-area fraction \\[dimensionless\\].
- `values`: two or more factors, or one tuple or vector of factors. Specify one factor per repeated course.

# Keywords

- `combine=:product`: gridspace composition rule.

# Returns

- A `FillFactor`, a tuple of `FillFactor` values, or the corresponding
  `Gridspace` when an explicit finite source is supplied.
"""
function FillFactor(η::Union{Real, AbstractGrid, Gridspace}; combine::Symbol = :product)
    return parameterize(FillFactor, FillFactor, (η,); combine)
end
FillFactor(factor::FillFactor; combine::Symbol = :product) = factor
function FillFactor(values::Union{Tuple, AbstractVector}; combine::Symbol = :product)
    return map(value -> FillFactor(value; combine), Tuple(values))
end
function FillFactor(first, second, remaining...; combine::Symbol = :product)
    return FillFactor((first, second, remaining...); combine)
end

# The compaction of repeated members: `nothing` keeps the natural circular members, a Bool
# selects whether a bounded formation compacts, and a fill factor sets the material-area
# fraction.
const Compaction = Union{Nothing, Bool, FillFactor}

# Compacted strips fill the declared area fraction of the course envelope.
function placements(pattern::Ring{Int, <:Real, <:Real}, item::Rectangle, factor::FillFactor)
    inner = max(zero(pattern.r), pattern.r - item.h / 2)
    outer = sqrt(
        inner^2 + 2pattern.n * area(item) / (pattern.span * factor.η)
    )
    angular_width = factor.η * pattern.span / pattern.n
    return _strip_placements(pattern, inner, outer, angular_width)
end

function capacity(pattern::Ring{<:Any, <:Real, <:Real}, item::Rectangle, factor::FillFactor)
    radial_half_extent = item.h / 2
    pattern.r > zero(pattern.r) || return 0
    factor.η <= inv(one(pattern.gap_frac) + pattern.gap_frac) || throw(DomainError(
        (factor.η, pattern.gap_frac),
        "fill factor cannot provide the requested member gap"
    ))
    inner = max(zero(pattern.r), pattern.r - radial_half_extent)
    outer = pattern.r + radial_half_extent
    envelope_area = pattern.span * (outer^2 - inner^2) / 2
    count = factor.η * envelope_area / (item.w * item.h)
    return max(0, floor(Int, count + 8eps(float(count))))
end
