module ClosureProbePlanted

struct Segment{T <: Real}
    p::T
end

struct Arc{T <: Real}
    r::T
end

struct EllipticArc{T <: Real}
    a::T
end

const Curve{T} = Union{Segment{T}, Arc{T}, EllipticArc{T}}

# `Arc` and `EllipticArc` lack a method of this owned generic.
measure(curve::Segment) = curve.p

# Every other member lacks a method at both positions of this pair.
pair(a::Arc, b::Arc) = a.r + b.r

# The elliptic arc lacks a method of this Base generic.
Base.reverse(curve::Segment) = curve
Base.reverse(curve::Arc) = curve

# The elliptic arc lacks a method beside this union of two members.
extent(curve::Union{Segment, Arc}) = zero(Float64)

end
