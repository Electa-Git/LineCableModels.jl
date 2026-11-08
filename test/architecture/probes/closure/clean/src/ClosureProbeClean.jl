module ClosureProbeClean

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

# One method per member.
measure(curve::Segment) = curve.p
measure(curve::Arc) = curve.r
measure(curve::EllipticArc) = curve.a

# The method on the union throws for every pair other than two arcs.
pair(a::Arc, b::Arc) = a.r + b.r
pair(a::Curve, b::Curve) = throw(DomainError((a, b), "a pair needs two arcs"))

# A Base generic with one method on the union.
Base.reverse(curve::Curve) = curve

end
