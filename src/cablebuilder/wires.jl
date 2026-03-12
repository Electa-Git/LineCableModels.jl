abstract type AbstractWire{T <: Real} end

# ==========================================
# CIRCULAR WIRE
# ==========================================
struct CircularWire{T <: Real} <: AbstractWire{T}
	r::T
end

CircularWire(r::T) where {T <: Real} = CircularWire{T}(r)
Base.convert(::Type{<:AbstractWire{T}}, c::CircularWire) where {T <: Real} =
	CircularWire{T}(convert(T, c.r))

struct CircularWireBuilder{T <: Real}
	r::T
end
@inline (b::CircularWireBuilder)() = CircularWire(b.r)

struct CircularWireSpec{T} <: AbstractSpec{CircularWireBuilder}
	r::T
end

# ==========================================
# RECTANGULAR WIRE
# ==========================================
struct RectangularWire{T <: Real} <: AbstractWire{T}
	w::T
	h::T
end

function RectangularWire(w, h)
	T = promote_type(typeof(w), typeof(h))
	return RectangularWire{T}(convert(T, w), convert(T, h))
end

Base.convert(::Type{<:AbstractWire{T}}, r::RectangularWire) where {T <: Real} =
	RectangularWire{T}(convert(T, r.w), convert(T, r.h))

struct RectangularWireBuilder{W <: Real, H <: Real}
	w::W
	h::H
end

# If you pass w and h backwards in your DSL, congratulations, your cable is now 
# physically impossible and your FEM mesh will look like a Picasso painting.
@inline (b::RectangularWireBuilder)() = RectangularWire(b.w, b.h)

struct RectangularWireSpec{W, H} <: AbstractSpec{RectangularWireBuilder}
	w::W
	h::H
end