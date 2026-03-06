
abstract type AbstractLayout end
struct Concentric <: AbstractLayout end
struct SectorShaped <: AbstractLayout end

abstract type AbstractShape{L <: AbstractLayout, T <: Real} end

import .Validation:
	Validation, required_fields, extra_rules, Normalized, Finite, Nonneg, Less, Positive


# Global accessors
@inline r_in(s::AbstractShape) = s.r_in
@inline r_ex(s::AbstractShape) = s.r_ex

include("solidcore.jl")
include("tubular.jl")
include("enclosure.jl")
include("circstranded.jl")
