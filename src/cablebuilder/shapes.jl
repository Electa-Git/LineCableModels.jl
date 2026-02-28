
abstract type AbstractLayout end
struct Concentric <: AbstractLayout end
struct SectorShaped <: AbstractLayout end

abstract type AbstractShape{L <: AbstractLayout, T <: Real} end

import .Validation:
	Validation, required_fields, extra_rules, Normalized, Finite, Nonneg, Less, Positive

include("solidcore.jl")
include("tubular.jl")
include("enclosure.jl")
