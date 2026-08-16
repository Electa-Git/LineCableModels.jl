# Export public API
export f₀, μ₀, ε₀, ρ₀, T₀, TOL, ΔTmax
export BASE_FLOAT, REALSCALAR, COMPLEXSCALAR

# General constants
"Base power system frequency, f₀ = 50.0 [Hz]."
const f₀ = 50.0
"Magnetic constant (vacuum permeability), μ₀ = 4π * 1e-7 [H/m]."
const μ₀ = 4π * 1e-7
"Electric constant (vacuum permittivity), ε₀ = 8.8541878128e-12 [F/m]."
const ε₀ = 8.8541878128e-12
"Annealed copper reference resistivity, ρ₀ = 1.724e-08 [Ω·m]."
const ρ₀ = 1.724e-08
"Base temperature for conductor properties, T₀ = 20.0 [°C]."
const T₀ = 20.0
"Maximum tolerance for temperature variations, ΔTmax = 150 [°C]."
const ΔTmax = 150.0
"Default tolerance for floating-point comparisons, TOL = 1e-6."
const TOL = 1e-6

# Dependency-free numeric interfaces. Optional scalar packages such as
# Measurements.jl participate through ordinary `Real` promotion in extensions.
const BASE_FLOAT = Float64
const REALSCALAR = Real
const COMPLEXSCALAR = Complex{T} where {T <: Real}

# Preserve the source definitions when a generic numerical kernel requests more
# precision than the public Float64 constants carry.
@inline _typed_ε₀(::Type{T}) where {T <: Real} = T(ε₀)
@inline _typed_ε₀(::Type{BigFloat}) = parse(BigFloat, "8.8541878128e-12")
@inline _typed_μ₀(::Type{T}) where {T <: Real} = T(μ₀)
@inline _typed_μ₀(::Type{BigFloat}) = BigFloat(4) * BigFloat(π) * parse(BigFloat, "1e-7")
