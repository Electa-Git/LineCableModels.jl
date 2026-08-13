abstract type Homogeneous <: EarthImpedanceFormulation end

struct Kernel{Tγ1, Tγ2, Tμ2}
	"Layer where the source conductor is placed."
	s::Int
	"Layer where the target conductor is placed."
	t::Int
	"Primary field propagation constant (0 = lossless, 1 = air, 2 = earth)."
	Γx::Int
	"Air propagation constant γ₁(jω, μ, σ, ε)."
	γ1::Tγ1
	"Earth propagation constant γ₂(jω, μ, σ, ε)."
	γ2::Tγ2
	"Earth magnetic-constant assumption μ₂(μ)."
	μ2::Tμ2
end

struct Papadopoulos{Tγ1, Tγ2, Tμ2} <: Homogeneous
	kernel::Kernel{Tγ1, Tγ2, Tμ2}
end

Papadopoulos(; s::Int = 2, t::Int = 2, Γx::Int = 2,
	γ1 = (jω, μ, σ, ε) -> sqrt(jω * μ * (σ + jω*ε)),
	γ2 = (jω, μ, σ, ε) -> sqrt(jω * μ * (σ + jω*ε)),
	μ2 = μ -> μ) =
	Papadopoulos(
		Kernel{typeof(γ1), typeof(γ2), typeof(μ2)}(s, t, Γx, γ1, γ2, μ2),
	)

get_description(::Papadopoulos) = "Papadopoulos"
from_kernel(f::Papadopoulos) = f.kernel


struct Pollaczek{Tγ1, Tγ2, Tμ2} <: Homogeneous
	kernel::Kernel{Tγ1, Tγ2, Tμ2}
end

Pollaczek(; s::Int = 2, t::Int = 2, Γx::Int = 0,
	γ1 = (jω, μ, σ, ε) -> jω * sqrt(μ * ε),
	γ2 = (jω, μ, σ, ε) -> sqrt(jω * μ * σ),
	μ2 = μ -> oftype(μ, μ₀)) =
	Pollaczek(
		Kernel{typeof(γ1), typeof(γ2), typeof(μ2)}(s, t, Γx, γ1, γ2, μ2),
	)

get_description(::Pollaczek) = "Pollaczek"
from_kernel(f::Pollaczek) = f.kernel

struct Carson{Tγ1, Tγ2, Tμ2} <: Homogeneous
	kernel::Kernel{Tγ1, Tγ2, Tμ2}
end

Carson(; s::Int = 1, t::Int = 1, Γx::Int = 0,
	γ1 = (jω, μ, σ, ε) -> jω * sqrt(μ * ε),
	γ2 = (jω, μ, σ, ε) -> sqrt(jω * μ * σ),
	μ2 = μ -> oftype(μ, μ₀)) =
	Carson(
		Kernel{typeof(γ1), typeof(γ2), typeof(μ2)}(s, t, Γx, γ1, γ2, μ2),
	)

get_description(::Carson) = "Carson"
from_kernel(f::Carson) = f.kernel


# ρ, ε, μ = ws.rho_g, ws.eps_g, ws.mu_g
#     f(h, d, @view(ρ[:,k]), @view(ε[:,k]), @view(μ[:,k]), ws.freq[k])

# Functor implementation for all homogeneous earth impedance formulations.
function (f::Homogeneous)(
	form::Symbol,
	h::AbstractVector{T},
	yij::T,
	rho_g::AbstractVector{T},
	eps_g::AbstractVector{T},
	mu_g::AbstractVector{T},
	jω::Complex{T},
) where {T <: REALSCALAR}
	Base.@nospecialize form
	return form === :self ? f(Val(:self), h, yij, rho_g, eps_g, mu_g, jω) :
		   form === :mutual ? f(Val(:mutual), h, yij, rho_g, eps_g, mu_g, jω) :
		   throw(ArgumentError("Unknown earth impedance form: $form"))
end

# function (f::Homogeneous)(
# 	h::AbstractVector{T},
# 	yij::T,
# 	rho_g::AbstractVector{T},
# 	eps_g::AbstractVector{T},
# 	mu_g::AbstractVector{T},
# 	jω::Complex{T},
# ) where {T <: REALSCALAR}
# 	return f(Val(:mutual), h, yij, rho_g, eps_g, mu_g, jω)
# end

function (f::Homogeneous)(
	::Val{:self},
	h::AbstractVector{T},
	yij::T,
	rho_g::AbstractVector{T},
	eps_g::AbstractVector{T},
	mu_g::AbstractVector{T},
	jω::Complex{T},
) where {T <: REALSCALAR}
	return f(Val(:mutual), h, yij, rho_g, eps_g, mu_g, jω)
end

@inline _not(s::Int) =
	(s == 1 || s == 2) ? (3 - s) :
	throw(ArgumentError("s must be 1 or 2"))

@inline _get_layer(z) =
	z > 0 ? 1 :
	(z < 0 ? 2 : throw(ArgumentError("Conductor at interface (h=0) is invalid")))

@noinline function _layer_mismatch(which::AbstractString, got::Int, expected::Int)
	throw(
		ArgumentError(
			"conductor $which is in layer $got but formulation expects layer $expected",
		),
	)
end

@inline function validate_layers!(f::Homogeneous, h)
	@boundscheck length(h) == 2 || throw(ArgumentError("h must have length 2"))
	ℓ1 = _get_layer(h[1])
	ℓ2 = _get_layer(h[2])
	(ℓ1 == f.s) || _layer_mismatch("i (h[1])", ℓ1, f.s)
	(ℓ2 == f.t) || _layer_mismatch("j (h[2])", ℓ2, f.t)
	return nothing
end

@inline function (f::Homogeneous)(
	::Val{:mutual},
	h::AbstractVector{T},
	yij::T,
	rho_g::AbstractVector{T},
	eps_g::AbstractVector{T},
	mu_g::AbstractVector{T},
	jω::Complex{T},
) where {T <: REALSCALAR}

	validate_layers!(f, h)

	s = f.s # index of source layer
	o = _not(s) # the other layer
	nL = length(rho_g)
	μ = similar(mu_g);
	σ = similar(rho_g);
	@inbounds for i in 1:nL
		μ[i] = (i == 1) ? mu_g[i] : f.μ2(mu_g[i]) # μ₂ for earth layers
		σ[i] = _to_σ(rho_g[i])
	end

	# construct propagation constants according to formulation assumptions
	γ = Vector{Complex{T}}(undef, nL)
	@inbounds for i in 1:nL
		γ[i] = (i == 1 ? f.γ1 : f.γ2)(jω, μ[i], σ[i], eps_g[i])
	end
	γ_s = γ[s];
	γ_o = γ[o]
	γs_2 = γ_s^2
	γo_2 = γ_o^2

	# kx from struct: 0:none, 1:air, 2:source layer
	kx_2 = if f.Γx == 0 # precalc squared
		zero(γs_2)
	else
		ℓ = (f.Γx == 1) ? 1 : s
		oftype(γs_2, (-jω^2) * μ[ℓ] * eps_g[ℓ])
	end

	# unpack geometry
	@inbounds hi, hj = abs(h[1]), abs(h[2])
	dij = hypot(yij, hi - hj)          # √(y^2 + (hi - hj)^2) - conductor-conductor
	Dij = hypot(yij, hi + hj)          # √(y^2 + (hi + hj)^2) - conductor-image

	# perfectly conducting earth term in Bessel form
	Λij = _bessel_diff(γ_s, dij, Dij)

	# precompute scalars for integrand
	μ_s = μ[s]
	μ_o = μ[o]
	H = hi + hj

	# Sij = 2 ∫_0^∞ Fij(λ) cos(yij λ) dλ
	# integrand = (λ) -> Fij(λ) * cos(yij * λ)
	@inline function integrand(λ::Float64)::Complex{T}
		as = sqrt(λ*λ + γs_2 + kx_2)
		ao = sqrt(λ*λ + γo_2 + kx_2)

		F = μ_o * exp(-as*H) / (as*μ_o + ao*μ_s)

		F * cos(yij*λ)
	end

	Sij, _ = quadgk(
		integrand,
		0.0,
		1.0;
		rtol = 1e-8,
		norm = z -> abs(complex(value(real(z)), value(imag(z)))),
	)
	Sij *= 2

	return (jω * μ_s / (2π)) * (Λij + Sij)
end
