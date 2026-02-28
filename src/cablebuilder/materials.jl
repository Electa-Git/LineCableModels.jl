
struct Material{T <: Real}
	"Electrical resistivity of the material \\[Ω·m\\]."
	rho::T
	"Relative permittivity \\[dimensionless\\]."
	eps_r::T
	"Relative permeability \\[dimensionless\\]."
	mu_r::T
	"Reference temperature for property evaluations \\[°C\\]."
	T0::T
	"Temperature coefficient of resistivity \\[1/°C\\]."
	alpha::T
	"Thermal resistivity \\[K·m/W\\]."
	rho_thermal::T
	"Maximum operating temperature \\[°C\\]."
	theta_max::T
	"Dielectric loss factor \\[dimensionless\\]."
	tan_delta::T
	"Solar absorption coefficient \\[dimensionless\\]."
	sigma_solar::T
end

function Material(
	rho::Real, eps_r::Real, mu_r::Real, T0::Real, alpha::Real,
	rho_thermal::Real, theta_max::Real, tan_delta::Real, sigma_solar::Real,
)
	# Promote all 9 inputs to a single common type T (e.g., mixed Float64 and Measurement)
	p = promote(rho, eps_r, mu_r, T0, alpha, rho_thermal, theta_max, tan_delta, sigma_solar)

	T_common = typeof(first(p))

	# Splat the perfectly aligned tuple into the strict, auto-generated constructor
	return Material{T_common}(p...)
end

import Base: convert

# Convert Material to a new precision T
function convert(::Type{Material{T}}, m::Material) where {T <: Real}
	return Material{T}(
		convert(T, m.rho),
		convert(T, m.eps_r),
		convert(T, m.mu_r),
		convert(T, m.T0),
		convert(T, m.alpha),
		convert(T, m.rho_thermal),
		convert(T, m.theta_max),
		convert(T, m.tan_delta),
		convert(T, m.sigma_solar),
	)
end

struct MaterialSpec{R, E, M, T0, A, RT, TM, TD, SS} <: AbstractSpec{Material}
	rho::R
	eps_r::E
	mu_r::M
	T0::T0
	alpha::A
	rho_thermal::RT
	theta_max::TM
	tan_delta::TD
	sigma_solar::SS
end

function MaterialSpec(;
	rho, eps_r, mu_r, T0 = 20.0, alpha = 0.0,
	rho_thermal = 0.0, theta_max = 90.0, tan_delta = 0.0, sigma_solar = 0.0,
)
	return MaterialSpec(
		Grid(rho), Grid(eps_r), Grid(mu_r), Grid(T0), Grid(alpha),
		Grid(rho_thermal), Grid(theta_max), Grid(tan_delta), Grid(sigma_solar),
	)
end

# The constructor to build a Spec directly from a materialized object
function MaterialSpec(m::Material)
	return MaterialSpec(
		Grid(m.rho), Grid(m.eps_r), Grid(m.mu_r), Grid(m.T0), Grid(m.alpha),
		Grid(m.rho_thermal), Grid(m.theta_max), Grid(m.tan_delta), Grid(m.sigma_solar),
	)
end

# The Diplomat (The exact tool for the job)
Base.convert(::Type{AbstractSpec{Material}}, m::Material) = MaterialSpec(m)
Base.convert(::Type{AbstractSpec{Material}}, m::AbstractSpec{Material}) = m