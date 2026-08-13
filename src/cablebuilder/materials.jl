@gridspace @relax struct Material{T <: Real}
	"Electrical resistivity of the material \\[Ω·m\\]."
	rho::T
	"Relative permittivity \\[dimensionless\\]."
	eps_r::T
	"Relative permeability \\[dimensionless\\]."
	mu_r::T
	"Reference temperature for property evaluations \\[°C\\]."
	T0::T = 20.0
	"Temperature coefficient of resistivity \\[1/°C\\]."
	alpha::T = 0.0
	"Thermal resistivity \\[K·m/W\\]."
	rho_thermal::T = 0.0
	"Maximum operating temperature \\[°C\\]."
	theta_max::T = 90.0
	"Dielectric loss factor \\[dimensionless\\]."
	tan_delta::T = 0.0
	"Solar absorption coefficient \\[dimensionless\\]."
	sigma_solar::T = 0.0
end