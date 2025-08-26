"""
	LineCableModels.Engine.EarthImpedance

# Dependencies

$(IMPORTS)

# Exports

$(EXPORTS)
"""
module EarthImpedance

# Export public API
export Papadopoulos

# Load common dependencies
using ...LineCableModels
include("../utils/commondeps.jl")

# Module-specific dependencies
using Measurements
using ...Utils
import ...LineCableModels: get_description
import ..Engine: EarthImpedanceFormulation

struct Papadopoulos <: EarthImpedanceFormulation end
get_description(::Papadopoulos) = "Papadopoulos (homogeneous earth)"

struct Carson <: EarthImpedanceFormulation end

function calc_self_impedance_papadopoulos(
    h::Vector{Measurement{Float64}},
    r::Vector{Measurement{Float64}},
    eps_g::Measurement{Float64},
    mu_g::Measurement{Float64},
    sigma_g::Measurement{Float64},
    f::Float64,
    con::Int,
    kx::Int=0,
)

    # Constants
    sig0 = 0.0
    eps0 = 8.8541878128e-12
    mu0 = 4 * pi * 1e-7
    w = 2 * pi * f  # Angular frequency

    # Define k_x based on input kx type
    # 0 = neglect propagation constant
    # 1 = use value of layer 1 (air)
    # 2 = use value of layer 2 (earth)
    k_x = if kx == 2
        ω -> ω * sqrt(mu_g * (eps_g - im * (sigma_g / ω)))
    elseif kx == 1
        ω -> ω * sqrt(mu0 * eps0)
    else
        ω -> ω * 0.0  # Default to zero
    end

    # Define gamma_0 and gamma_1
    gamma_0 = ω -> sqrt(im * ω * mu0 * (sig0 + im * ω * eps0))
    gamma_1 = ω -> sqrt(im * ω * mu_g * (sigma_g + im * ω * eps_g))

    # Define a_0 and a_1
    a_0 = (λ, ω) -> sqrt(λ^2 + gamma_0(ω)^2 + k_x(ω)^2)
    a_1 = (λ, ω) -> sqrt(λ^2 + gamma_1(ω)^2 + k_x(ω)^2)

    # Initialize Zg_self matrix (complex numbers)
    Zg_self = zeros(Complex{Measurement{Float64}}, con, con)

    for k ∈ 1:con
        if h[k] < 0  # Only process if h(k) < 0 (as per the original MATLAB logic)

            # Define the function zz
            zz =
                (a0, a1, hi, hj, λ, mu0, mu_g, ω, y) -> (
                    (mu_g * ω * exp(-a1 * abs(hi - hj + 1e-3)) * cos(λ * y) * 0.5im) /
                    (a1 * pi) -
                    (
                        mu_g *
                        ω *
                        exp(-a1 * (hi - hj)) *
                        cos(λ * y) *
                        (a0 * mu_g + a1 * mu0 * sign(hi)) *
                        0.5im
                    ) / (a1 * pi * (a0 * mu_g + a1 * mu0))
                )

            # Define zfun based on lambda and omega (as in the MATLAB code)
            zfun = λ -> begin
                a0 = a_0(λ, w)
                a1 = a_1(λ, w)
                zz(a0, a1, h[k], h[k], λ, mu0, mu_g, w, r[k])
            end

            # Perform the numerical integration (over complex numbers)
            Js, _ = quadgk(zfun, 0.0, Inf; rtol=1e-6)

            # Store the result (which is complex) in Zg_self
            Zg_self[k, k] = Js
        end
    end

    return Zg_self
end

"""
	calc_mutual_impedance_papadopoulos(h, d, eps_g, mu_g, sigma_g, f, con; kx=0)

Calculates the mutual earth impedance between conductors using the Papadopoulos formula,
considering the properties of the ground and the frequency of the system.

# Arguments
- `h`: Vector of vertical distances to ground for conductors.
- `d`: Matrix of distances between conductors.
- `eps_g`: Relative permittivity of the earth.
- `mu_g`: Permeability of the earth.
- `sigma_g`: Conductivity of the earth.
- `f`: Frequency of the system (in Hz).
- `con`: Number of conductors in the system.
- `kx`: An optional flag for selecting propagation constant:
  - `0`: Default, no propagation constant.
  - `1`: Propagation constant for air.
  - `2`: Propagation constant for earth.

# Returns
- `Matrix{Complex{Measurement{Float64}}}`: The mutual earth impedance matrix for the given conductors.

# Category: Earth-return impedances and admittances

"""
function calc_mutual_impedance_papadopoulos(
    h::Vector{Measurement{Float64}},
    d::Matrix{Measurement{Float64}},
    eps_g::Measurement{Float64},
    mu_g::Measurement{Float64},
    sigma_g::Measurement{Float64},
    f::Float64,
    con::Int,
    kx::Int=0,
)

    # Constants
    sig0 = 0.0
    eps0 = 8.8541878128e-12
    mu0 = 4 * pi * 1e-7
    w = 2 * pi * f  # Angular frequency

    # Define k_x based on input kx type
    # 0 = neglect propagation constant
    # 1 = use value of layer 1 (air)
    # 2 = use value of layer 2 (earth)
    k_x = if kx == 2
        ω -> ω * sqrt(mu_g * (eps_g - im * (sigma_g / ω)))
    elseif kx == 1
        ω -> ω * sqrt(mu0 * eps0)
    else
        ω -> ω * 0.0  # Default to zero
    end

    # Define gamma_0 and gamma_1
    gamma_0 = ω -> sqrt(im * ω * mu0 * (sig0 + im * ω * eps0))
    gamma_1 = ω -> sqrt(im * ω * mu_g * (sigma_g + im * ω * eps_g))

    # Define a_0 and a_1
    a_0 = (λ, ω) -> sqrt(λ^2 + gamma_0(ω)^2 + k_x(ω)^2)
    a_1 = (λ, ω) -> sqrt(λ^2 + gamma_1(ω)^2 + k_x(ω)^2)

    # Initialize Zg_mutual matrix (complex numbers)
    Zg_mutual = zeros(Complex{Measurement{Float64}}, con, con)

    # Mutual Impedance
    for x ∈ 1:con
        for y ∈ x+1:con
            if x != y
                h1 = h[x]
                h2 = h[y]

                if h1 < 0 && h2 < 0
                    # Define the function zz
                    zz =
                        (a0, a1, hi, hj, λ, mu0, mu_g, ω, y) -> (
                            (
                                mu_g *
                                ω *
                                exp(-a1 * abs(hi - hj + 1e-3)) *
                                cos(λ * y) *
                                0.5im
                            ) / (a1 * pi) -
                            (
                                mu_g *
                                ω *
                                exp(a1 * (hi + hj)) *
                                cos(λ * y) *
                                (a0 * mu_g + a1 * mu0 * sign(hi)) *
                                0.5im
                            ) / (a1 * pi * (a0 * mu_g + a1 * mu0))
                        )

                    # Define zfun based on lambda and omega
                    zfun = λ -> begin
                        a0 = a_0(λ, w)
                        a1 = a_1(λ, w)
                        zz(a0, a1, h1, h2, λ, mu0, mu_g, w, d[x, y])
                    end

                    # Perform the numerical integration (over complex numbers)
                    Jm, _ = quadgk(zfun, 0.0, Inf; rtol=1e-6)

                    # Store the result (which is complex) in Zg_mutual
                    Zg_mutual[x, y] = Jm
                    Zg_mutual[y, x] = Zg_mutual[x, y]
                end
            end
        end
    end

    return Zg_mutual
end



# This is the "functor" part. We're making the Carson object callable.
# It takes the arguments needed for ANY earth impedance calculation.
function (f::Carson)(ω::Float64, ρ_earth::Float64, h_i::Float64, h_j::Float64, x_ij::Float64)
    # ... implementation of Carson's formula here ...
    # This function will compute and return the self/mutual impedance Z_ij
    # based on the inputs.
    return z_ij
end

end # module EarthImpedance