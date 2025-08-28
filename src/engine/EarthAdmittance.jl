"""
	LineCableModels.Engine.EarthAdmittance

# Dependencies

$(IMPORTS)

# Exports

$(EXPORTS)
"""
module EarthAdmittance

# Export public API
# export 

# Load common dependencies
using ...LineCableModels
include("../utils/commondeps.jl")

# Module-specific dependencies
using Measurements
using ...Utils
import ...LineCableModels: get_description
import ..Engine: EarthAdmittanceFormulation

struct Papadopoulos <: EarthAdmittanceFormulation end
get_description(::Papadopoulos) = "Papadopoulos (homogeneous earth)"

function calc_self_potential_coeff_papadopoulos(
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
    w = 2 * pi * f

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

    # Define gamma and a functions
    gamma_0 = ω -> sqrt(im * ω * mu0 * (sig0 + im * ω * eps0))
    gamma_1 = ω -> sqrt(im * ω * mu_g * (sigma_g + im * ω * eps_g))
    a_0 = (λ, ω) -> sqrt(λ^2 + gamma_0(ω)^2 + k_x(ω)^2)
    a_1 = (λ, ω) -> sqrt(λ^2 + gamma_1(ω)^2 + k_x(ω)^2)

    Pg_self = zeros(Complex{Measurement{Float64}}, con, con)
    TOL = 1e-3

    for k ∈ 1:con
        if h[k] < 0
            yy =
                (a0, a1, gamma0, gamma1, hi, hj, λ, mu0, mu_g, ω, y) ->
                    (
                        1.0 / gamma1^2 *
                        mu_g *
                        ω *
                        exp(-a1 * abs(hi - hj + TOL)) *
                        cos(λ * y) *
                        0.5im
                    ) / (a1 * pi) -
                    (
                        1.0 / gamma1^2 *
                        mu_g *
                        ω *
                        exp(a1 * (hi + hj)) *
                        cos(λ * y) *
                        (a0 * mu_g + a1 * mu0 * sign(hi)) *
                        0.5im
                    ) / (a1 * pi * (a0 * mu_g + a1 * mu0)) +
                    (
                        a1 * 1.0 / gamma1^2 *
                        mu0 *
                        mu_g^2 *
                        ω *
                        exp(a1 * (hi + hj)) *
                        cos(λ * y) *
                        (sign(hi) - 1.0) *
                        (gamma0^2 - gamma1^2) *
                        0.5im
                    ) / (
                        pi *
                        (a0 * gamma1^2 * mu0 + a1 * gamma0^2 * mu_g) *
                        (a0 * mu_g + a1 * mu0)
                    )

            yfun =
                λ -> yy(
                    a_0(λ, w),
                    a_1(λ, w),
                    gamma_0(w),
                    gamma_1(w),
                    h[k],
                    h[k],
                    λ,
                    mu0,
                    mu_g,
                    w,
                    r[k],
                )
            Qs, _ = quadgk(yfun, 0, Inf, rtol=1e-6)
            Pg_self[k, k] = (im * w * Qs)
        end
    end

    return Pg_self
end

"""
	calc_mutual_potential_coeff_papadopoulos(h, d, eps_g, mu_g, sigma_g, f, con; kx=0)

Calculates the mutual earth potential coefficient between conductors using the Papadopoulos
formula, considering the properties of the ground and the frequency of the system.

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
- `Matrix{Complex{Measurement{Float64}}}`: The mutual earth potential coefficient matrix for the given conductors.

# Category: Earth-return impedances and admittances

"""
function calc_mutual_potential_coeff_papadopoulos(
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
    w = 2 * pi * f

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

    # Define gamma and a functions
    gamma_0 = ω -> sqrt(im * ω * mu0 * (sig0 + im * ω * eps0))
    gamma_1 = ω -> sqrt(im * ω * mu_g * (sigma_g + im * ω * eps_g))
    a_0 = (λ, ω) -> sqrt(λ^2 + gamma_0(ω)^2 + k_x(ω)^2)
    a_1 = (λ, ω) -> sqrt(λ^2 + gamma_1(ω)^2 + k_x(ω)^2)

    Pg_mutual = zeros(Complex{Measurement{Float64}}, con, con)
    TOL = 1e-3

    # Mutual potential coefficient
    for x ∈ 1:con
        for y ∈ x+1:con
            if x != y
                h1 = h[x]
                h2 = h[y]
                if abs(h2 - h1) < TOL
                    h2 += TOL
                end

                if h1 < 0 && h2 < 0
                    yy =
                        (a0, a1, gamma0, gamma1, hi, hj, λ, mu0, mu_g, ω, y) ->
                            (
                                1.0 / gamma1^2 *
                                mu_g *
                                ω *
                                exp(-a1 * abs(hi - hj)) *
                                cos(λ * y) *
                                0.5im
                            ) / (a1 * pi) -
                            (
                                1.0 / gamma1^2 *
                                mu_g *
                                ω *
                                exp(a1 * (hi + hj)) *
                                cos(λ * y) *
                                (a0 * mu_g + a1 * mu0 * sign(hi)) *
                                0.5im
                            ) / (a1 * pi * (a0 * mu_g + a1 * mu0)) +
                            (
                                a1 * 1.0 / gamma1^2 *
                                mu0 *
                                mu_g^2 *
                                ω *
                                exp(a1 * (hi + hj)) *
                                cos(λ * y) *
                                (sign(hi) - 1.0) *
                                (gamma0^2 - gamma1^2) *
                                0.5im
                            ) / (
                                pi *
                                (a0 * gamma1^2 * mu0 + a1 * gamma0^2 * mu_g) *
                                (a0 * mu_g + a1 * mu0)
                            )

                    yfun =
                        λ -> yy(
                            a_0(λ, w),
                            a_1(λ, w),
                            gamma_0(w),
                            gamma_1(w),
                            h1,
                            h2,
                            λ,
                            mu0,
                            mu_g,
                            w,
                            d[x, y],
                        )
                    Qm, _ = quadgk(yfun, 0, Inf, rtol=1e-3)
                    Pg_mutual[x, y] = (im * w * Qm)
                    Pg_mutual[y, x] = Pg_mutual[x, y]
                end
            end
        end
    end

    return Pg_mutual
end

end # module EarthAdmittance