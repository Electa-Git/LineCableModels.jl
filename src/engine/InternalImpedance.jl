"""
	LineCableModels.Engine.InternalImpedance

# Dependencies

$(IMPORTS)

# Exports

$(EXPORTS)
"""
module InternalImpedance

# Export public API
export ScaledBessel

# Load common dependencies
using ...LineCableModels
include("../utils/commondeps.jl")

# Module-specific dependencies
using Measurements
using ...Utils
import ...LineCableModels: get_description
import ..Engine: InternalImpedanceFormulation

struct ScaledBessel <: InternalImpedanceFormulation end
get_description(::ScaledBessel) = "Scaled Bessel (Schelkunoff)"


macro uncertain_bessel(expr::Expr)
    f = esc(expr.args[1]) # Function name
    order = expr.args[2]  # First argument (order), no need to escape this
    a = esc(expr.args[3]) # Second argument (complex number with uncertainties)

    return :(Measurements.result(
        $f($order, Measurements.value($a)),
        vcat(
            Calculus.gradient(
                x -> real($f($order, complex(x[1], x[2]))),
                [reim(Measurements.value($a))...],
            ),
            Calculus.gradient(
                x -> imag($f($order, complex(x[1], x[2]))),
                [reim(Measurements.value($a))...],
            ),
        ),
        $a,
    ))
end

function calc_outer_skin_effect_impedance(
    radius_ex::Measurement{T},
    radius_in::Measurement{T},
    sigma_c::Measurement{T},
    mur_c::Measurement{T},
    f::T;
    SimplifiedFormula=false,
) where {T<:Real}

    # Constants
    m0 = 4 * pi * 1e-7
    mu_c = m0 * mur_c
    TOL = 1e-6
    omega = 2 * pi * f

    # Calculate the reciprocal of the skin depth
    m = sqrt(im * omega * mu_c * sigma_c)

    # Approximated skin effect
    if radius_in == 0
        radius_in = eps()  # Avoid division by zero
    end

    if SimplifiedFormula
        if radius_in < TOL
            cothTerm = coth(m * radius_ex * 0.733)
        else
            cothTerm = coth(m * (radius_ex - radius_in))
        end
        Z1 = (m / sigma_c) / (2 * pi * radius_ex) * cothTerm

        if radius_in < TOL
            Z2 = 0.3179 / (sigma_c * pi * radius_ex^2)
        else
            Z2 = 1 / (sigma_c * 2 * pi * radius_ex * (radius_in + radius_ex))
        end
        zin = Z1 + Z2
    else
        # More detailed solution with Bessel functions and uncertainty
        w_out = m * radius_ex
        w_in = m * radius_in

        s_in = exp(abs(real(w_in)) - w_out)
        s_out = exp(abs(real(w_out)) - w_in)
        sc = s_in / s_out  # Should be applied to all besseli() involving w_in

        # Bessel function terms with uncertainty handling using the macro
        N =
            @uncertain_bessel(besselix(0, w_out)) * @uncertain_bessel(besselkx(1, w_in)) +
            sc *
            @uncertain_bessel(besselkx(0, w_out)) *
            @uncertain_bessel(besselix(1, w_in))

        D =
            @uncertain_bessel(besselix(1, w_out)) * @uncertain_bessel(besselkx(1, w_in)) -
            sc *
            @uncertain_bessel(besselkx(1, w_out)) *
            @uncertain_bessel(besselix(1, w_in))

        # Final impedance calculation
        zin = (im * omega * mu_c / (2 * pi)) * (1 / w_out) * (N / D)
    end

    return zin
end

function calc_inner_skin_effect_impedance(
    radius_ex::Measurement{T},
    radius_in::Measurement{T},
    sigma_c::Measurement{T},
    mur_c::Measurement{T},
    f::T;
    SimplifiedFormula=false,
) where {T<:Real}

    # Constants
    m0 = 4 * pi * 1e-7
    mu_c = m0 * mur_c
    omega = 2 * pi * f

    # Calculate the reciprocal of the skin depth
    m = sqrt(im * omega * mu_c * sigma_c)

    # Approximated skin effect
    if radius_in == 0
        radius_in = eps()  # Avoid division by zero
    end

    if SimplifiedFormula
        cothTerm = coth(m * (radius_ex - radius_in))
        Z1 = (m / sigma_c) / (2 * pi * radius_in) * cothTerm
        Z2 = 1 / (2 * pi * radius_in * (radius_in + radius_ex) * sigma_c)
        zin = Z1 + Z2
    else
        # More detailed solution with Bessel functions and uncertainty
        w_out = m * radius_ex
        w_in = m * radius_in

        s_in = exp(abs(real(w_in)) - w_out)
        s_out = exp(abs(real(w_out)) - w_in)
        sc = s_in / s_out  # Should be applied to all besselix() involving w_in

        # Bessel function terms with uncertainty handling using the macro
        N =
            sc *
            (@uncertain_bessel besselix(0, w_in)) *
            (@uncertain_bessel besselkx(1, w_out)) +
            (@uncertain_bessel besselkx(0, w_in)) * (@uncertain_bessel besselix(1, w_out))

        D =
            (@uncertain_bessel besselix(1, w_out)) * (@uncertain_bessel besselkx(1, w_in)) -
            sc *
            (@uncertain_bessel besselkx(1, w_out)) *
            (@uncertain_bessel besselix(1, w_in))

        # Final impedance calculation
        zin = (im * omega * mu_c / (2 * pi)) * (1 / w_in) * (N / D)
    end

    return zin
end

function calc_mutual_skin_effect_impedance(
    radius_ex::Measurement{T},
    radius_in::Measurement{T},
    sigma_c::Measurement{T},
    mur_c::Measurement{T},
    f::T;
    SimplifiedFormula=false,
) where {T<:Real}

    # Constants
    m0 = 4 * pi * 1e-7
    mu_c = m0 * mur_c
    omega = 2 * pi * f

    # Calculate the reciprocal of the skin depth
    m = sqrt(im * omega * mu_c * sigma_c)

    # Approximated skin effect
    if radius_in == 0
        radius_in = eps()  # Avoid division by zero
    end

    if SimplifiedFormula
        cschTerm = csch(m * (radius_ex - radius_in))
        zm = m / (sigma_c * pi * (radius_in + radius_ex)) * cschTerm
    else
        # More detailed solution with Bessel functions and uncertainty
        w_out = m * radius_ex
        w_in = m * radius_in

        s_in = exp(abs(real(w_in)) - w_out)
        s_out = exp(abs(real(w_out)) - w_in)
        sc = s_in / s_out  # Should be applied to all besselix() involving w_in

        # Bessel function terms with uncertainty handling using the macro
        D =
            (@uncertain_bessel besselix(1, w_out)) * (@uncertain_bessel besselkx(1, w_in)) -
            sc *
            (@uncertain_bessel besselkx(1, w_out)) *
            (@uncertain_bessel besselix(1, w_in))

        # Final mutual impedance calculation
        zm = 1 / (2 * pi * radius_ex * radius_in * sigma_c * D * s_out)
    end

    return zm
end

# # This is a hypothetical internal function
# function _compute_earth_impedance_matrix!(
#     Z_earth::Matrix,
#     problem::CableSystem,
#     formulation::EarthImpedanceFormulation, # Dispatches on the abstract type!
#     ω::Float64
# )
#     ρ_earth = problem.soil.resistivity

#     for i in 1:problem.num_conductors
#         for j in 1:i
#             # Get conductor positions
#             h_i = problem.conductors[i].y
#             h_j = problem.conductors[j].y
#             x_ij = abs(problem.conductors[i].x - problem.conductors[j].x)

#             # Call the formulation object directly!
#             # Julia will dispatch to the correct implementation (Carson, Deri, etc.)
#             Z_earth[i, j] = formulation(ω, ρ_earth, h_i, h_j, x_ij)
#             Z_earth[j, i] = Z_earth[i, j] # Assuming symmetry
#         end
#     end

#     return Z_earth
# end

# struct AnalyticalFormulation <: AbstractFormulation
#     # --- Fields for different calculation parts ---
#     internal_impedance::InternalImpedanceFormulation
#     earth_return::EarthImpedanceFormulation
#     # You could add more here, e.g., admittance_model, temperature_correction, etc.

#     # --- Keyword constructor for usability ---
#     function AnalyticalFormulation(;
#         internal_impedance::InternalImpedanceFormulation,
#         earth_return::EarthImpedanceFormulation
#     )
#         new(internal_impedance, earth_return)
#     end
# end

# function compute!(problem::CableSystem, formulation::AnalyticalFormulation, opts::Options, ω::Float64)

#     # 1. Compute Internal Impedance
#     # The specific method (e.g., rigorous Bessel vs. approximate) is determined
#     # by the type of formulation.internal_impedance
#     Z_internal = _compute_internal_impedance_matrix!(
#         problem, 
#         formulation.internal_impedance, # Pass the specific object
#         ω
#     )

#     # 2. Compute Earth-Return Impedance
#     # The specific method (Carson vs. Deri) is determined by the
#     # type of formulation.earth_return
#     Z_earth = _compute_earth_impedance_matrix!(
#         problem, 
#         formulation.earth_return, # Pass the specific object
#         ω
#     )

#     # 3. Combine results
#     Z_phase = Z_internal + Z_earth

#     # ... compute Y, sequence components, etc. ...

#     return line_parameters
# end

end # module InternalImpedance

