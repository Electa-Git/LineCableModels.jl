"""
$(TYPEDEF)

Abstract marker for an earth-return method implemented by an external
calculation backend but not by the built-in analytical solver.
"""
abstract type ReferenceEarthImpedance <: EarthImpedanceFormulation end

"""
$(TYPEDEF)

Select the Deri-Semlyen overhead-line approximation.
"""
struct Deri <: ReferenceEarthImpedance end

"""
$(TYPEDEF)

Select the Wedepohl underground-cable approximation.
"""
struct Wedepohl <: ReferenceEarthImpedance end

"""
$(TYPEDEF)

Select the Saad underground-cable approximation.
"""
struct Saad <: ReferenceEarthImpedance end

"""
$(TYPEDEF)

Select the Ametani aerial-to-underground mutual-impedance method.
"""
struct Ametani <: ReferenceEarthImpedance end

"""
$(TYPEDEF)

Select the Lucca aerial-to-underground mutual-impedance method.
"""
struct Lucca <: ReferenceEarthImpedance end

"""
$(TYPEDEF)

Select direct numerical integration for `Placement`, which is either
`:overhead` or `:underground`.
"""
struct DirectNumericalIntegration{Placement} <: ReferenceEarthImpedance
    function DirectNumericalIntegration(placement::Symbol)
        placement in (:overhead, :underground) || throw(ArgumentError(
            "direct numerical integration placement must be :overhead or :underground",
        ))
        return new{placement}()
    end
end

description(::Deri) = "Deri-Semlyen"
description(::Wedepohl) = "Wedepohl"
description(::Saad) = "Saad"
description(::Ametani) = "Ametani"
description(::Lucca) = "Lucca"
function description(::DirectNumericalIntegration{:overhead})
    "Direct numerical integration (overhead)"
end
function description(::DirectNumericalIntegration{:underground})
    "Direct numerical integration (underground)"
end

function (formulation::ReferenceEarthImpedance)(args...)
    throw(ErrorException(
        "$(description(formulation)) earth impedance is not implemented for the analytical backend",
    ))
end
