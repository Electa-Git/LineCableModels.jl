"""
$(TYPEDEF)

Select the Deri-Semlyen overhead-line earth-impedance formulation.
"""
struct Deri <: EarthImpedanceFormulation end

"""
$(TYPEDEF)

Select the Wedepohl underground-cable earth-impedance formulation.
"""
struct Wedepohl <: EarthImpedanceFormulation end

"""
$(TYPEDEF)

Select the Saad underground-cable earth-impedance formulation.
"""
struct Saad <: EarthImpedanceFormulation end

"""
$(TYPEDEF)

Select the Ametani aerial-to-underground mutual-impedance formulation.
"""
struct Ametani <: EarthImpedanceFormulation end

"""
$(TYPEDEF)

Select the Lucca aerial-to-underground mutual-impedance formulation.
"""
struct Lucca <: EarthImpedanceFormulation end

description(::Deri) = "Deri-Semlyen"
description(::Wedepohl) = "Wedepohl"
description(::Saad) = "Saad"
description(::Ametani) = "Ametani"
description(::Lucca) = "Lucca"
