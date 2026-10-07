"""
$(TYPEDSIGNATURES)

Identify the Wedepohl-Wilcox analytical approximations for the inner, outer,
and transfer surface impedances of a round conductor, in Ω/m. A solid cylinder
requires only the outer coefficient. An annulus uses all three coefficients.

This scientific identity is registered for backend dispatch. The owned coaxial backend has
no expression for it, so selecting it there fails before evaluation.
PSCAD's line-constants program uses these approximations for conductor surfaces.

Reference: L. M. Wedepohl and D. J. Wilcox, “Transient Analysis of Underground
Power-Transmission Systems: System-Model and Wave-Propagation Characteristics,”
*Proceedings of the IEE*, 120, 253–260, 1973. DOI: 10.1049/piee.1973.0056.
PSCAD 5.1 help reproduces the surface approximations in *Deriving System Y and Z
Matrices*, Eqs. (8-8), (8-9), and (8-11)-(8-13).
"""
function description(::Type{<:Formula{:wedepohl1973}}; compact::Bool = false)
    compact ? "Wedepohl" : "Wedepohl-Wilcox round-conductor surface impedances (1973)"
end

formulation_options(::FormulaMethod{<:Formula{:wedepohl1973}, typeof(internal_impedance)}) =
    FormulationOptions()

:wedepohl1973
