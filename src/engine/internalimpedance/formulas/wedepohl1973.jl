"""
$(TYPEDSIGNATURES)

Identify the Wedepohl-Wilcox analytical approximations for the inner, outer,
and transfer surface impedances of a round conductor, in Ω/m. A solid cylinder
requires only the outer coefficient; an annulus uses all three coefficients.

This scientific identity is registered for backend dispatch. Evaluation by the
owned coaxial backend is not yet implemented and has no numerical fallback.
PSCAD's line-constants program uses these approximations for conductor surfaces.

Reference: L. M. Wedepohl and D. J. Wilcox, “Transient Analysis of Underground
Power-Transmission Systems: System-Model and Wave-Propagation Characteristics,”
*Proceedings of the IEE*, 120, 253–260, 1973. DOI: 10.1049/piee.1973.0056.
PSCAD 5.1 help reproduces the surface approximations in *Deriving System Y and Z
Matrices*, Eqs. (8-8), (8-9), and (8-11)–(8-13).
"""
function description(::Type{<:Formula{:wedepohl1973}}; compact::Bool = false)
    compact ? "Wedepohl" : "Wedepohl-Wilcox round-conductor surface impedances (1973)"
end

function (formula::Formula{:wedepohl1973})(
        r_in::T, r_ex::T, rho_c::T, mur_c::T, jω::Complex{T}) where {T <: Real}
    throw(ArgumentError("internal_impedance :wedepohl1973: not yet implemented for the coaxial backend"))
end

function internal_impedance(::Formula{:wedepohl1973},
        kind::Union{Val{:inner}, Val{:outer}, Val{:transfer}}, functor, workspace)
    throw(ArgumentError("internal_impedance :wedepohl1973 ($kind): not yet implemented for the coaxial backend"))
end

formulation_options(::FormulaMethod{<:Formula{:wedepohl1973}, typeof(internal_impedance)}) =
    FormulationOptions()

:wedepohl1973
