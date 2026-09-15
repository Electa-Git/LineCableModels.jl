"""
$(TYPEDSIGNATURES)

**Identification.** Longitudinal magnetic impedance of one concentric
insulation region in Ametani's single-core cable formulation.

**Reference.** A. Ametani, “A General Formulation of Impedance and Admittance
of Cables,” *IEEE Transactions on Power Apparatus and Systems*, PAS-99(3),
902–910, 1980.
"""
function description(::Type{<:Formula{:ametani1980}}; compact::Bool=false)
    compact ? "Ametani" : "Ametani coaxial-insulation magnetic impedance (1980)"
end

"Delegate the normalized Ametani route to the package default implementation."
insulation_impedance(
        ::Val{:ametani1980}, args...
) = insulation_impedance(Val(:default), args...)

computation_options(::FormulaMethod{:ametani1980, typeof(insulation_impedance)}) = (;)

:ametani1980
