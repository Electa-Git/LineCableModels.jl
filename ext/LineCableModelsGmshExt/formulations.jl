"""
Record the consumed FEM material laws and selected field assumptions.
"""
function formulation_record(formulation::LineCableModelsFEM)
    return merge((
        schema_version = 7,
        assumptions = (
            impedance = "Axial current-driven Maxwell equations; Z = -U/I + Gamma^2 P",
            admittance = "Coupled prescribed-Gamma Maxwell A_z/A_t/phi equations; axial current supplies normalized leakage; vertical path voltage includes A_t/Gamma; Y = inv(P)",
            earth = "Horizontal air and one semi-infinite soil; soil constitutive properties evaluated at each frequency",
            propagation = "Prescribed complex Gamma; exact A_t/Gamma and phi/Gamma variables with a regular zero limit; all Gamma^2 terms retained",
            semicon_domain = "Passive material region, without electrical terminal ownership",
            enclosure = "Supported enclosures are represented by their material and terminal domains"
        )
    ),NamedTuple(formulation))
end
