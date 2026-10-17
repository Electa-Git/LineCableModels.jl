"""
Record the consumed FEM material laws and selected field assumptions.
"""
function formulation_record(formulation::LineCableModelsFEM)
    return merge((
        schema_version = 5,
        assumptions = (
            impedance = "Axial current-driven A_z/u_r finite-element equations",
            admittance = _quasi_full(formulation.options.data.physics) ?
                "Coupled first-order Maxwell A_z/A_t/phi equations; axial current supplies normalized leakage; vertical path voltage includes A_t/Gamma; Y = inv(P)" :
                "Scalar electrodynamic Helmholtz equation in surrounding media; equipotential terminals with unit transverse-current excitation; Y = inv(P)",
            earth = "Horizontal air and one semi-infinite soil; soil constitutive properties evaluated at each frequency",
            propagation = _quasi_full(formulation.options.data.physics) ?
                "Gamma -> 0 with A_t/Gamma and phi/Gamma retained; conduction and displacement; one coupled factorization" :
                "Gamma = 0; medium diffusion and displacement retained; independent Z/P blocks in one factorization",
            semicon_domain = "Passive material region, without electrical terminal ownership",
            enclosure = "Supported enclosures are represented by their material and terminal domains"
        )
    ),NamedTuple(formulation))
end
