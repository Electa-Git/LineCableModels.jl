const PSCAD_FORMULATIONS = Dict(
    :overhead_deri_carson_lossless => (
        pscad_parameter = :EarthForm2,
        pscad_value = 0,
        earth_impedance = "Carson",
        earth_admittance = "Electrostatic images",
        insulation_admittance = "Lossless insulation (ideal dielectric)"
    ),
    :underground_wedepohl_pollaczek_lossless => (
        pscad_parameter = :EarthForm,
        pscad_value = 0,
        earth_impedance = "Pollaczek",
        earth_admittance = "Pollaczek",
        insulation_admittance = "Lossless insulation (ideal dielectric)"
    )
)

function formulation_spec(name::Symbol)
    haskey(PSCAD_FORMULATIONS, name) || throw(ArgumentError(
        "unsupported PSCAD gauntlet formulation $(repr(name))",
    ))
    return PSCAD_FORMULATIONS[name]
end

function validate_formulation(case::GauntletCase)
    expected = formulation_spec(case.pscad_formulation)
    actual = case.formulation
    description(actual.earth_impedance) == expected.earth_impedance ||
        throw(ArgumentError("case earth-impedance formulation contradicts its PSCAD mapping"))
    description(actual.earth_admittance) == expected.earth_admittance ||
        throw(ArgumentError("case earth-admittance formulation contradicts its PSCAD mapping"))
    description(actual.insulation_admittance) == expected.insulation_admittance ||
        throw(ArgumentError(
            "case insulation-admittance formulation contradicts its PSCAD mapping",
        ))
    return expected
end
