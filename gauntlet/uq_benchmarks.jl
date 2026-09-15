const UQ_MONTE_CARLO_TRIALS = 512

function uq_inner_formulation()
    return Formulation(
        earth_impedance = :default,
        earth_admittance = :default,
        shunt_model = :coaxial,
        insulation_admittance = formula(:default),
        options = (
            kron_reduction = false,
            reduce_bundle = false,
            ideal_transposition = false
        )
    )
end

# Measurement budget only; comparisons do not carry acceptance thresholds.
uq_timing_settings() = (performance=(samples=3, seconds=20.0),)
