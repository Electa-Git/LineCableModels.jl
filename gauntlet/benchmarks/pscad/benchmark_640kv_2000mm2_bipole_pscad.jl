# Ordinary benchmark declaration. Construction performs no solve.
(; frequencies = nothing,
    reference_options = (;),
    candidate_options = (;),
    variation = NoVariation()) -> begin
    frequencies === nothing ||
        (variation=compose_variations(variation, ExactOverrides(; frequencies)))
    model=load_case(:cable_640kv_2000mm2_bipole; variation)
    reference_formulation=Formulation(
        :pscad; earth_impedance = :WedepohlWilcox1973
    )
    candidate_formulation=Formulation(
        earth_impedance = :Pollaczek1926,
        earth_admittance = :default,
        insulation_admittance = formula(:default),
        options = (
            kron_reduction = false,
            reduce_bundle = false,
            ideal_transposition = false
        )
    )
    return benchmark_definition(
        :benchmark_640kv_2000mm2_bipole_pscad, :cable_640kv_2000mm2_bipole, :pscad, @__FILE__, model,
        BenchmarkCalculation(:pscad, model.problem, reference_formulation; options = reference_options),
        BenchmarkCalculation(:lcm, model.problem, candidate_formulation; options = candidate_options),
        (; quantities = (:Z, :Y, :R, :L, :G, :C)), (;))
end
