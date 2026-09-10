# Ordinary benchmark declaration. Construction performs no solve.
(; frequencies = nothing,
    reference_options = (;),
    candidate_options = (;),
    variation = NoVariation()) -> begin
    frequencies === nothing ||
        (variation=compose_variations(variation, ExactOverrides(; frequencies)))
    model=load_case(:two_insulated_wires; variation)
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
        :benchmark_two_insulated_wires_pscad, :two_insulated_wires, :pscad, @__FILE__, model,
        BenchmarkCalculation(:pscad, model.problem, reference_formulation; options = reference_options),
        BenchmarkCalculation(:lcm, model.problem, candidate_formulation; options = candidate_options),
        (; quantities = (:Z, :Y, :R, :L, :G, :C)), (;))
end
