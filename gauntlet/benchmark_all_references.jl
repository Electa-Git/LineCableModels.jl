(; frequencies = nothing,
    pscad_reference_options = (;),
    fem_reference_options = (;),
    monte_carlo_reference_options = (;),
    candidate_options = (;),
    variation = NoVariation()) -> begin
    definitions = BenchmarkDefinition[]
    for case_id in PSCAD_CATALOGUE_CASE_IDS
        push!(definitions, benchmark_definition(
            _catalogue_benchmark_id(case_id, :pscad);
            frequencies,
            reference_options = pscad_reference_options,
            candidate_options,
            variation,
        ))
    end
    for case_id in CATALOGUE_CASE_IDS
        push!(definitions, benchmark_definition(
            _catalogue_benchmark_id(case_id, :fem);
            frequencies,
            reference_options = fem_reference_options,
            candidate_options,
            variation,
        ))
    end
    for case_id in CATALOGUE_CASE_IDS
        push!(definitions, benchmark_definition(
            _catalogue_benchmark_id(case_id, :uq);
            frequencies,
            reference_options = monte_carlo_reference_options,
            candidate_options,
            variation,
        ))
    end
    return definitions
end
