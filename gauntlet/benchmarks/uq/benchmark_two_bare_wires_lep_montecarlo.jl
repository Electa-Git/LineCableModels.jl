using Measurements

(; frequencies = nothing,
    reference_options = (;),
    candidate_options = (;),
    variation = NoVariation()) ->
    _catalogue_uq_benchmark(
        :two_bare_wires,
        @__FILE__;
        frequencies,
        reference_options,
        candidate_options,
        variation,
    )
