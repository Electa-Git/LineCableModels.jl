using Measurements

(; frequencies = nothing,
    reference_options = (;),
    candidate_options = (;),
    variation = NoVariation()) ->
    _catalogue_uq_benchmark(
        :cable_30kv_na2xs2y_630mm2_trefoil,
        @__FILE__;
        frequencies,
        reference_options,
        candidate_options,
        variation,
    )
