using Measurements

(; frequencies = nothing,
    reference_options = (;),
    candidate_options = (;),
    variation = NoVariation()) ->
    _catalogue_uq_benchmark(
        :cable_18kv_1000mm2_trefoil_homogenized,
        @__FILE__;
        frequencies,
        reference_options,
        candidate_options,
        variation,
    )
