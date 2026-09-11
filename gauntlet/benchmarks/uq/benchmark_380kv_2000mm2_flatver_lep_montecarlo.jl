using Measurements

(; frequencies = nothing,
    reference_options = (;),
    candidate_options = (;),
    variation = NoVariation()) ->
    _catalogue_uq_benchmark(
        :cable_380kv_2000mm2_flatver,
        @__FILE__;
        frequencies,
        reference_options,
        candidate_options,
        variation,
    )
