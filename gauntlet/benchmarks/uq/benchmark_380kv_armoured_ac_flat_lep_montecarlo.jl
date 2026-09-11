using Measurements

(; frequencies = nothing,
    reference_options = (;),
    candidate_options = (;),
    variation = NoVariation()) ->
    _catalogue_uq_benchmark(
        :cable_380kv_armoured_ac_flat,
        @__FILE__;
        frequencies,
        reference_options,
        candidate_options,
        variation,
    )
