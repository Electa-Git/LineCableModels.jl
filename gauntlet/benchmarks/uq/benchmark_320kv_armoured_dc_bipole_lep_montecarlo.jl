using Measurements

(; frequencies = nothing,
    reference_options = (;),
    candidate_options = (;),
    variation = NoVariation()) ->
    _catalogue_uq_benchmark(
        :cable_320kv_armoured_dc_bipole,
        @__FILE__;
        frequencies,
        reference_options,
        candidate_options,
        variation,
    )
