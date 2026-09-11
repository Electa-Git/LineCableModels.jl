using Measurements

(; frequencies = nothing,
    reference_options = (;),
    candidate_options = (;),
    variation = NoVariation()) ->
    _catalogue_uq_benchmark(
        :cable_525kv_land_no_armour_ac_flat,
        @__FILE__;
        frequencies,
        reference_options,
        candidate_options,
        variation,
    )
