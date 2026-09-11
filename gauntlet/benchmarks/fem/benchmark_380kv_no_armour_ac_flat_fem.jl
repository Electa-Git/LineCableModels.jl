using Gmsh

(; frequencies = nothing,
    reference_options = (;),
    candidate_options = (;),
    variation = NoVariation()) ->
    _catalogue_fem_benchmark(
        :cable_380kv_no_armour_ac_flat,
        @__FILE__;
        frequencies,
        reference_options,
        candidate_options,
        variation,
    )
