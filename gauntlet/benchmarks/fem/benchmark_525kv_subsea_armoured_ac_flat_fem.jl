using Gmsh

(; frequencies = nothing,
    reference_options = (;),
    candidate_options = (;),
    fem_options = (;),
    variation = NoVariation()) ->
    _catalogue_fem_benchmark(
        :cable_525kv_subsea_armoured_ac_flat,
        @__FILE__;
        frequencies,
        reference_options,
        candidate_options,
        fem_options,
        variation,
    )
