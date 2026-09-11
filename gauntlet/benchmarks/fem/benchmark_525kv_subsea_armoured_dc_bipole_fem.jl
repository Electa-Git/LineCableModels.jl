using Gmsh

(; frequencies = nothing,
    reference_options = (;),
    candidate_options = (;),
    variation = NoVariation()) ->
    _catalogue_fem_benchmark(
        :cable_525kv_subsea_armoured_dc_bipole,
        @__FILE__;
        frequencies,
        reference_options,
        candidate_options,
        variation,
    )
