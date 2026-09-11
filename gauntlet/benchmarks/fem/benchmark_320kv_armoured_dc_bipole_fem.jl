using Gmsh

(; frequencies = nothing,
    reference_options = (;),
    candidate_options = (;),
    fem_options = (;),
    variation = NoVariation()) ->
    _catalogue_fem_benchmark(
        :cable_320kv_armoured_dc_bipole,
        @__FILE__;
        frequencies,
        reference_options,
        candidate_options,
        fem_options,
        variation,
    )
