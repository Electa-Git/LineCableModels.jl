using Gmsh

(; frequencies = nothing,
    reference_options = (;),
    candidate_options = (;),
    fem_options = (;),
    variation = NoVariation()) ->
    _catalogue_fem_benchmark(
        :cable_525kv_1600mm2_bipole,
        @__FILE__;
        frequencies,
        reference_options,
        candidate_options,
        fem_options,
        variation,
    )
