using Gmsh

(; frequencies = nothing,
    reference_options = (;),
    candidate_options = (;),
    variation = NoVariation()) ->
    _catalogue_fem_benchmark(
        :cable_220kv_eaxecew_1x2500_252_trefoil,
        @__FILE__;
        frequencies,
        reference_options,
        candidate_options,
        variation,
    )
