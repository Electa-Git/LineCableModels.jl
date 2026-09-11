using Gmsh

(; frequencies = nothing,
    reference_options = (;),
    candidate_options = (;),
    fem_options = (;),
    variation = NoVariation()) ->
    _catalogue_fem_benchmark(
        :cable_30kv_na2xs2y_630mm2_trefoil,
        @__FILE__;
        frequencies,
        reference_options,
        candidate_options,
        fem_options,
        variation,
    )
