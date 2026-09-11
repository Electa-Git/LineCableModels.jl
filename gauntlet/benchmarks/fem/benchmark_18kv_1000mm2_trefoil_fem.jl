using Gmsh

(; frequencies = nothing,
    reference_options = (;),
    candidate_options = (;),
    variation = NoVariation()) ->
    _catalogue_fem_benchmark(
        :cable_18kv_1000mm2_trefoil,
        @__FILE__;
        frequencies,
        reference_options,
        candidate_options,
        variation,
    )
