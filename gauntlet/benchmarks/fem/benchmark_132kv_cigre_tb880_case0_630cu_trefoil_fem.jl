using Gmsh

(; frequencies = nothing,
    reference_options = (;),
    candidate_options = (;),
    variation = NoVariation()) ->
    _catalogue_fem_benchmark(
        :cable_132kv_cigre_tb880_case0_630cu_trefoil,
        @__FILE__;
        frequencies,
        reference_options,
        candidate_options,
        variation,
    )
