using Gmsh

(; frequencies = nothing,
    reference_options = (;),
    candidate_options = (;),
    variation = NoVariation()) ->
    _catalogue_fem_benchmark(
        :solid_1000mm2_single,
        @__FILE__;
        frequencies,
        reference_options,
        candidate_options,
        variation,
    )
