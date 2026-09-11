using Gmsh

(; frequencies = nothing,
    reference_options = (;),
    candidate_options = (;),
    fem_options = (;),
    variation = NoVariation()) ->
    _catalogue_fem_benchmark(
        :cable_380kv_2000mm2_flatver,
        @__FILE__;
        frequencies,
        reference_options,
        candidate_options,
        fem_options,
        variation,
    )
