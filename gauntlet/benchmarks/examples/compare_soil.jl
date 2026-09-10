# CLI example: all numerical inputs are materialized before scheduling.
(;
    frequencies = 10.0 .^ range(-1, 7; length = 101)) -> begin
    model=load_case(:two_insulated_wires; variation = ExactOverrides(; frequencies))
    options=(reduce_bundle = false, kron_reduction = false, ideal_transposition = false)
    benchmark_definition(model; id=:compare_soil, collection=:manual,
        source_file=@__FILE__, reference=Formulation(;options),
        formulations=Formulation(earth_impedance=Grid((:default,:Pollaczek1926));options))
end
