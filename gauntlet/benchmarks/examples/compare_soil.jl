# CLI example: all numerical inputs are materialized before scheduling.
(;
    frequencies = 10.0 .^ range(-1, 7; length = 101)) -> begin
    model=load_case(:two_insulated_wires; variation = ExactOverrides(; frequencies))
    options=(reduce_bundle = false, kron_reduction = false, ideal_transposition = false)
    reference=BenchmarkCalculation(:default, model.problem, Formulation(; options))
    candidate=BenchmarkCalculation(
        :pollaczek, model.problem, Formulation(earth_impedance = :Pollaczek1926; options))
    benchmark_definition(
        :compare_soil, model.id, :manual, @__FILE__, model, reference, candidate,
        (; quantities = (:Z, :Y, :R, :L, :G, :C)), (;))
end
