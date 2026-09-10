using Gmsh

# Fresh native FEM against the two registered analytical earth formulations.
(; frequencies = nothing,
    reference_options = (;),
    candidate_options = (;),
    variation = NoVariation()) -> begin
    model = load_case(:two_bare_wires; variation = ExactOverrides(; frequencies))
    for design in model.problem.system.designs
        regions = design.geometry.regions
        @assert length(regions) == 1 "The bare-wire benchmark must have no coating."
        @assert only(regions).source.material.kind === :conductor
    end
    physical = (reduce_bundle = false, kron_reduction = false,
        ideal_transposition = false)
    reference = Formulation(:LineCableModelsFEM; options = physical,
        fem_options = (mesh_policy = :remesh, keep_run_directory = true,
            gmsh_verbosity = 3, getdp_verbosity = 4,
            frequency_workers = 2, solver_threads = 1))
    candidates = Formulation(
        earth_impedance = Grid((:default, :Xue2018)),
        earth_admittance = Grid((:default, :Xue2018));
        combine = :zip, options = physical)
    @info "Two bare copper wires: FEM / proposed / Xue2018" frequencies
    benchmark_definition(model; id = :benchmark_two_bare_wires_fem,
        source_file = @__FILE__, collection = :fem,
        reference = BenchmarkCalculation(:fem, model.problem, reference;
            options = (trace = true, verbosity = (default = 1,))),
        formulations = candidates,
        report = BenchmarkTableDefinition(
            quantities = (:Z, :Y), bands = (:all, :dc, :harmonic, :narrow, :wide)))
end
