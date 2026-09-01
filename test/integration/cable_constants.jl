@testitem "Engine / cable constants / earth-free assemblies and physical invariants" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestFixtures,
    TestNumerics
] begin
    design=TestFixtures.mv_cable_design()
    constants=CableConstants(design)
    row=only(constants)

    @test row.core === :core
    @test row.R > 0
    @test row.L > 0
    @test row.C > 0
    @test row.G >= 0
    @test constants.frequency == 50.0
    @test basis(constants) === :pul
    @test resistance(constants) === constants.R
    @test inductance(constants) === constants.L
    @test capacitance(constants) === constants.C
    @test conductance(constants) === constants.G

    packed_material=Material(:conductor, 1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
    packed_dielectric=Material(:insulator, 1e14, 2.3)
    packed_design=build(
        CableDesign,
        "hexagonally-packed-core",
        terminal(
            :packed_core,
            stranded(
                packed_material;
                shape = Disk(0.5e-3),
                layers = 2,
                lay = (LayRatio(15), LayRatio(11))
            )
        ),
        insulation(packed_dielectric; t = 1e-3)
    )
    packed_constants=CableConstants(packed_design; frequency = 50.0)
    @test all(isfinite,
        Iterators.flatten(
            (packed_constants.R, packed_constants.L,
            packed_constants.C, packed_constants.G)
        ))

    problem=CableConstantsProblem(design)
    formulation=CableConstantsFormulation()
    @test problem.design === design
    blueprint=@inferred LineCableModels.Engine.flatten(
        LineCableModelsCoaxial(), design
    )
    @test blueprint isa LineCableModels.Engine.CableBlueprint{Float64}
    @test blueprint.assembly_ranges == [1:3]
    @test getproperty.(blueprint.conductors, :terminal) == design.terminal_order
    @test fieldnames(typeof(blueprint)) == (
        :cable_id, :conductors, :dielectrics, :dielectric_ranges,
        :assembly_ranges
    )
    @test :frequency ∉ fieldnames(typeof(blueprint))
    @test compute(problem, formulation) == constants
    @test @inferred(compute(problem, formulation)) == constants
    @test_throws ArgumentError compute(problem, formulation; options = (trace = true,))
    @test_throws MethodError CableConstants(design; position = at(x = 0, y = -2))
    @test_throws MethodError CableConstants(design; earth_props = Earth(rho = 200))

    homogeneous=LineCableModels.homogenize(design)
    source_components=LineCableModels.DataModel.flatten(design, 50.0)
    reduced_components=LineCableModels.DataModel.flatten(
        homogeneous,
        50.0
    )
    @test getproperty.(reduced_components, :name) ==
          getproperty.(source_components, :name)
    for (source, reduced) in zip(source_components, reduced_components)
        @test reduced.conductor.resistance ≈ source.conductor.resistance
        @test reduced.conductor.gmr ≈ source.conductor.gmr
        @test reduced.dielectric.shunt_capacitance ≈
              source.dielectric.shunt_capacitance
        @test reduced.dielectric.shunt_conductance ≈
              source.dielectric.shunt_conductance
        @test reduced.dielectric.material.mu_r ≈ source.dielectric.material.mu_r
    end
    flattening_formulation=CableConstantsFormulation(
        insulation_admittance = formula(:Ametani2004)
    )
    source_at_flattening_frequency=CableConstants(
        design;
        formulation = flattening_formulation,
        frequency = 50.0
    )
    homogeneous_constants=CableConstants(
        homogeneous;
        formulation = flattening_formulation,
        frequency = 50.0
    )
    @test homogeneous_constants.R ≈ source_at_flattening_frequency.R
    @test homogeneous_constants.L ≈ source_at_flattening_frequency.L
    @test homogeneous_constants.C ≈ source_at_flattening_frequency.C
    @test homogeneous_constants.G ≈ source_at_flattening_frequency.G

    sixty_hertz=CableConstants(design; frequency = 60.0)
    @test all(isfinite,
        Iterators.flatten(
            (sixty_hertz.R, sixty_hertz.L, sixty_hertz.C, sixty_hertz.G)
        ))
    @test sixty_hertz.frequency == 60.0
    @test_throws DomainError CableConstantsProblem(design; frequency = 1e6)

    hot=CableConstants(design; temperature = 80.0)
    uncorrected=CableConstants(
        design;
        temperature = 80.0,
        formulation = CableConstantsFormulation(
            options = (temperature_correction = false,)
        )
    )
    @test hot.R[1] > constants.R[1]
    @test uncorrected.R ≈ constants.R

    lossless_dielectric=CableConstants(
        design;
        frequency = 50.0,
        formulation = CableConstantsFormulation(
            insulation_admittance = formula(:Gustavsen2013)
        )
    )
    lossy_dielectric=CableConstants(design; frequency = 50.0)
    @test lossless_dielectric.G[1] < lossy_dielectric.G[1]
    @test lossless_dielectric.C[1] ≈ lossy_dielectric.C[1]

    traced_system=build(
        LineCableSystem,
        design,
        Pose2(0.0, -1.0);
        connections = Dict(:core=>1, :sheath=>0, :jacket=>0)
    )
    traced_problem=LineParametersProblem(
        traced_system;
        earth_props = Earth(rho = 100.0),
        frequencies = [50.0, 60.0]
    )
    traced=compute(
        traced_problem,
        Formulation(options = (ideal_transposition = false,));
        options = (trace = true,)
    )
    primitive=details(traced).trace
    for (frequency_index, (frequency, expected)) in enumerate((
        (50.0, constants),
        (60.0, sixty_hertz)
    ))
        local_Z=kronify(
            Matrix(primitive.Zin[:, :, frequency_index]),
            [1, 0, 0]
        )[1, 1]
        local_Y=im*2π*frequency*inv(
            Matrix(primitive.Pin[:, :, frequency_index])
        )
        @test only(expected).R ≈ real(local_Z)
        @test only(expected).L ≈ imag(local_Z)/(2π*frequency)
        @test only(expected).G ≈ real(local_Y[1, 1])
        @test only(expected).C ≈ imag(local_Y[1, 1])/(2π*frequency)
    end

    copper=Material(:conductor, 1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
    dielectric=Material(:insulator, 1e14, 2.3)
    bare_design=build(
        CableDesign,
        "bare-core",
        Group(:bare_core, Region(:bare_metal, Disk(0.01), copper))
    )
    @test_throws ArgumentError CableConstants(bare_design)

    function unsheathed_member(core_name)
        Stack(
            Group(core_name,
                Region(Symbol(core_name, :_metal), Disk(0.01), copper)),
            Region(Symbol(core_name, :_insulation),
                Annulus(0.01, 0.02), dielectric)
        )
    end
    unsheathed_design=build(
        CableDesign,
        "unsheathed-core",
        unsheathed_member(:unsheathed_core)
    )
    unsheathed_problem=CableConstantsProblem(unsheathed_design; frequency = 50.0)
    @test unsheathed_problem.design === unsheathed_design
    unsheathed=only(compute(unsheathed_problem, CableConstantsFormulation()))
    @test unsheathed.core === :unsheathed_core
    @test unsheathed.R > 0
    @test unsheathed.L > 0
    @test unsheathed.C > 0
    @test unsheathed.G >= 0

    unsheathed_pair=build(
        CableDesign,
        "two-unsheathed-assemblies",
        assembly(
            at(unsheathed_member(:left_unsheathed_core), -0.03, 0.0),
            at(unsheathed_member(:right_unsheathed_core), 0.03, 0.0)
        )
    )
    unsheathed_pair_problem=CableConstantsProblem(
        unsheathed_pair;
        frequency = 50.0
    )
    @test unsheathed_pair_problem.design === unsheathed_pair
    unsheathed_pair_constants=compute(
        unsheathed_pair_problem,
        CableConstantsFormulation()
    )
    @test unsheathed_pair_constants.cores ==
          [:left_unsheathed_core, :right_unsheathed_core]
    @test all(>(0), unsheathed_pair_constants.R)
    @test all(>(0), unsheathed_pair_constants.L)
    @test all(>(0), unsheathed_pair_constants.C)
    @test all(>=(0), unsheathed_pair_constants.G)

    eccentric_member=Stack(
        Group(:eccentric_core, Region(:eccentric_metal, Disk(0.01), copper)),
        Region(:eccentric_insulation, Annulus(0.01, 0.02), dielectric)
    )
    common_return=Group(
        :common_return,
        Region(:common_return_metal, Annulus(0.05, 0.051), copper)
    )
    eccentric_design=build(
        CableDesign,
        "eccentric-return",
        assembly(at(eccentric_member, -0.02, 0.0), at(common_return, 0.0, 0.0))
    )
    @test_throws ArgumentError CableConstants(eccentric_design; frequency = 50.0)

    function coax_member(core_name, sheath_name)
        Stack(
            Group(core_name, Region(Symbol(core_name, :_metal), Disk(0.01), copper)),
            Region(Symbol(core_name, :_insulation), Annulus(0.01, 0.02), dielectric),
            Group(sheath_name,
                Region(Symbol(sheath_name, :_metal), Annulus(0.02, 0.021), copper))
        )
    end
    left=coax_member(:left_core, :left_sheath)
    right=coax_member(:right_core, :right_sheath)
    multicore=build(CableDesign, "two-assembly",
        assembly(at(left, -0.03, 0.0), at(right, 0.03, 0.0)))
    combined=CableConstants(multicore; frequency = 50.0)
    left_constants=CableConstants(build(CableDesign, "left", left); frequency = 50.0)
    right_constants=CableConstants(build(CableDesign, "right", right); frequency = 50.0)
    @test combined.cores == [:left_core, :right_core]
    @test combined.R ≈ [only(left_constants).R, only(right_constants).R]
    @test combined.L ≈ [only(left_constants).L, only(right_constants).L]
    @test combined.C ≈ [only(left_constants).C, only(right_constants).C]
    @test combined.G ≈ [only(left_constants).G, only(right_constants).G]

    multicore_system=build(
        LineCableSystem,
        multicore,
        Pose2(0.0, -1.0);
        connections = Dict(
            :left_core=>1,
            :left_sheath=>0,
            :right_core=>2,
            :right_sheath=>0
        )
    )
    multicore_parameters=compute(
        LineParametersProblem(
            multicore_system;
            earth_props = Earth(rho = 100.0),
            frequencies = [50.0]
        );
        options = (trace = true,)
    )
    @test size(multicore_parameters.Z) == (2, 2, 1)
    @test size(multicore_parameters.Y) == (2, 2, 1)
    @test all(isfinite, multicore_parameters.Z)
    @test all(isfinite, multicore_parameters.Y)
    @test details(multicore_parameters).trace.cable_map == [1, 1, 2, 2]

    function many_conductor_design()
        parts=AbstractCablePart[
            Group(:terminal_1, Region(:metal_1, Disk(0.005), copper))
        ]
        outer_radius=0.005
        for index in 2:6
            dielectric_outer=outer_radius+0.0005
            push!(parts, Region(
                Symbol(:dielectric_, index-1),
                Annulus(outer_radius, dielectric_outer),
                dielectric
            ))
            conductor_outer=dielectric_outer+0.0003
            push!(parts, Group(
                Symbol(:terminal_, index),
                Region(
                    Symbol(:metal_, index),
                    Annulus(dielectric_outer, conductor_outer),
                    copper
                )
            ))
            outer_radius=conductor_outer
        end
        return build(CableDesign, "six-conductor-coax", Stack(parts))
    end
    many_layer_design=many_conductor_design()
    many_layer_constants=@inferred CableConstants(
        many_layer_design;
        frequency = 50.0
    )
    @test many_layer_constants.cores == [:terminal_1]
    @test all(isfinite, Iterators.flatten((
        many_layer_constants.R,
        many_layer_constants.L,
        many_layer_constants.C,
        many_layer_constants.G
    )))

    designs=build(
        CableDesign,
        Grid(("constants-a", "constants-b")),
        design.root
    )
    constant_space=CableConstants(designs)
    @test constant_space isa Gridspace{CableConstants}
    @test length(constant_space) == 2
    @test all(value -> value isa CableConstants, constant_space)
    @test first(constant_space) == CableConstants(first(designs))

    problem_space=CableConstantsProblem(designs; frequency = 50.0)
    @test problem_space isa Gridspace{CableConstantsProblem}
    selected_problem=first(LineCableModels.ParametricBuilder.points(problem_space))
    @test selected_problem isa LineCableModels.Gridpoint{CableConstantsProblem}
    calculated=compute(
        ParametricProblem(problem_space),
        Combinatorial(CableConstantsFormulation())
    )
    @test length(calculated) == 2
    @test all(value -> value isa CableConstants, calculated)
    @test first(calculated) == CableConstants(first(designs); frequency = 50.0)

    sampled=compute(
        ParametricProblem(CableConstantsProblem(
            build(CableDesign, Grid(("sampled",)), design.root);
            frequency = 50.0
        )),
        MonteCarlo(
            CableConstantsFormulation();
            trials = 2,
            seed = 11,
            return_samples = true,
            return_histograms = true
        )
    )
    @test sampled isa MonteCarloResult{<:CableConstants}
    @test size(only(samples(sampled)).R) == (1, 2)
    @test keys(only(statistics(sampled))) == (:R, :L, :C, :G)
end

@testitem "Fixtures / factories / mutable state is never shared" tags=[:integration] setup=[
    TestFixtures,
] begin
    first_design=TestFixtures.mv_cable_design()
    second_design=TestFixtures.mv_cable_design()
    @test first_design !== second_design
    @test first_design.root !== second_design.root
    @test first_design.geometry.regions !== second_design.geometry.regions

    first_system=TestFixtures.three_phase_system()
    second_system=TestFixtures.three_phase_system()
    @test first_system !== second_system
    @test first_system.designs !== second_system.designs
    @test first_system.positions !== second_system.positions
    @test ncables(first_system) == ncables(second_system) == 3
    @test nphases(first_system) == nphases(second_system) == 3

    first_parameters=TestFixtures.two_conductor_results()
    second_parameters=TestFixtures.two_conductor_results()
    @test first_parameters !== second_parameters
    @test first_parameters.Z.values !== second_parameters.Z.values
    first_parameters.Z.values[1]=complex(Inf, Inf)
    @test all(isfinite, second_parameters.Z.values)

    first_monte_carlo=TestFixtures.cable_monte_carlo_result()
    second_monte_carlo=TestFixtures.cable_monte_carlo_result()
    @test first_monte_carlo !== second_monte_carlo
    first_samples=only(samples(first_monte_carlo))
    second_samples=only(samples(second_monte_carlo))
    @test first_samples.R !== second_samples.R
    first_samples.R[1, 1]=Inf
    @test all(isfinite, second_samples.R)
end
