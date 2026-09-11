@testitem "Gauntlet / case variations and correlation" tags=[:gauntlet_toolkit] setup=[
    GauntletSupport
] begin
    using Test
    using Measurements
    using LineCableModels
    using .GauntletSupport.Gauntlet

    function nested_groups(part)
        groups=Group[]
        if part isa Group
            push!(groups, part)
            append!(groups, nested_groups(part.item))
        elseif part isa Stack
            for item in part.items
                append!(groups, nested_groups(item))
            end
        elseif part isa Enclosure
            append!(groups, nested_groups(part.item))
            part.wall===nothing||append!(groups, nested_groups(part.wall))
        end
        return groups
    end

    first_load=load_case(:cable_18kv_1000mm2_trefoil)
    first_load.nominal_problem.system.connections[1][1]=99
    second_load=load_case(:cable_18kv_1000mm2_trefoil)
    @test second_load.nominal_problem.system.connections[1][1] == 1

    single_frequency=load_case(
        :cable_18kv_1000mm2_trefoil;
        variation = ExactOverrides(frequencies = [50.0])
    )
    @test single_frequency.expected_size == (9, 9, 1)
    @test single_frequency.problem.frequencies == [50.0]
    @test single_frequency.selected_parameters == [:frequencies]

    frequency_grid=load_case(
        :two_bare_wires;
        variation = ParameterGrids(frequencies = Grid(([50.0], [60.0])))
    )
    @test length(frequency_grid.problem) == 2
    @test [only(problem.frequencies) for problem in frequency_grid.problem] == [50.0, 60.0]

    @test_throws ArgumentError load_case(
        :two_bare_wires;
        variation = ExactOverrides(unknown_parameter = 1.0)
    )
    @test_throws ArgumentError load_case(
        :two_bare_wires;
        variation = RelativeStandardUncertainty(10.0; tags = (:not_a_tag,))
    )
    @test_throws ArgumentError load_case(
        :two_bare_wires;
        variation = ExactOverrides(frequencies = Grid(([50.0],)))
    )
    @test_throws ArgumentError load_case(:unknown_case)

    uncertain=load_case(
        :cable_18kv_1000mm2_trefoil;
        variation = RelativeStandardUncertainty(
            10.0; tags = (:geometry, :cable_layer)
        )
    )
    expected_selected=[
        :core_strand_diameter,
        :core_ring_lay_ratio_1,
        :core_ring_lay_ratio_2,
        :core_ring_lay_ratio_3,
        :core_ring_lay_ratio_4,
        :screen_wire_diameter,
        :screen_wire_lay_ratio,
        :semicon_tape_thickness,
        :inner_semicon_thickness,
        :insulation_thickness,
        :outer_semicon_thickness,
        :copper_tape_thickness,
        :copper_tape_width,
        :copper_tape_lay_ratio,
        :water_blocking_thickness,
        :aluminum_tape_thickness,
        :pe_face_thickness,
        :jacket_thickness
    ]
    @test uncertain.selected_parameters == expected_selected
    @test length(uncertain.problem) == 1
    for id in expected_selected
        descriptor=only(getproperty(uncertain.sources, id))
        nominal=getproperty(uncertain.definition.parameters, id).nominal
        @test LineCableModels.nominal(descriptor) == nominal
        @test LineCableModels.uncertainty(descriptor) ≈ 0.1abs(nominal)
    end
    @test uncertain.sources.formation_clearance_ratio == 0.2
    @test uncertain.sources.screen_wires == 49
    @test uncertain.sources.earth_rho == 100.0

    nominal_system=uncertain.nominal_problem.system
    nominal_radius=outer_radius(first(nominal_system.designs))
    nominal_spacing=hypot(
        nominal_system.positions[1].x-nominal_system.positions[2].x,
        nominal_system.positions[1].y-nominal_system.positions[2].y
    )
    @test nominal_radius ≈ 34.575e-3
    @test nominal_spacing ≈ 2.2nominal_radius
    @test nominal_spacing - 2nominal_radius ≈ 0.2nominal_radius

    expanded_armor=load_case(
        :cable_525kv_1600mm2_bipole;
        variation = ExactOverrides(armor_wire_diameter = 1.2*5.827e-3)
    )
    expanded_design=first(expanded_armor.problem.system.designs)
    armor_group=only(filter(
        item->item isa Group&&item.name===:armor,
        nested_groups(expanded_design.root)
    ))
    armor=only(filter(
        value->value.name===:armor,
        LineCableModels.DataModel.flatten(
            expanded_design, 50.0
        )
    )).conductor
    armor_wire_radius=armor_group.item.primitive.r
    @test armor.num_wires == 68
    @test armor.r_in >=
          armor_wire_radius / sinpi(1 / armor.num_wires) - armor_wire_radius
    expanded_bedding=only(filter(
        source->source.source.tag===:sheath_bedding,
        expanded_design.geometry.regions
    ))
    @test thickness(expanded_bedding.primitive) > 3.0e-3

    nominal_armor=load_case(:cable_525kv_1600mm2_bipole)
    nominal_design=first(nominal_armor.problem.system.designs)
    nominal_bedding=only(filter(
        source->source.source.tag===:sheath_bedding,
        nominal_design.geometry.regions
    ))
    nominal_unbuffered_outer_radius=63.22185e-3+5.827e-3+10.0e-3
    @test thickness(nominal_bedding.primitive) ≈
          3.0e-3 + 0.2nominal_unbuffered_outer_radius

    materialized=only(uncertain.problem)
    uncertain_system=materialized.system
    uncertain_radius=outer_radius(first(uncertain_system.designs))
    uncertain_spacing=hypot(
        uncertain_system.positions[1].x-uncertain_system.positions[2].x,
        uncertain_system.positions[1].y-uncertain_system.positions[2].y
    )
    @test uncertain_spacing - 2uncertain_radius ≈ 0.2uncertain_radius

    geometry=first(materialized.system.designs).geometry.regions
    first_tape=thickness(only(filter(
        source->source.source.tag===:core_semicon_tape_inner,
        geometry
    )).primitive)
    second_tape=thickness(only(filter(
        source->source.source.tag===:core_semicon_tape_outer,
        geometry
    )).primitive)
    @test Measurements.cov(first_tape, second_tape) > 0
    @test Measurements.cov(first_tape, second_tape) ≈
          Measurements.uncertainty(first_tape)^2

    uq_model=load_case(
        :cable_132kv_630mm2_flathor;
        variation = RelativeStandardUncertainty(
            10.0; tags = (:geometry, :cable_layer)
        )
    )
    uq_selected=[
        :core_strand_diameter,
        :core_lay_ratio,
        :screen_wire_diameter,
        :screen_lay_ratio,
        :semicon_tape_thickness,
        :inner_semicon_thickness,
        :insulation_thickness,
        :outer_semicon_thickness,
        :copper_tape_thickness,
        :copper_tape_width,
        :copper_tape_lay_ratio,
        :water_blocking_thickness,
        :aluminum_tape_thickness,
        :jacket_thickness
    ]
    @test uq_model.selected_parameters == uq_selected
    @test length(uq_model.problem) == 1
    for parameter in values(uq_model.definition.parameters)
        source=getproperty(uq_model.sources, parameter.id)
        if parameter.id in uq_selected
            descriptor=only(source)
            @test LineCableModels.nominal(descriptor) == parameter.nominal
            @test LineCableModels.uncertainty(descriptor) ≈
                  0.1abs(parameter.nominal)
        else
            @test isequal(source, parameter.nominal)
        end
    end
end


@testitem "Gauntlet / UQ moment comparison contract" tags=[:gauntlet_toolkit] setup=[
    GauntletSupport
] begin
    using Test
    using LineCableModels
    using LineCableModels.Engine
    using .GauntletSupport.Gauntlet

    @test parentmodule(MomentResult) === GauntletSupport.Gauntlet
    @test MomentResult === GauntletSupport.MomentResult

    frequencies_value=[1.0, 10.0]
    ports=["a", "b"]
    values=map((R = 1.0, L = 2.0, C = 3.0, G = 4.0)) do scale
        (
            mean = fill(scale, 2, 2, 2),
            std = fill(scale/10, 2, 2, 2)
        )
    end
    reference=MomentResult(values, frequencies_value, :pul, PhaseDomain, ports)
    equal_comparison=compare(reference, reference)
    tolerance=(
        mean = map(_->(absolute = 0.0, relative = 0.0), values),
        std = map(_->(absolute = 0.0, relative = 0.0), values)
    )
    @test reference ≈ reference
    @test moment_comparison_passes(equal_comparison, tolerance)
    @test all(iszero, equal_comparison.errors.R.mean.absolute)

    changed_values=merge(values, (
        R = (mean = copy(values.R.mean), std = copy(values.R.std)),
    ))
    changed_values.R.mean[1, 2, :].=1.5
    changed=MomentResult(changed_values, frequencies_value, :pul, PhaseDomain, ports)
    changed_comparison=compare(reference, changed)
    @test changed_comparison.errors.R.mean.absolute[1, 2] == 0.5
    @test !moment_comparison_passes(changed_comparison, tolerance)

    zero_values=map(values) do product
        (mean = zeros(size(product.mean)), std = zeros(size(product.std)))
    end
    small_values=map(zero_values) do product
        (mean = fill(1.0e-12, size(product.mean)), std = copy(product.std))
    end
    zero_reference=MomentResult(zero_values, frequencies_value, :pul, PhaseDomain, ports)
    small_candidate=MomentResult(small_values, frequencies_value, :pul, PhaseDomain, ports)
    floor_tolerance=(
        mean = map(_->(absolute = 1.0e-11, relative = 0.0), values),
        std = map(_->(absolute = 0.0, relative = 0.0), values)
    )
    @test all(ismissing, compare(zero_reference, small_candidate).errors.R.mean.relative)
    @test moment_comparison_passes(
        compare(zero_reference, small_candidate), floor_tolerance
    )

    @test_throws ArgumentError compare(
        reference,
        MomentResult(values, [1.0, 11.0], :pul, PhaseDomain, ports)
    )
    @test_throws ArgumentError compare(
        reference,
        MomentResult(values, frequencies_value, :total, PhaseDomain, ports)
    )
    @test_throws ArgumentError compare(
        reference,
        MomentResult(values, frequencies_value, :pul, PhaseDomain, reverse(ports))
    )
    wrong_shape=merge(values, (
        R = (mean = zeros(1, 1, 2), std = zeros(1, 1, 2)),
    ))
    @test_throws DimensionMismatch compare(
        reference,
        MomentResult(wrong_shape, frequencies_value, :pul, PhaseDomain, ports)
    )
end


@testitem "Gauntlet / fixed-seed Monte Carlo reproducibility" tags=[:gauntlet_toolkit] setup=[
    GauntletSupport
] begin
    using Test
    using LineCableModels
    using LineCableModels.Engine
    using .GauntletSupport.Gauntlet

    model=load_case(
        :cable_132kv_630mm2_flathor;
        variation = compose_variations(
            ExactOverrides(frequencies = [50.0]),
            RelativeStandardUncertainty(
                10.0; tags = (:geometry, :cable_layer)
            )
        )
    )
    inner=Formulation(
        earth_impedance = :Pollaczek1926,
        earth_admittance = :default,
        insulation_admittance = formula(:default),
        options = (
            kron_reduction = false,
            reduce_bundle = false,
            ideal_transposition = false
        )
    )
    formulation=MonteCarlo(
        inner;
        trials = 8,
        seed = 0x51eed,
        distribution = :normal,
        return_samples = false,
        return_histograms = false
    )
    problem=ParametricProblem(model.problem)
    first_result=compute(problem, formulation)
    second_result=compute(problem, formulation)
    @test first_result.root_seed == second_result.root_seed == UInt64(0x51eed)
    @test first_result.trial_counts == second_result.trial_counts == [8]
    @test statistics(first_result) == statistics(second_result)
    @test samples(first_result) === nothing
    @test histograms(first_result) === nothing
end


@testitem "Gauntlet / formula comparison dispatch" tags=[:gauntlet_toolkit] setup=[
    GauntletSupport
] begin
    using Test
    using LineCableModels
    using LineCableModels.Engine
    using .GauntletSupport.Gauntlet

    previous_mode=get(ENV, "LINECABLEMODELS_GAUNTLET_MODE", nothing)
    try
        ENV["LINECABLEMODELS_GAUNTLET_MODE"]="live"
        model=load_case(
            :two_insulated_wires;
            variation = ExactOverrides(frequencies = [50.0])
        )
        pollaczek=Formulation(
            earth_impedance = :Pollaczek1926,
            earth_admittance = :default,
            insulation_admittance = formula(:default),
            options = (kron_reduction = false, reduce_bundle = false)
        )
        saad=Formulation(
            earth_impedance = :Saad1996,
            earth_admittance = :default,
            insulation_admittance = formula(:default),
            options = (kron_reduction = false, reduce_bundle = false)
        )
        benchmark=benchmark_definition(
            :benchmark_owned_formula_fixture,
            model.id,
            :owned,
            @__FILE__,
            model,
            BenchmarkCalculation(
                :pollaczek, model.problem, pollaczek
            ),
            BenchmarkCalculation(
                :saad, model.problem, saad
            ),
            (;),
            (;)
        )
        outcome=run_benchmark(benchmark)
        @test outcome.reference isa LineParameters
        @test outcome.candidate isa LineParameters
        @test outcome.passes === nothing
        @test outcome.metadata.calculations.reference.id === :pollaczek
        @test outcome.metadata.calculations.candidate.id === :saad
        @test size(Z(outcome.reference)) == (2, 2, 1)
    finally
        previous_mode===nothing ?
        delete!(ENV, "LINECABLEMODELS_GAUNTLET_MODE") :
        (ENV["LINECABLEMODELS_GAUNTLET_MODE"]=previous_mode)
    end
end


@testitem "Gauntlet / local performance comparison" tags=[:gauntlet_toolkit] setup=[
    GauntletSupport,
    TestFixtures
] begin
    using Test
    using LineCableModels
    using .GauntletSupport.Gauntlet
    using .TestFixtures

    problem=TestFixtures.line_parameters_problem()
    timing=benchmark_local(
        (; problem, formulation = LineCableModels.Formulation()); samples = 1, seconds = 1
    )
    @test timing.samples == 1
    @test timing.minimum_seconds >= 0
    @test timing.median_seconds >= timing.minimum_seconds
    @test timing.bytes >= 0
    @test timing.allocations >= 0
    @test timing.environment.julia_version == string(VERSION)
    @test timing.environment.cpu == Sys.CPU_NAME
    @test timing.environment.cpu_threads == Sys.CPU_THREADS
    @test timing.environment.blas_threads == GauntletSupport.Gauntlet.BLAS.get_num_threads()

    tolerance=(median_time_ratio = 1.2, bytes_ratio = 1.05, allocations_ratio = 1.05)
    diagnostic=performance_comparison(timing, timing, tolerance)
    @test diagnostic.comparable == !gauntlet_instrumented()
    @test diagnostic.passes === (gauntlet_instrumented() ? nothing : true)

    # Exercise comparison arithmetic with declared inputs, independently of
    # this test process's timing/coverage instrumentation. The actual measured
    # record above remains subject to the declared comparison checks.
    accepted=(;
        timing...,
        median_seconds = 1.0,
        bytes = 1000,
        allocations = 100
    )
    current=(; accepted..., median_seconds = 1.1, bytes = 1010, allocations = 101)
    compared=performance_comparison(accepted, current, tolerance; instrumented = false)
    @test compared.comparable
    @test compared.passes
    @test compared.ratios == (median_time = 1.1, bytes = 1.01, allocations = 1.01)
    for (field, value) in ((:median_seconds, 1.3), (:bytes, 1060), (:allocations, 106))
        slower=merge(current, NamedTuple{(field,)}((value,)))
        @test !performance_comparison(accepted, slower, tolerance; instrumented = false).passes
    end
    no_allocations=(; accepted..., bytes = 0, allocations = 0)
    equal=performance_comparison(no_allocations, no_allocations, tolerance; instrumented = false)
    @test equal.passes
    @test equal.ratios.bytes == equal.ratios.allocations == 1.0
    allocated=performance_comparison(no_allocations, current, tolerance; instrumented = false)
    @test !allocated.passes
    @test isinf(allocated.ratios.bytes) && isinf(allocated.ratios.allocations)
    other_environment=(;
        accepted...,
        environment = (; accepted.environment..., julia_version = "different")
    )
    diagnostic=performance_comparison(other_environment, current, tolerance; instrumented = false)
    @test !diagnostic.comparable
    @test diagnostic.passes === nothing
    for (field,
        value) in (
        (:cpu, "different CPU"), (:blas_threads, timing.environment.blas_threads+1))
        mismatched=merge(
            accepted, (;
                environment = merge(accepted.environment, NamedTuple{(field,)}((value,)))))
        comparison=performance_comparison(mismatched, current, tolerance; instrumented = false)
        @test !comparison.comparable
        @test comparison.passes === nothing
    end

    instrumented=performance_comparison(
        accepted,
        current,
        tolerance;
        instrumented = true
    )
    @test !instrumented.comparable
    @test instrumented.passes === nothing
end
