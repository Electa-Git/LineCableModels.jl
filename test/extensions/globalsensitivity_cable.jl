@testitem "GlobalSensitivity / Sobol / real cable-domain smoke" tags=[:extension] setup=[
    EngineTestSupport, UseEngineSupport] begin
    using GlobalSensitivity
    import LineCableModels.ParametricBuilder as PB

    copper = PB.Material(
        kind = :conductor,
        rho = 1.7241e-8,
        eps_r = 1.0,
        mu_r = 1.0,
        T0 = 20.0,
        alpha = 0.00393
    )
    dielectric = PB.Material(
        kind = :insulator,
        rho = Grid(2.0e11, 10.0),
        eps_r = Grid(2.4, 5.0)
    )
    outer_dielectric = PB.Material(
        kind = :insulator,
        rho = 1.0e14,
        eps_r = 2.3
    )
    function build_design(resolved_dielectric)
        core_outer = 0.010
        insulation_outer = 0.017
        sheath_outer = 0.018
        jacket_outer = 0.020
        return build(
            CableDesign,
            "sobol-cable-smoke",
            Stack(AbstractCablePart[
                Group(
                    :core,
                    Region(
                        :core_conductor,
                        Annulus(0.0, core_outer),
                        copper
                    )
                ),
                Region(
                    :core_insulation,
                    Annulus(core_outer, insulation_outer),
                    resolved_dielectric
                ),
                Group(
                    :sheath,
                    Region(
                        :sheath_conductor,
                        Annulus(insulation_outer, sheath_outer),
                        copper
                    )
                ),
                Region(
                    :outer_insulation,
                    Annulus(sheath_outer, jacket_outer),
                    outer_dielectric
                )
            ])
        )
    end

    design = Gridspace{CableDesign}(build_design, (dielectric,))
    system = build(
        LineCableSystem,
        design,
        Pose2(0.0, -1.0, 0.0);
        connections = Dict(:core => 1, :sheath => 2),
        system_id = "sobol-cable-smoke",
        line_length = 1000.0
    )
    problem_space = LineParametersProblem(
        system;
        temperature = 20.0,
        earth_props = PB.Earth(rho = 100.0, eps_r = 10.0),
        frequencies = [50.0]
    )
    problem = ParametricProblem(problem_space, (trace = true,))
    inner = Formulation(;
        insulation_admittance = formula(:Marti2001),
        options = (
            reduce_bundle = false,
            kron_reduction = false,
            ideal_transposition = true
        )
    )

    point = only(collect(PB.points(problem.space)))
    descriptors = PB.uncertainties(point)
    @test length(descriptors) == 2
    nominal_problem = PB.realize(point, map(PB.nominal, descriptors))
    direct_before = compute(nominal_problem, inner; options = problem.options)

    formulation = Sensitivity(
        inner,
        GlobalSensitivity.Sobol(),
        (@observe(G[1, 1, 1]), @observe(C[1, 1, 1]));
        samples = 16,
        distribution = :uniform,
        input_labels = ("dielectric resistivity", "relative permittivity"),
        max_evaluations = 64,
        max_output_values = 128,
        options = (retain_details = true,)
    )
    first_result = compute(problem, formulation)
    second_result = compute(problem, formulation)
    @test first_result.values == second_result.values
    @test details(first_result) == details(second_result)

    product = only(first_result)
    @test product.evaluations == 64
    @test product.inputs == ["dielectric resistivity", "relative permittivity"]
    @test product.active_indices == [1, 2]
    @test length(product.first_order) == length(product.total_order) == 2
    @test size.(product.first_order) == ((2,), (2,))
    @test size.(product.total_order) == ((2,), (2,))
    @test all(values -> all(isfinite, values), product.first_order)
    @test all(values -> all(isfinite, values), product.total_order)

    records = only(details(first_result).evaluations)
    @test length(records) == product.evaluations
    @test all(record -> keys(record) == (:trace,), records)
    @test all(record -> size(record.trace.Z) == (2, 2, 1), records)
    @test all(record -> size(record.trace.P) == (2, 2, 1), records)

    direct_after = compute(nominal_problem, inner; options = problem.options)
    @test direct_after.Z.values == direct_before.Z.values
    @test direct_after.Y.values == direct_before.Y.values
    @test details(direct_after) == details(direct_before)
end
