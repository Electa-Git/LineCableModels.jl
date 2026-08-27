@testitem "DataModel–Engine / operating temperature is problem-owned and isolated" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestFixtures,
    TestNumerics
] begin
    using LinearAlgebra

    system=TestFixtures.three_phase_system()
    earth=EarthModel(100.0, 10.0, 1.0)
    frequencies=[50.0, 500.0]
    problem_at(temperature) = LineParametersProblem(
        system; temperature, earth_props = earth, frequencies
    )
    formulation=Formulation(:analytical; options = (ideal_transposition = false,))
    uncorrected=Formulation(
        Val(:analytical);
        options = (ideal_transposition = false, temperature_correction = false)
    )

    design=first(system.cables).design_data
    eager_before=[(
                      component.conductor_props.rho,
                      component.conductor_group.resistance,
                      component.insulator_group.shunt_capacitance,
                      component.insulator_group.shunt_conductance
                  ) for component in design.components]
    cold_problem=problem_at(20.0)
    hot_problem=problem_at(90.0)
    cold=compute(cold_problem, formulation)
    hot=compute(hot_problem, formulation)
    repeated=compute(cold_problem, formulation)

    @test hot_problem.temperature == 90.0
    @test all(!hasfield(typeof(layer), :temperature)
    for component in design.components
    for layer in (
        component.conductor_group.layers...,
        component.insulator_group.layers...
    ))
    @test sum(diag(real.(hot.Z.values[:, :, 1]))) >
          sum(diag(real.(cold.Z.values[:, :, 1])))
    @test repeated.Z.values == cold.Z.values
    @test repeated.Y.values == cold.Y.values

    cold_input=LineCableModels.Engine.AnalyticalInput(cold_problem)
    hot_input=LineCableModels.Engine.AnalyticalInput(hot_problem)
    cold_rho=LineCableModels.Engine._operating_resistivity(
        cold_input, cold_problem, formulation
    )
    hot_rho=LineCableModels.Engine._operating_resistivity(
        hot_input, hot_problem, formulation
    )
    @test cold_rho ≈ cold_input.rho0_cond
    @test hot_rho ≈
          hot_input.rho0_cond .* [temperature_factor(alpha, 90.0, reference)
           for (alpha, reference) in zip(hot_input.alpha_cond, hot_input.T0_cond)]

    cold_base=compute(cold_problem, uncorrected)
    hot_base=compute(hot_problem, uncorrected)
    @test hot_base.Z.values == cold_base.Z.values
    @test hot_base.Y.values == cold_base.Y.values
    @test eager_before == [(
               component.conductor_props.rho,
               component.conductor_group.resistance,
               component.insulator_group.shunt_capacitance,
               component.insulator_group.shunt_conductance
           ) for component in design.components]
end

@testitem "DataModel–Engine / immutable baseline parity and transposition delta" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestFixtures,
    TestNumerics
] begin
    using JSON3
    using LinearAlgebra

    baseline=JSON3.read(read(
        joinpath(
            pkgdir(LineCableModels), "test", "fixtures", "reference",
            "datamodel_engine_baseline.json"
        ),
        String))
    @test baseline["baseline_sha"] == "a71bdfe1ac832f27a0c88b1d02596194aac46ec7"

    tensor(node) = reshape(
        complex.(Float64.(node["real"]), Float64.(node["imag"])),
        Tuple(Int.(node["size"]))
    )
    design=TestFixtures.mv_cable_design()
    constants=compute(CableConstantsProblem(design), Formulation())
    @test TestNumerics.isapprox_scaled(constants.R, baseline["cable_constants"]["R"])
    @test TestNumerics.isapprox_scaled(constants.L, baseline["cable_constants"]["L"])
    @test TestNumerics.isapprox_scaled(constants.C, baseline["cable_constants"]["C"])

    for (component, expected) in zip(design.components, baseline["components"])
        @test TestNumerics.isapprox_scaled(
            component.conductor_group.resistance,
            expected["conductor_group"]["resistance"]
        )
        @test TestNumerics.isapprox_scaled(
            component.conductor_props.rho, expected["conductor_rho"]
        )
        @test TestNumerics.isapprox_scaled(
            component.insulator_group.shunt_capacitance,
            expected["insulator_group"]["capacitance"]
        )
        @test TestNumerics.isapprox_scaled(
            component.insulator_group.shunt_conductance,
            expected["insulator_group"]["conductance"]
        )
    end

    problem=TestFixtures.line_parameters_problem(
        frequencies = Float64.(baseline["frequencies"]),
    )
    ordinary=compute(problem, Formulation(
        Val(:analytical); options = (ideal_transposition = false,)
    ))
    transposed=compute(problem, Formulation(
        Val(:analytical); options = (ideal_transposition = true,)
    ))

    # The baseline applied the flag to opposite matrix families. These crossed
    # references preserve every unchanged value while recording the approved fix.
    old_false=baseline["formulations"]["reduced_option_false"]
    old_true=baseline["formulations"]["reduced_option_true"]
    @test TestNumerics.isapprox_scaled(ordinary.Z.values, tensor(old_true["Z"]))
    @test TestNumerics.isapprox_scaled(ordinary.Y.values, tensor(old_false["Y"]))
    @test TestNumerics.isapprox_scaled(transposed.Z.values, tensor(old_false["Z"]))
    @test TestNumerics.isapprox_scaled(transposed.Y.values, tensor(old_true["Y"]))

    for values in (transposed.Z.values, transposed.Y.values)
        for frequency in axes(values, 3)
            slice=@view values[:, :, frequency]
            @test all(diag(slice) .≈ first(diag(slice)))
            off_diagonal=[slice[row, column]
                          for column in axes(slice, 2), row in axes(slice, 1)
                          if row!=column]
            @test all(off_diagonal .≈ first(off_diagonal))
        end
    end
end

@testitem "DataModel–Engine / scalar precision and external Bessel boundary" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestFixtures,
    TestNumerics
] begin
    import LineCableModels.ParametricBuilder as PB
    using SpecialFunctions

    function compact_problem(::Type{T}) where {T <: AbstractFloat}
        conductor=PB.Material(
            kind = :conductor,
            rho = T(0.01), eps_r = one(T), mu_r = one(T), T0 = T(20), alpha = T(0.004)
        )
        dielectric=PB.Material(
            kind = :insulator,
            rho = T(100), eps_r = T(2), mu_r = one(T), T0 = T(20), alpha = zero(T)
        )
        design=PB.CableBuilder(
            "precision-cable",
            PB.Conductor.Solid(:core; radius = T(0.10), material = conductor),
            PB.Insulator.Tubular(:core; thickness = T(0.05), material = dielectric)
        )
        specification=PB.SystemBuilder(
            "precision-system",
            design,
            PB.at(x = zero(T), y = -one(T), phases = :core=>1);
            length = T(100),
            temperature = T(20),
            earth = PB.Earth(rho = T(100), eps_r = T(10), mu_r = one(T)),
            frequencies = T[50]
        )
        return specification
    end

    half=compact_problem(Float16)
    @test eltype(half) === Float32
    @test eltype(first(half.system.cables).design_data) === Float32
    @test eltype(compute(half, Formulation())) === ComplexF32

    for T in (Float32, Float64, BigFloat)
        problem=compact_problem(T)
        result=compute(problem, Formulation())
        @test eltype(problem) === T
        @test eltype(result) === Complex{T}
        @test all(isfinite, result.Z.values)
        @test all(isfinite, result.Y.values)
    end

    for order in 0:1
        value=Complex{BigFloat}(1, 2)
        float_value=ComplexF64(value)
        @test ComplexF64(LineCableModels.Engine.special_besselix(order, value)) ≈
              SpecialFunctions.besselix(order, float_value)
        @test ComplexF64(LineCableModels.Engine.special_besselkx(order, value)) ≈
              SpecialFunctions.besselkx(order, float_value)
    end
end

@testitem "Engine / scoped verbosity preserves the caller logger" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestFixtures
] begin
    using Logging
    import NLsolve

    options=LineCableModels.computation_options(Val(AnalyticalFormulation),
        (verbosity = (
            default = 0, LineCableModels = 2, NLsolve = 1, QuadGK = 0
        ),))
    @test verbosity(options, :LineCableModels) == 2
    @test verbosity(options, :NLsolve) == 1
    @test verbosity(options, :unlisted) == 0

    sink=IOBuffer()
    logger=LineCableModels.Engine.ConsoleVerbosityLogger(
        SimpleLogger(sink, Logging.Debug), options.verbosity
    )
    @test Logging.shouldlog(logger, Logging.Debug, LineCableModels, :test, :debug)
    @test Logging.shouldlog(logger, Logging.Info, NLsolve, :test, :info)
    @test !Logging.shouldlog(logger, Logging.Debug, NLsolve, :test, :debug)
    @test !Logging.shouldlog(logger, Logging.Info, Base, :test, :info)
    @test Logging.shouldlog(logger, Logging.Warn, Base, :test, :warn)

    caller=current_logger()
    compute(TestFixtures.line_parameters_problem(), Formulation())
    @test current_logger() === caller
    @test_throws ArgumentError LineCableModels.computation_options(
        Val(AnalyticalFormulation), (
            verbosity = (LineCableModels = 1,),
        ))
    @test_throws ArgumentError LineCableModels.computation_options(
        Val(AnalyticalFormulation), (
            verbosity = (default = 3,),
        ))
end
