@testitem "Extensions / Measurements boundary / unloaded uncertainty APIs fail explicitly" tags=[
    :extension,
    :core_only
] begin
    import LineCableModels

    @test Base.get_extension(LineCableModels, :LineCableModelsMeasurementsExt) === nothing
    @test !isdefined(LineCableModels, :Utils)
    @test LineCableModels.nominal(1.0) == 1.0
    @test LineCableModels.nominal(1.0 + 2.0im) == 1.0 + 2.0im
    @test LineCableModels.standard_uncertainty(1.0) == 0.0

    encoded_measurement=Dict(
        "__type__"=>"Measurement",
        "value"=>Dict("__type__"=>"Float", "value"=>1.0),
        "uncertainty"=>Dict("__type__"=>"Float", "value"=>0.1)
    )
    @test_throws ArgumentError LineCableModels.ImportExport._deserialize_value(
        encoded_measurement
    )
end

@testitem "Measurements / natural promotion preserves covariance" tags=[:extension] setup=[
    DataModelTestSupport, UseDataModelSupport, TestNumerics] begin
    using Measurements
    import Measurements: derivative

    @testset "Owned constructors preserve primitive dependencies" begin
        x = measurement(2.0, 0.1)
        material = Material(x, 1.0f0, 1, 20, 0.004)

        @test eltype(material) === Measurement{Float64}
        @test derivative(material.rho, x) ≈ 1.0
        @test nominal(material.rho) == value(x)
        @test standard_uncertainty(material.rho) == uncertainty(x)
    end

    @testset "Mixed BaseParams calls retain primitive sensitivities" begin
        r_ex = measurement(0.02, 0.001)
        rho = measurement(1.7241e-8, 1.0e-9)
        resistance = tubular_resistance(0.0, r_ex, rho, 0.0, 20.0, 20.0)

        @test !iszero(derivative(resistance, r_ex))
        @test !iszero(derivative(resistance, rho))
    end
end

@testitem "DataModel / Measurements / covariance through complete assembly" tags=[:extension] setup=[
    DataModelTestSupport, UseDataModelSupport, TestNumerics] begin
    using Measurements
    import Measurements: derivative
    copper_props=Material(1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
    insulator_props=Material(1.0e14, 2.3, 1.0, 20.0, 0.0)
    semicon_props=Material(1000.0, 1000.0, 1.0, 20.0, 0.0)

    diameter=measurement(0.02, 0.001)
    insulation_thickness=measurement(0.01, 0.001)
    semicon_thickness=measurement(0.005, 0.0005)

    core=Tubular(0.0, diameter/2, copper_props)
    insulators=InsulatorGroup(
        Insulator(core.r_ex, core.r_ex+insulation_thickness, insulator_props),
    )

    @test insulators.r_in === core.r_ex
    @test iszero(uncertainty(core.r_ex - insulators.r_in))

    insulators=add!(insulators, Semicon(
        insulators.r_ex,
        insulators.r_ex+semicon_thickness,
        semicon_props
    ))

    expected_radius=diameter/2+insulation_thickness+semicon_thickness
    @test value(insulators.r_ex) ≈ value(expected_radius)
    @test uncertainty(insulators.r_ex) ≈ uncertainty(expected_radius)
    @test derivative(insulators.r_ex, diameter) ≈ 0.5
    @test derivative(insulators.r_ex, insulation_thickness) ≈ 1.0
    @test derivative(insulators.r_ex, semicon_thickness) ≈ 1.0

    conductors=ConductorGroup(core)
    component=CableComponent("core", conductors, insulators)
    design=CableDesign("uq-path", component)
    position=CablePosition(design, 0.0, -1.0, Dict("core"=>1))
    system=LineCableSystem("uq-path", 1000.0, position)
    assembled_radius=system.cables[1].design_data.components[1].insulator_group.r_ex

    @test system.cables[1].design_data === design
    @test derivative(assembled_radius, diameter) ≈ 0.5
    @test derivative(assembled_radius, insulation_thickness) ≈ 1.0
    @test derivative(assembled_radius, semicon_thickness) ≈ 1.0
end

@testitem "Measurements / external boundaries and numerical kernels" tags=[:extension] setup=[
    EngineTestSupport, UseEngineSupport, TestNumerics] begin
    using Measurements
    using SpecialFunctions
    import LineCableModels

    extension_module=Base.get_extension(LineCableModels, :LineCableModelsMeasurementsExt)
    @test extension_module !== nothing

    uncertain=measurement(-20.0, 1.0)
    @test value(uncertain) == -20.0
    @test uncertainty(uncertain) == 1.0
    @test LineCableModels.nominal(uncertain) == -20.0
    @test LineCableModels.standard_uncertainty(uncertain) == 1.0

    clipped=LineCableModels.Engine._clip_field(measurement(1.0e-12, 2.0e-12), 1.0e-9)
    @test value(clipped) == 0.0
    @test uncertainty(clipped) == 0.0
    retained=LineCableModels.Engine._clip_field(measurement(2.0, 0.25), 1.0e-9)
    @test value(retained) == 2.0
    @test uncertainty(retained) == 0.25

    encoded=LineCableModels.ImportExport._serialize_value(uncertain)
    @test encoded["__type__"] == "Measurement"
    decoded=LineCableModels.ImportExport._deserialize_value(encoded)
    @test value(decoded) == value(uncertain)
    @test uncertainty(decoded) == uncertainty(uncertain)
    @test LineCableModels.ImportExport.stringify(uncertain) == "-20 ± 1"

    argument=complex(measurement(1.25, 0.01), measurement(0.5, 0.02))
    nominal=complex(1.25, 0.5)
    for function_value in (besselix, besselkx)
        propagated=function_value(0, argument)
        @test value(real(propagated)) ≈ real(function_value(0, nominal)) rtol=1.0e-11
        @test value(imag(propagated)) ≈ imag(function_value(0, nominal)) rtol=1.0e-11
        @test uncertainty(real(propagated)) > 0
        @test uncertainty(imag(propagated)) > 0
    end

    singleton=extension_module._joint_coordinates(reshape([3.0, 4.0], 2, 1))
    @test value.(singleton) == [3.0, 4.0]
    @test all(iszero, uncertainty.(singleton))

    coordinates=extension_module._joint_coordinates([1.0 2.0 3.0; 2.0 4.0 6.0])
    @test value.(coordinates) ≈ [2.0, 4.0]
    @test uncertainty(coordinates[1]) ≈ 1.0
    @test uncertainty(coordinates[2] - 2 * coordinates[1]) ≤ 10 * eps(Float64)
end

@testitem "Measurements / Monte Carlo / retained and marginal reconstruction" tags=[:extension] setup=[
    EngineTestSupport, UseEngineSupport, TestNumerics] begin
    using Measurements
    import LineCableModels.ParametricBuilder: AbsoluteUncertainty, UncertainValue

    samples=CableConstants(
        [1.0, 2.0, 3.0],
        [10.0, 20.0, 30.0],
        [100.0, 200.0, 300.0]
    )
    nominal=CableConstants(2.0, 20.0, 200.0)
    deviations=CableConstants(1.0, 10.0, 100.0)
    retained=MonteCarloResult(
        nominal,
        nothing,
        samples,
        nothing,
        UncertainValue(nominal, deviations, AbsoluteUncertainty()),
        3,
        0.95,
        0.02,
        :normal,
        UInt64(9),
        (hash = "measurement-retained",)
    )
    correlated=Measurements.measurement(retained)
    @test value(correlated.R) ≈ 2.0
    @test uncertainty(correlated.R) ≈ 1.0
    @test uncertainty(correlated.L - 10 * correlated.R) ≤ 10 * eps(Float64)

    marginal=MonteCarloResult(
        nominal,
        nothing,
        nothing,
        nothing,
        UncertainValue(nominal, deviations, AbsoluteUncertainty()),
        3,
        0.95,
        0.02,
        :normal,
        UInt64(10),
        (hash = "measurement-marginal",)
    )
    reconstructed=@test_logs (
        :warn,
        r"discarding output correlations"
    ) Measurements.measurement(marginal)
    @test value(reconstructed.R) == nominal.R
    @test value(reconstructed.L) == nominal.L
    @test value(reconstructed.C) == nominal.C
    @test uncertainty(reconstructed.R) == deviations.R
    @test Measurements.cov(reconstructed.R, reconstructed.L) == 0.0

    frequencies=[50.0, 100.0]
    resistance=reshape([1.0, 2.0], 1, 1, :)
    inductance=reshape([1.0e-3, 2.0e-3], 1, 1, :)
    conductance=reshape([1.0e-6, 2.0e-6], 1, 1, :)
    capacitance=reshape([1.0e-9, 2.0e-9], 1, 1, :)
    angular=reshape(2π .* frequencies, 1, 1, :)
    line_nominal=LineParameters(
        complex.(resistance, inductance .* angular),
        complex.(conductance, capacitance .* angular),
        frequencies
    )
    line_deviations=RLCG(
        fill(0.1, 1, 1, 2),
        fill(1.0e-4, 1, 1, 2),
        fill(1.0e-10, 1, 1, 2),
        fill(1.0e-7, 1, 1, 2)
    )
    line_marginal=MonteCarloResult(
        line_nominal,
        nothing,
        nothing,
        nothing,
        UncertainValue(line_nominal, line_deviations, AbsoluteUncertainty()),
        3,
        0.95,
        0.02,
        :normal,
        UInt64(11),
        (hash = "line-marginal",)
    )
    line_reconstructed=@test_logs (
        :warn,
        r"discarding output correlations"
    ) Measurements.measurement(line_marginal)
    @test value(real(Z(line_reconstructed, 1, 1, 1))) == resistance[1]
    @test uncertainty(real(Z(line_reconstructed, 1, 1, 1))) == line_deviations.R[1]
    @test uncertainty(imag(Y(line_reconstructed, 1, 1, 2))) ≈
          angular[2] * line_deviations.C[2]
end
