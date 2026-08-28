@testitem "Extensions / Measurements boundary / unloaded uncertainty APIs fail explicitly" tags=[
    :extension,
    :core_only
] begin
    import LineCableModels

    @test Base.get_extension(LineCableModels, :LineCableModelsMeasurementsExt) === nothing
    @test !isdefined(LineCableModels, :Utils)
    @test LineCableModels.nominal(1.0) == 1.0
    @test LineCableModels.nominal(1.0 + 2.0im) == 1.0 + 2.0im
    @test LineCableModels.uncertainty(1.0) == 0.0

    encoded_measurement=Dict(
        "__type__"=>"Measurement",
        "value"=>Dict("__type__"=>"Float", "value"=>1.0),
        "uncertainty"=>Dict("__type__"=>"Float", "value"=>0.1)
    )
    @test_throws ArgumentError LineCableModels.ImportExport.deserialize_value(
        encoded_measurement
    )
end

@testitem "Measurements / natural promotion preserves covariance" tags=[:extension] setup=[
    DataModelTestSupport, UseDataModelSupport, TestNumerics] begin
    using Measurements
    import Measurements: derivative

    @testset "Owned constructors preserve primitive dependencies" begin
        x = measurement(2.0, 0.1)
        material = Material(:conductor, x, 1.0f0, 1, 20, 0.004)

        @test eltype(material) === Measurement{Float64}
        @test derivative(material.rho, x) ≈ 1.0
        @test nominal(material.rho) == value(x)
        @test uncertainty(material.rho) == uncertainty(x)
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
    copper_props=Material(:conductor, 1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
    insulator_props=Material(:insulator, 1.0e14, 2.3, 1.0, 20.0, 0.0)
    semicon_props=Material(:semicon, 1000.0, 1000.0, 1.0, 20.0, 0.0)

    diameter=measurement(0.02, 0.001)
    insulation_thickness=measurement(0.01, 0.001)
    semicon_thickness=measurement(0.005, 0.0005)

    core=Region(:core_metal, Disk(diameter/2), copper_props)
    insulation=Region(:insulation, Shell(insulation_thickness), insulator_props)
    semicon=Region(:semicon, Shell(semicon_thickness), semicon_props)
    root=Stack(Group(:core, core), insulation, semicon)
    resolved=resolve(EmptyBoundary(), root)
    insulation_shape=resolved.regions[2].primitive
    semicon_shape=resolved.regions[3].primitive

    @test r_in(insulation_shape) === r_ex(resolved.regions[1].primitive)
    @test iszero(uncertainty(
        r_ex(resolved.regions[1].primitive) - r_in(insulation_shape)
    ))

    expected_radius=diameter/2+insulation_thickness+semicon_thickness
    @test value(r_ex(semicon_shape)) ≈ value(expected_radius)
    @test uncertainty(r_ex(semicon_shape)) ≈ uncertainty(expected_radius)
    @test derivative(r_ex(semicon_shape), diameter) ≈ 0.5
    @test derivative(r_ex(semicon_shape), insulation_thickness) ≈ 1.0
    @test derivative(r_ex(semicon_shape), semicon_thickness) ≈ 1.0

    design=build(CableDesign, "uq-path", root)
    system=build(
        LineCableSystem,
        design,
        Pose2(0.0, -1.0, 0.0);
        connections = Dict(:core=>1),
        system_id = "uq-path",
        line_length = 1000.0
    )
    assembled_radius=r_ex(system.geometry[end].primitive)

    @test system.designs[1] === design
    @test derivative(assembled_radius, diameter) ≈ 0.5
    @test derivative(assembled_radius, insulation_thickness) ≈ 1.0
    @test derivative(assembled_radius, semicon_thickness) ≈ 1.0
end

@testitem "Measurements / external boundaries and numerical kernels" tags=[:extension] setup=[
    EngineTestSupport, UseEngineSupport, TestNumerics] begin
    using Measurements
    using Statistics
    using SpecialFunctions
    import LineCableModels

    extension_module=Base.get_extension(LineCableModels, :LineCableModelsMeasurementsExt)
    @test extension_module !== nothing
    @test any(
        method -> method.module === extension_module,
        methods(LineCableModels.Grammar.detach)
    )
    @test any(
        method -> method.module === extension_module,
        methods(LineCableModels.ReportBuilder.encode_cell)
    )

    uncertain=measurement(-20.0, 1.0)
    @test value(uncertain) == -20.0
    @test uncertainty(uncertain) == 1.0
    @test LineCableModels.nominal(uncertain) == -20.0
    @test LineCableModels.uncertainty(uncertain) == 1.0

    clipped=LineCableModels.Grammar.detach(
        measurement(eps(Float64)/2, eps(Float64)/4),
        1.0,
        true
    )
    @test value(clipped) == 0.0
    @test uncertainty(clipped) == 0.0
    clipped_array=LineCableModels.Grammar.detach(
        [measurement(eps(Float64)/2, eps(Float64)/4)],
        1.0,
        true
    )
    @test value(only(clipped_array)) == 0.0
    @test uncertainty(only(clipped_array)) == 0.0
    retained=LineCableModels.Grammar.detach(
        measurement(2.0, 0.25),
        1.0,
        true
    )
    @test value(retained) == 2.0
    @test uncertainty(retained) == 0.25

    encoded=LineCableModels.ImportExport.serialize_value(uncertain)
    @test encoded["__type__"] == "Measurement"
    decoded=LineCableModels.ImportExport.deserialize_value(encoded)
    @test value(decoded) == value(uncertain)
    @test uncertainty(decoded) == uncertainty(uncertain)
    @test LineCableModels.ReportBuilder.encode_cell(
        LineCableModels.ReportBuilder.XLSXReportDefinition(),
        uncertain
    ) == "-20 ± 1"

    argument=complex(measurement(1.25, 0.01), measurement(0.5, 0.02))
    nominal=complex(1.25, 0.5)
    for function_value in (besselix, besselkx)
        propagated=function_value(0, argument)
        @test value(real(propagated)) ≈ real(function_value(0, nominal)) rtol=1.0e-11
        @test value(imag(propagated)) ≈ imag(function_value(0, nominal)) rtol=1.0e-11
        @test uncertainty(real(propagated)) > 0
        @test uncertainty(imag(propagated)) > 0
    end

    @test !isdefined(extension_module, :_joint_coordinates)
end

@testitem "Measurements / Monte Carlo / no implicit completed-result conversion" tags=[:extension] setup=[
    EngineTestSupport, UseEngineSupport, TestNumerics] begin
    using Measurements
    using Statistics
    formulation=MonteCarlo(Formulation(); trials = 3, seed = 9,
        return_samples = true)
    values=[CableConstants(2.0, 20.0, 200.0)]
    stats=[(
        R = SampleSummary(2.0, 1.0, 1.0, 1.1, 2.0, 2.9, 3.0, 3),
        L = SampleSummary(20.0, 10.0, 10.0, 11.0, 20.0, 29.0, 30.0, 3),
        C = SampleSummary(200.0, 100.0, 100.0, 110.0, 200.0, 290.0, 300.0, 3)
    )]
    sample_values=[(
        R = [1.0, 2.0, 3.0],
        L = [10.0, 20.0, 30.0],
        C = [100.0, 200.0, 300.0]
    )]
    completed=MonteCarloResult(
        formulation, values, stats, sample_values, nothing,
        UInt64(9), UInt64[9], [3]
    )

    @test !applicable(Measurements.measurement, completed)
    retained=only(samples(completed))
    sample_matrix=vcat(
        reshape(retained.R, 1, :),
        reshape(retained.L, 1, :),
        reshape(retained.C, 1, :)
    )
    means=vec(mean(sample_matrix; dims = 2))
    centered=sample_matrix .- means
    covariance=centered*transpose(centered)/(size(sample_matrix, 2)-1)
    coordinates=Measurements.correlated_values(means, covariance)
    @test value.(coordinates) ≈ [2.0, 20.0, 200.0]
    @test uncertainty(coordinates[2] - 10 * coordinates[1]) ≤ 20 * eps(Float64)
end
