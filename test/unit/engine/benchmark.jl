@testitem "Engine / line parameters / RMS benchmark" tags=[:unit] setup=[EngineTestSupport] begin
    using Test
    using DataFrames
    using LineCableModels
    using LineCableModels.Engine

    frequencies_value=[1.0, 10.0]
    impedance=Array{ComplexF64, 3}(undef, 2, 2, 2)
    admittance=similar(impedance)
    impedance[:, :, 1]=[1+2im 2+0im; 3+0im 4+1im]
    impedance[:, :, 2]=[2+3im 4+0im; 6+0im 8+2im]
    admittance[:, :, 1]=[2im 1im; 1im 4im]
    admittance[:, :, 2]=[4im 2im; 2im 8im]
    reference=LineParameters(PhaseDomain, impedance, admittance, frequencies_value)

    identical=compare(reference, reference)
    @test identical.Z.absolute == zeros(2, 2)
    @test identical.Z.relative == zeros(2, 2)
    @test identical.Y.absolute == zeros(2, 2)
    @test identical.Y.relative == zeros(2, 2)

    candidate=LineParameters(
        PhaseDomain,
        impedance .+ reshape(ComplexF64[1, 2, 3, 4], 2, 2, 1),
        admittance .+ reshape(ComplexF64[2im, 3im, 4im, 5im], 2, 2, 1),
        frequencies_value
    )
    comparison=compare(reference, candidate)
    comparison_observables=observables(
        comparison,
        (
            (Z, absolute_error),
            (Z, relative_error),
            (Y, absolute_error),
            (Y, relative_error)
        )
    )
    @test comparison_observables isa Tuple
    @test comparison_observables[1].values ==
          1_000comparison.Z.absolute
    @test comparison_observables[2].values ==
          comparison.Z.relative
    @test comparison_observables[3].values ==
          1_000comparison.Y.absolute
    @test comparison_observables[4].values ==
          comparison.Y.relative
    @test comparison_observables[1].values !==
          comparison.Z.absolute
    @test comparison_observables[4].values !==
          comparison.Y.relative
    @test comparison.Z.absolute ≈ [1.0 3.0; 2.0 4.0]
    @test comparison.Z.relative[1, 1] ≈
          sqrt(2 / sum(abs2, impedance[1, 1, :]))
    @test comparison.Z.relative[2, 2] ≈
          sqrt(32 / sum(abs2, impedance[2, 2, :]))
    @test comparison.Y.absolute ≈ [2.0 4.0; 3.0 5.0]
    @test comparison.Y.relative[1, 1] ≈
          sqrt(8 / sum(abs2, admittance[1, 1, :]))

    zero_reference=LineParameters(
        PhaseDomain,
        zeros(ComplexF64, 2, 2, 2),
        zeros(ComplexF64, 2, 2, 2),
        frequencies_value
    )
    @test compare(zero_reference, zero_reference).Z.relative == zeros(2, 2)
    @test all(isinf, compare(zero_reference, candidate).Z.relative)

    different_frequency=LineParameters(
        PhaseDomain, impedance, admittance, [1.0, 11.0]
    )
    @test_throws ArgumentError compare(reference, different_frequency)

    different_size=LineParameters(
        PhaseDomain,
        zeros(ComplexF64, 3, 3, 2),
        zeros(ComplexF64, 3, 3, 2),
        frequencies_value
    )
    @test_throws DimensionMismatch compare(reference, different_size)

    total=LineParameters(
        PhaseDomain, impedance, admittance, frequencies_value; basis = :total
    )
    @test_throws ArgumentError compare(reference, total)

    modal_source=LineParameters(
        PhaseDomain,
        cat([1+2im 0; 0 4+1im], [2+3im 0; 0 8+2im]; dims = 3),
        cat([2im 0; 0 4im], [4im 0; 0 8im]; dims = 3),
        frequencies_value
    )
    modal=compute(
        ModalTransformationProblem(modal_source),
        ModalTransformationFormulation(:Fortescue)
    )
    @test_throws ArgumentError compare(reference, modal)

    displayed=LineParametersBenchmark(
        RMSError([2.0e-20 1.0e-5; 2.0e-5 3.0e-5], [0.5 0.01; 0.02 0.03]),
        RMSError([3.0e-20 2.0e-8; 3.0e-8 4.0e-8], [0.4 0.001; 0.002 0.003])
    )
    table=DataFrame(displayed)
    @test names(table) == ["row", "column", "ΔZ", "εZ", "ΔY", "εY"]
    z_noise=only(eachrow(table[(table.row .== 1) .& (table.column .== 1), :]))
    y_signal=only(eachrow(table[(table.row .== 1) .& (table.column .== 2), :]))
    @test z_noise.ΔZ == 0.0
    @test z_noise.εZ == 0.5
    @test y_signal.εY == 0.001
    @test displayed.Z.relative[1, 1] == 0.5
    @test DataFrame(displayed; clip = false).ΔZ[1] ≈ 2.0e-17
    observed=LineCableModels.ReportBuilder.observation_columns(table)
    @test keys(observed) == (:ΔZ, :εZ, :ΔY, :εY)
    @test LineCableModels.Units.label(observed.ΔZ.unit) == "Ω/km"
    @test LineCableModels.Units.label(observed.εZ.unit) == ""
end
