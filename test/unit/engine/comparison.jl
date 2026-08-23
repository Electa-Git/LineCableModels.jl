@testitem "Engine / line parameters / RMS comparison" tags=[:unit] setup=[EngineTestSupport] begin
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

    modal=LineParameters(ModalDomain, impedance, admittance, frequencies_value)
    @test_throws ArgumentError compare(reference, modal)

    displayed=LineParametersBenchmark(
        RMSError([2.0e-20 1.0e-5; 2.0e-5 3.0e-5], [0.5 0.01; 0.02 0.03]),
        RMSError([3.0e-20 2.0e-8; 3.0e-8 4.0e-8], [0.4 0.001; 0.002 0.003])
    )
    table=DataFrame(displayed; zero_atol = (Z = 1.0e-9, Y = 1.0e-12))
    z_noise=only(eachrow(table[(table.quantity .== :Z) .& (table.row .== 1) .& (table.column .== 1), :]))
    y_signal=only(eachrow(table[(table.quantity .== :Y) .& (table.row .== 1) .& (table.column .== 2), :]))
    @test ismissing(z_noise.rms_relative)
    @test z_noise.relative_status === :below_absolute_floor
    @test y_signal.rms_relative == 0.001
    @test y_signal.relative_status === :reported
    @test displayed.Z.relative[1, 1] == 0.5
    @test metadata(table, "zero_atol") == (Z = 1.0e-9, Y = 1.0e-12)
    @test isequal(
        DataFrame(displayed; zero_atol = (Y = 1.0e-12, Z = 1.0e-9)),
        table
    )
    @test_throws ArgumentError DataFrame(displayed; zero_atol = (Z = 0.0,))
    @test_throws ArgumentError DataFrame(
        displayed; zero_atol = (Z = -1.0, Y = 0.0)
    )
    @test_throws ArgumentError DataFrame(
        displayed; zero_atol = (Z = 0.0, Y = Inf)
    )
end
