@testitem "Engine / line parameters / RMS comparison" tags=[:unit] setup=[EngineTestSupport] begin
    using Test
    using LineCableModels
    using LineCableModels.Engine

    frequencies_value=[1.0, 10.0]
    impedance=reshape(ComplexF64[1 + 2im, 2 + 3im], 1, 1, :)
    admittance=reshape(ComplexF64[2im, 4im], 1, 1, :)
    reference=LineParameters(PhaseDomain, impedance, admittance, frequencies_value)

    identical=compare(reference, reference)
    @test identical.Z == RMSError(0.0, 0.0)
    @test identical.Y == RMSError(0.0, 0.0)

    candidate=LineParameters(
        PhaseDomain,
        impedance .+ (1+0im),
        admittance .+ (0+2im),
        frequencies_value
    )
    comparison=compare(reference, candidate)
    @test comparison.Z.absolute ≈ 1.0
    @test comparison.Z.relative ≈ sqrt(2 / sum(abs2, impedance))
    @test comparison.Y.absolute ≈ 2.0
    @test comparison.Y.relative ≈ sqrt(8 / sum(abs2, admittance))

    zero_reference=LineParameters(
        PhaseDomain,
        zeros(ComplexF64, 1, 1, 2),
        zeros(ComplexF64, 1, 1, 2),
        frequencies_value
    )
    @test compare(zero_reference, zero_reference).Z.relative == 0.0
    @test isinf(compare(zero_reference, candidate).Z.relative)

    different_frequency=LineParameters(
        PhaseDomain, impedance, admittance, [1.0, 11.0]
    )
    @test_throws ArgumentError compare(reference, different_frequency)

    different_size=LineParameters(
        PhaseDomain,
        zeros(ComplexF64, 2, 2, 2),
        zeros(ComplexF64, 2, 2, 2),
        frequencies_value
    )
    @test_throws DimensionMismatch compare(reference, different_size)

    total=LineParameters(
        PhaseDomain, impedance, admittance, frequencies_value; basis = :total
    )
    @test_throws ArgumentError compare(reference, total)

    modal=LineParameters(ModalDomain, impedance, admittance, frequencies_value)
    @test_throws ArgumentError compare(reference, modal)
end
