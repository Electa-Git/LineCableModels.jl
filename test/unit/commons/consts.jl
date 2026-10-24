@testitem "Commons / vacuum permittivity is evaluated in the requested type" tags=[:unit, :commons] begin
    using LineCableModels.Commons: vacuum_permittivity

    # Rational arithmetic keeps the decimal definition exact.
    @test vacuum_permittivity(Rational{BigInt}) == 88541878128 // big(10)^22
    for T in (Float32, Float64, BigFloat)
        @test vacuum_permittivity(T) isa T
    end
    # The Float64 evaluation of 88541878128 × 10⁻²² is one ulp below the decimal literal.
    @test vacuum_permittivity(Float64) == prevfloat(8.8541878128e-12)
end
