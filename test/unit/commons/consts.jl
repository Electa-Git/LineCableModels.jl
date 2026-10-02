@testitem "Commons / vacuum constants are evaluated in the requested type" tags=[:unit] begin
    using LineCableModels.Commons: vacuum_permittivity, vacuum_permeability

    # Rational arithmetic keeps the decimal definition exact.
    @test vacuum_permittivity(Rational{BigInt}) == 88541878128 // big(10)^22
    for T in (Float32, Float64, BigFloat)
        @test vacuum_permittivity(T) isa T
        @test vacuum_permeability(T) isa T
    end
    @test vacuum_permeability(BigFloat) ≈ 4 * big(π) / big(10)^7 rtol = 4eps(BigFloat)
    # The Float64 evaluation of 88541878128 × 10⁻²² is one ulp below the decimal literal.
    @test vacuum_permittivity(Float64) == prevfloat(8.8541878128e-12)
    @test vacuum_permeability(Float64) == 4π * 1e-7
    @test inv(sqrt(vacuum_permittivity(BigFloat) * vacuum_permeability(BigFloat))) ≈
          299792458 rtol = 1e-9
end
