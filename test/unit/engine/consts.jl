@testitem "Engine / vacuum permeability is evaluated in the requested type" tags=[:unit, :engine] begin
    using LineCableModels.Commons: vacuum_permittivity
    using LineCableModels.Engine: vacuum_permeability

    for T in (Float32, Float64, BigFloat)
        @test vacuum_permeability(T) isa T
    end
    @test vacuum_permeability(BigFloat) ≈ 4 * big(π) / big(10)^7 rtol = 4eps(BigFloat)
    @test vacuum_permeability(Float64) == 4π * 1e-7
    @test inv(sqrt(vacuum_permittivity(BigFloat) * vacuum_permeability(BigFloat))) ≈
          299792458 rtol = 1e-9
end
