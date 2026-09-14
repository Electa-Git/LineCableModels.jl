@testitem "Comparison / indexed observations retain numerical-zero unavailability" tags=[:unit] begin
    using LineCableModels.Engine: absolute_error, relative_error
    a = reshape(ComplexF64[1, 2, 3, 4], 2, 2, 1)
    z = zeros(ComplexF64, 2, 2, 1)
    reference = LineParameters(PhaseDomain, a, z, [50.0])
    candidate = LineParameters(PhaseDomain, 2a, z, [50.0])
    result = LineCableModels.Engine.compare(reference, candidate)
    @test basis(result) === :pul
    @test observe(result, Z, absolute_error, 2, 1) == 2.0
    @test observe(result, Z, relative_error, 2, 1) == 1.0
    @test observe(result, Y, absolute_error, 2, 1) == 0.0
    @test ismissing(observe(result, Y, relative_error, 2, 1))
end
