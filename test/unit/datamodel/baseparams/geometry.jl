@testitem "BaseParams / geometry / ring and helical invariants" tags=[:unit] begin
    const BP = LineCableModels.DataModel.BaseParams
    using LinearAlgebra: norm

    for T in (Float32, Float64, BigFloat)
        setprecision(BigFloat, 128) do
            radius = T(0.5e-3)
            inner, outer, ratio = T(0.01), T(0.02), T(12)
            mean_diameter, pitch_length, overlength = @inferred BP.helix(inner, outer, ratio)
            @test (mean_diameter, pitch_length, overlength) isa NTuple{3,T}
            @test mean_diameter == inner + outer
            @test pitch_length == ratio * mean_diameter
            @test overlength ≈ sqrt(one(T) + (T(pi) / ratio)^2)
            @test all(isapprox.(BP.helix(2inner, 2outer, ratio),
                (2mean_diameter, 2pitch_length, overlength)))
            @test BP.helix(inner, outer, zero(T)) == (mean_diameter, zero(T), one(T))

            centre = (T(0.02), T(-0.03))
            positions = BP.wire_coordinates(6, radius, radius; C=centre)
            @test eltype(positions) === Tuple{T,T}
            @test length(positions) == 6
            @test all(position -> hypot(position[1] - centre[1], position[2] - centre[2]) ≈ 2radius,
                positions)
            for index in eachindex(positions)
                neighbour = positions[mod1(index + 1, length(positions))]
                @test norm(collect(positions[index]) - collect(neighbour)) ≈ 2radius
            end
            @test BP.wire_coordinates(6, radius, radius, centre) == positions
            @test BP.wire_coordinates(1, radius, inner; C=centre) == [centre]
            @test isempty(BP.wire_coordinates(0, radius, inner))
            @test BP.wire_coordinates(6, radius, radius)[1] == (2radius, zero(T))

            turns = T(4)
            correction = BP.solenoid_factor(turns, inner, outer)
            @test correction isa T
            @test correction > one(T)
            @test BP.solenoid_factor(zero(T), inner, outer) == one(T)
            @test BP.solenoid_factor(2turns, inner, outer) - 1 ≈ 4(correction - 1)
            @test BP.solenoid_factor(turns / 2, 2inner, 2outer) ≈ correction
        end
    end
    @test BP.helix(Float32(0.01), big"0.02", 12) isa NTuple{3,BigFloat}
    @test_throws DomainError BP.helix(-0.01, 0.02, 12)
    @test_throws DomainError BP.helix(0.02, 0.01, 12)
    @test_throws DomainError BP.helix(0.01, 0.02, -12)
    @test_throws DomainError BP.helix(0.01, Inf, 12)
    @test_throws DomainError BP.wire_coordinates(-1, 0.001, 0.01)
    @test_throws DomainError BP.wire_coordinates(6, 0.0, 0.01)
    @test_throws DomainError BP.wire_coordinates(6, 0.001, -0.01)
    @test_throws DomainError BP.wire_coordinates(6, 0.001, 0.01; C=(NaN, 0.0))
    @test_throws DomainError BP.solenoid_factor(-4.0, 0.01, 0.02)
    @test_throws DomainError BP.solenoid_factor(4.0, -0.01, 0.02)
    @test_throws DomainError BP.solenoid_factor(4.0, 0.02, 0.01)
end
