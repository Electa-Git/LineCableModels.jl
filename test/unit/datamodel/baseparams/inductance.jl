@testitem "BaseParams / inductance / GMR and equivalent permeability" tags=[:unit] setup=[
    BaseParamsTestSupport,
    UseBaseParamsSupport
] begin
    for T in (Float32, Float64, BigFloat)
        radius=tubular_gmr(T(0.02), T(0.01), T(1.5))
        @test equivalent_mu(radius, T(0.02), T(0.01)) ≈ T(1.5)
    end

    @test equivalent_mu(tubular_gmr(0.02, 0.0, 1.0), 0.02, 0.0) ≈ 1.0
    @test equivalent_mu(0.02, 0.02, 0.02) ≈ 0.0
    @test_throws DomainError equivalent_mu(-0.01, 0.02, 0.01)
    @test_throws DomainError equivalent_mu(0.015, 0.01, 0.02)

    lay_radius=0.01
    wire_radius=0.001
    wire_count=12
    coordinates=[(
                     lay_radius*cos(2π*index/wire_count),
                     lay_radius*sin(2π*index/wire_count)
                 )
                 for index in 0:(wire_count - 1)]
    @test strand_gmr(coordinates, wire_radius, 1.0) ≈
          strand_gmr(lay_radius, wire_count, wire_radius, 1.0)
    @test_throws ArgumentError strand_gmr(Tuple{Float64, Float64}[], wire_radius, 1.0)
    @test_throws ArgumentError strand_gmr([(0.0, 0.0), (0.0, 0.0)], wire_radius, 1.0)

    first_gmr=0.004
    second_gmr=0.012
    distance=0.009
    combined=equivalent_gmr(first_gmr, 2.0, second_gmr, 2.0, distance)
    @test combined ≈ first_gmr^0.25*second_gmr^0.25*sqrt(distance)
    @test equivalent_gmr(
        Float32(first_gmr),
        Float32(2),
        Float32(second_gmr),
        Float32(2),
        Float32(distance)
    ) isa Float32
    @test_throws DomainError equivalent_gmr(0.0, 1.0, second_gmr, 1.0, distance)
    @test_throws DomainError equivalent_gmr(first_gmr, 0.0, second_gmr, 1.0, distance)
    @test_throws DomainError equivalent_gmr(first_gmr, 1.0, second_gmr, 1.0, 0.0)
end
