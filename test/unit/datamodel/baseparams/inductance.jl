@testitem "BaseParams / inductance / tubular and equivalent permeability" tags=[:unit] setup=[
    BaseParamsTestSupport,
    UseBaseParamsSupport
] begin
    μ0=4π*1e-7
    inductance_value=tubular_inductance(0.01, 0.02, 1.0)
    @test inductance_value ≈ μ0/(2π)*log(2)
    @test tubular_inductance(0.01, 0.02, 2.0) ≈ 2inductance_value
    @test tubular_inductance(0.01, 0.02, 0.0) == 0.0
    @test_throws DomainError tubular_inductance(0.0, 0.02, 1.0)
    @test_throws DomainError tubular_inductance(0.02, 0.01, 1.0)
    @test_throws DomainError tubular_inductance(0.01, 0.02, -1.0)

    for T in (Float32, Float64, BigFloat)
        result=tubular_inductance(T(0.01), T(0.02), one(T))
        @test result isa T
        radius=tubular_gmr(T(0.02), T(0.01), T(1.5))
        @test equivalent_mu(radius, T(0.02), T(0.01)) ≈ T(1.5)
    end

    uncertain=tubular_inductance(
        measurement(0.01, 1e-5),
        measurement(0.02, 1e-5),
        measurement(1.0, 1e-3)
    )
    @test uncertain isa Measurement{Float64}
    @test uncertainty(uncertain) > 0

    @test equivalent_mu(tubular_gmr(0.02, 0.0, 1.0), 0.02, 0.0) ≈ 1.0
    @test equivalent_mu(0.02, 0.02, 0.02) ≈ 0.0
    @test_throws DomainError equivalent_mu(-0.01, 0.02, 0.01)
    @test_throws DomainError equivalent_mu(0.015, 0.01, 0.02)
end

@testitem "BaseParams / inductance / trifoil estimate" tags=[:unit] setup=[
    BaseParamsTestSupport,
    UseBaseParamsSupport
] begin
    args=(0.010, 0.015, 1.72e-8, 1.0, 0.020, 0.025, 2.82e-8, 1.0, 0.100)
    reference=trifoil_inductance(args...; rho_e = 100.0, f = 50.0)
    @test reference ≈ 1.573964832699787e-7
    @test trifoil_inductance(args...; rho_e = 100.0, f = 60.0) < reference

    better_screen=Base.setindex(args, args[7]/10, 7)
    magnetic_core=Base.setindex(args, 2args[4], 4)
    @test trifoil_inductance(better_screen...; rho_e = 100.0, f = 50.0) < reference
    @test trifoil_inductance(magnetic_core...; rho_e = 100.0, f = 50.0) > reference

    args32=Float32.(args)
    result32=trifoil_inductance(args32...)
    @test result32 isa Float32
    @test result32 ≈ Float32(reference) rtol=sqrt(eps(Float32))

    uncertain_args=Base.setindex(args, measurement(args[9], 1e-4), 9)
    uncertain=trifoil_inductance(uncertain_args...; rho_e = 100.0, f = 50.0)
    @test uncertain isa Measurement{Float64}
    @test uncertainty(uncertain) > 0
end
