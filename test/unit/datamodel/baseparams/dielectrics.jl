@testitem "BaseParams / dielectrics / coaxial kernels" tags=[:unit] setup=[
    BaseParamsTestSupport,
    UseBaseParamsSupport
] begin
    ε0=8.8541878128e-12
    capacitance=shunt_capacitance(0.01, 0.02, 2.3)
    conductance=shunt_conductance(0.01, 0.02, 1e9)
    @test capacitance ≈ 2π*ε0*2.3/log(2)
    @test conductance ≈ 2π/(1e9*log(2))
    @test equivalent_eps(capacitance, 0.02, 0.01) ≈ 2.3
    @test equivalent_conductivity(conductance, 0.01, 0.02) ≈ 1e-9
    @test shunt_capacitance(0.01, 0.02, 3.0) > capacitance
    @test shunt_conductance(0.01, 0.02, 1e8) > conductance

    for call in (
        ()->shunt_capacitance(0.0, 0.02, 2.3),
        ()->shunt_capacitance(0.02, 0.01, 2.3),
        ()->shunt_capacitance(0.01, 0.02, -1.0),
        ()->shunt_conductance(0.01, 0.02, 0.0),
        ()->equivalent_eps(-1e-10, 0.02, 0.01),
        ()->equivalent_conductivity(-1e-9, 0.01, 0.02)
    )
        @test_throws DomainError call()
    end

    for T in (Float32, Float64, BigFloat)
        C=shunt_capacitance(T(0.01), T(0.02), T(2.3))
        G=shunt_conductance(T(0.01), T(0.02), T(1e9))
        @test C isa T
        @test G isa T
        @test equivalent_eps(C, T(0.02), T(0.01)) isa T
        @test equivalent_conductivity(G, T(0.01), T(0.02)) isa T
    end

    uncertain=shunt_capacitance(
        measurement(0.01, 1e-4),
        measurement(0.02, 1e-4),
        measurement(2.3, 1e-2)
    )
    @test uncertain isa Measurement{Float64}
    @test uncertainty(uncertain) > 0

    omega=2π*50
    first_admittance=complex(2e-9, omega*3e-10)
    second_admittance=complex(7e-9, omega*8e-10)
    expected=inv(inv(first_admittance)+inv(second_admittance))
    combined=series_shunt_admittance(2e-9, 3e-10, 7e-9, 8e-10, omega)
    @test combined.conductance ≈ real(expected)
    @test combined.capacitance ≈ imag(expected)/omega
    @test series_shunt_admittance(
        Float32(2e-9),
        Float32(3e-10),
        Float32(7e-9),
        Float32(8e-10),
        Float32(omega)
    ).conductance isa Float32
    @test_throws DomainError series_shunt_admittance(
        2e-9,
        3e-10,
        7e-9,
        8e-10,
        0.0
    )
end
