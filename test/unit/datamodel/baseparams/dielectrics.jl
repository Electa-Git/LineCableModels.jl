@testitem "BaseParams / dielectrics / coaxial kernels" tags=[:unit, :measurements] setup=[
    UseBaseParamsSupport
] begin
    # Closed-form coaxial capacitance and conductance of a 10-20 mm interval with
    # relative permittivity 2.3 and resistivity 1e9 ohm-m.
    ε0=8.8541878128e-12
    capacitance=2π*ε0*2.3/log(2)
    conductance=2π/(1e9*log(2))
    @test equivalent_eps(capacitance, 0.02, 0.01) ≈ 2.3
    @test equivalent_conductivity(conductance, 0.01, 0.02) ≈ 1e-9
    @test equivalent_eps(capacitance, 0.03, 0.01) > 2.3
    @test equivalent_conductivity(conductance, 0.01, 0.03) > 1e-9

    for call in (
        ()->equivalent_eps(-1e-10, 0.02, 0.01),
        ()->equivalent_eps(1e-10, 0.01, 0.02),
        ()->equivalent_eps(1e-10, 0.02, 0.0),
        ()->equivalent_conductivity(-1e-9, 0.01, 0.02),
        ()->equivalent_conductivity(1e-9, 0.02, 0.01)
    )
        @test_throws DomainError call()
    end

    for T in (Float32, Float64, BigFloat)
        @test equivalent_eps(T(capacitance), T(0.02), T(0.01)) isa T
        @test equivalent_conductivity(T(conductance), T(0.01), T(0.02)) isa T
    end

    uncertain=equivalent_eps(
        measurement(capacitance, 1e-13),
        measurement(0.02, 1e-4),
        measurement(0.01, 1e-4)
    )
    @test uncertain isa Measurement{Float64}
    @test uncertainty(uncertain) > 0
end
