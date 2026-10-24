@testitem "BaseParams / numeric matrix / precision and exactness" tags=[:unit, :datamodel] setup=[
    UseBaseParamsSupport,
    TestNumerics
] begin
    for T in (Float32, Float64, BigFloat)
        setprecision(BigFloat, 128) do
            r_in=T(0.01)
            r_ex=T(0.02)
            rho=T(1.7241e-8)

            resistance=tubular_resistance(
                r_in,
                r_ex,
                rho
            )
            permittivity=equivalent_eps(T(1e-10), r_ex, r_in)
            conductivity=equivalent_conductivity(T(1e-9), r_in, r_ex)

            @test resistance isa T
            @test permittivity isa T
            @test conductivity isa T

            @test TestNumerics.isapprox_scaled(
                resistance,
                rho / (T(π) * (r_ex^2 - r_in^2))
            )
            @test TestNumerics.isapprox_scaled(
                conductivity,
                T(1e-9) * log(r_ex / r_in) / (T(2) * T(π))
            )
            @test permittivity > zero(T)
        end
    end

    @test parallel(1 // 2, 1 // 3) == 1 // 5
    complex_result=parallel(ComplexF32(1, 2), ComplexF32(3, 4))
    @test complex_result isa ComplexF32

    normalized=equivalent_eps(big"1e-10", 0.02, Float32(0.01))
    @test normalized isa BigFloat
end

@testitem "BaseParams / numeric matrix / cross-precision convergence" tags=[:unit, :datamodel] setup=[
    UseBaseParamsSupport
] begin
    # Independent electrostatic Dirichlet solution v=log(b/r)/log(b/a): integrate normal
    # displacement on the inner circular boundary. The capacitance of the interval with
    # relative permittivity 3 gives back 3.
    capacitance=setprecision(BigFloat,256) do
        2big(pi)*big"8.8541878128e-12"*3/log(big".01"/big".005")
    end

    approximations=(
        equivalent_eps(Float32(capacitance), Float32(.01), Float32(.005)),
        equivalent_eps(Float64(capacitance), .01, .005),
        setprecision(BigFloat, 128) do
            equivalent_eps(BigFloat(capacitance), big".01", big".005")
        end
    )
    errors=map(value->abs(BigFloat(value)-3)/3, approximations)
    @test errors[2] <= errors[1]
    @test errors[3] <= errors[2]
    @test errors[1] < sqrt(eps(Float32))
    @test errors[2] < sqrt(eps(Float64))
end
