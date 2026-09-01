@testitem "BaseParams / numeric matrix / precision and exactness" tags=[:unit] setup=[
    BaseParamsTestSupport,
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
            capacitance=shunt_capacitance(r_in, r_ex, T(2.3))
            conductance=shunt_conductance(r_in, r_ex, T(1e9))

            @test resistance isa T
            @test capacitance isa T
            @test conductance isa T

            @test TestNumerics.isapprox_scaled(
                resistance,
                rho / (T(π) * (r_ex^2 - r_in^2))
            )
            @test TestNumerics.isapprox_scaled(
                conductance,
                T(2) * T(π) / (T(1e9) * log(r_ex / r_in))
            )
            @test capacitance > zero(T)
        end
    end

    @test parallel(1 // 2, 1 // 3) == 1 // 5
    complex_result=parallel(ComplexF32(1, 2), ComplexF32(3, 4))
    @test complex_result isa ComplexF32

    normalized=shunt_capacitance(Float32(0.01), 0.02, big"2.3")
    @test normalized isa BigFloat
end

@testitem "BaseParams / numeric matrix / cross-precision convergence" tags=[:unit] setup=[
    BaseParamsTestSupport,
    UseBaseParamsSupport
] begin
    using TOML

    reference=TOML.parsefile(joinpath(
        pkgdir(LineCableModels),
        "test",
        "fixtures",
        "reference",
        "coaxial_capacitance.toml"
    ))
    oracle=setprecision(BigFloat, reference["precision_bits"]) do
        parse(BigFloat, reference["value"])
    end

    approximations=(
        shunt_capacitance(Float32(0.01), Float32(0.02), Float32(2.3)),
        shunt_capacitance(0.01, 0.02, 2.3),
        setprecision(BigFloat, 128) do
            shunt_capacitance(big"0.01", big"0.02", big"2.3")
        end
    )
    errors=map(value->abs(BigFloat(value)-oracle)/abs(oracle), approximations)
    @test errors[2] <= errors[1]
    @test errors[3] <= errors[2]
    @test errors[1] < sqrt(eps(Float32))
    @test errors[2] < sqrt(eps(Float64))
end
