@testitem "Engine / internal impedance / solid and hollow working precision" tags=[:unit] begin
    const II = LineCableModels.Engine.InternalImpedance
    formula = II.Formula(:default)
    @test II.formula_id(formula) === :Schelkunoff1934

    # Float64 uses SpecialFunctions; Complex{BigFloat} uses the package's
    # precision-preserving Bessel implementation. Compare the complete surface
    # impedances rather than testing an otherwise unused numerical helper.
    setprecision(BigFloat, 128) do
        for inner in (0.0, 0.008), frequency in (50.0, 1000.0)
            reference = formula(inner, 0.01, 1.7241e-8, 1.0, complex(0.0, 2pi * frequency))
            precise = @inferred formula(BigFloat(inner), BigFloat(0.01),
                BigFloat(1.7241e-8), BigFloat(1.0), complex(big"0", 2big(pi) * frequency))
            for interaction in (Val(:inner), Val(:outer), Val(:mutual))
                actual = @inferred precise(interaction)
                @test actual isa Complex{BigFloat}
                @test precision(real(actual)) == 128
                @test precision(imag(actual)) == 128
                @test isfinite(actual)
                @test actual ≈ reference(interaction) rtol=2e-12 atol=0
            end
            if iszero(inner)
                @test iszero(precise(Val(:inner)))
                @test iszero(precise(Val(:mutual)))
            end
        end

        # Independent DC and low-frequency internal-inductance limits of a
        # solid round wire; the test frequency is explicitly positive.
        radius = big"0.01"
        resistivity = big"1.7241e-8"
        omega = 2big(pi) * big"1e-6"
        impedance = formula(BigFloat(0), radius, resistivity, BigFloat(1),
            complex(big"0", omega))(Val(:outer))
        @test real(impedance) ≈ resistivity / (big(pi) * radius^2) rtol=1e-12
        @test imag(impedance) / omega ≈ big"5e-8" rtol=1e-12
    end
end
