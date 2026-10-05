@testitem "Commons / nominal and uncertainty / deterministic, complex and array values" tags=[:unit, :commons] begin
    # A minimal uncertain number with its own scalar nominal value and uncertainty.
    struct UncertainProbe <: Real
        center::Float64
        sigma::Float64
    end
    LineCableModels.nominal(value::UncertainProbe) = value.center
    LineCableModels.uncertainty(value::UncertainProbe) = value.sigma

    for value in (3, 2.5f0, 1.5, 4.0 + 5.0im)
        @test nominal(value) === value
    end
    @test nominal(complex(UncertainProbe(1.0, 0.1), UncertainProbe(2.0, 0.2))) ===
          complex(1.0, 2.0)
    values = [UncertainProbe(1.0, 0.1), UncertainProbe(2.0, 0.2)]
    nominals = nominal(values)
    @test nominals isa Vector{Float64}
    @test nominals == [1.0, 2.0]
    @test values == [UncertainProbe(1.0, 0.1), UncertainProbe(2.0, 0.2)]

    @test uncertainty(3) === 0
    @test uncertainty(2.5f0) === 0.0f0
    @test uncertainty(1.5) === 0.0
    @test uncertainty(complex(UncertainProbe(1.0, 0.1), UncertainProbe(2.0, 0.2))) ===
          complex(0.1, 0.2)
    sigmas = uncertainty(values)
    @test sigmas isa Vector{Float64}
    @test sigmas == [0.1, 0.2]
    @test_throws MethodError uncertainty("not numeric")
end
