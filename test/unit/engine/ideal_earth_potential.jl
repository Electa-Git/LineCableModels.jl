@testitem "Engine / ideal external potential / images, zero cases and types" tags=[:unit, :engine] begin
    const E = LineCableModels.Engine
    function evaluate(selected, arguments...)
        return @inferred selected(arguments...)
    end
    for T in (Float32, Float64, BigFloat)
        μ0, ε0 = 4T(pi)*T(10)^(-7), parse(T, "8.8541878128e-12")
        s = complex(zero(T), 100T(pi))
        formula = E.EarthAdmittance.Formula(:ideal)
        for (kind, heights, x, layers, reference) in (
            (:self, (10, 10), 0, (1, 1),
            "136627009046.11958786519963276511094116675099612472"),
            (:mutual, (10, 12), 0.75, (1, 1),
            "41930305761.703608584723801139014535060706292407094"),
            (:self, (-1, -1), 0, (2, 2), "0"),
            (:mutual, (-1, -2), 0.75, (2, 2), "0"),
            (:mutual, (10, -2), 0.75, (1, 2), "0"),
            (:mutual, (-2, 10), 0.75, (2, 1), "0"))
            pair = E.EarthPair(1, kind === :self ? 1 : 2, T.(heights), T(x), layers;
                radius = kind === :self ? T(1)/100 : nothing)
            value = evaluate(formula, T[Inf, 100], [ε0, 10ε0], [μ0, μ0], s, pair)
            @test value isa Complex{T}
            @test iszero(imag(value))
            if reference == "0"
                @test iszero(value)
            elseif T !== Float32
                @test real(value) ≈ parse(T, reference) rtol=5e-14
            end
        end
    end
end

@testitem "Engine / ideal potential / geometric uncertainty" tags=[:unit, :measurements] begin
    using Measurements
    const E=LineCableModels.Engine
    ε0, μ0=8.8541878128e-12, 4pi*1e-7
    function potential(radius)
        T=typeof(radius)
        pair=E.EarthPair(1, 1, (T(10), T(10)), T(0), (1, 1); radius)
        return @inferred E.EarthAdmittance.Formula(:ideal)(T[Inf, 100], T[ε0, 10ε0],
            T[μ0, μ0], Complex{T}(100pi*im), pair)
    end
    radius=measurement(0.01, 0.0001)
    coefficient=potential(radius)
    @test Measurements.derivative(real(coefficient), radius) ≈ -1/(2pi*ε0*0.01)
    @test iszero(uncertainty(real(coefficient-potential(radius))))
    @test uncertainty(real(coefficient-potential(measurement(0.01, 0.0001)))) > 0
    @test iszero(uncertainty(imag(coefficient)))
end
