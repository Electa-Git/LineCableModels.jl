using LineCableModels: nominal
using SpecialFunctions: besselix, besselkx, besseljx, besselk

@testset "ForwardDiff extension" begin
    h = 1.0e-6
    both = (0.7 + 0.4im, -0.7 + 0.4im, 0.3 - 1.1im, 2.5 + 2.5im)
    right = (0.7 + 0.4im, 0.3 - 1.1im, 2.5 + 2.5im)
    @testset "the partials of $f against central differences" for (f, points) in (
        (besselix, both), (besseljx, both), (besselkx, right), (besselk, right))
        for ν in (0, 1), z in points

            dx = ForwardDiff.derivative(t -> f(ν, z + t), 0.0)
            dy = ForwardDiff.derivative(t -> f(ν, z + im * t), 0.0)
            @test dx ≈ (f(ν, z + h) - f(ν, z - h)) / 2h rtol = 1.0e-6
            @test dy ≈ (f(ν, z + im * h) - f(ν, z - im * h)) / 2h rtol = 1.0e-6
            @test nominal(f(ν, complex(ForwardDiff.Dual(real(z), 1.0), imag(z)))) == f(ν, z)
        end
    end

    @testset "the second derivative of $f through a nested dual" for f in (
        besselix, besseljx, besselkx, besselk)
        d2 = ForwardDiff.derivative(
            t -> ForwardDiff.derivative(s -> f(0, complex(s, 0.3s)), t), 0.8)
        d1 = (ForwardDiff.derivative(s -> f(0, complex(s, 0.3s)), 0.8 + h) -
              ForwardDiff.derivative(s -> f(0, complex(s, 0.3s)), 0.8 - h)) / 2h
        @test d2 ≈ d1 rtol = 1.0e-6
    end

    @testset "the nominal value of a dual" begin
        nested = ForwardDiff.Dual(ForwardDiff.Dual(0.8, 1.0), ForwardDiff.Dual(1.0, 0.0))
        @test nominal(nested) === 0.8
        @test nominal(complex(ForwardDiff.Dual(0.8, 1.0), ForwardDiff.Dual(0.6, 0.0))) ===
              complex(0.8, 0.6)
    end
end
