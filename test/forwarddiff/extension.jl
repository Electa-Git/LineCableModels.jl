using LineCableModels: nominal
using ForwardDiff: Dual, partials
using Measurements: Measurements, measurement
using SpecialFunctions: besseli, besselj, besselix, besselkx, besseljx, besselk

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

    @testset "the partial of a conductor's outer surface at a harmonic order in its radius" begin
        # Measured against a central difference of step 1e-7 r_ex: 3.3e-8. The value lies
        # within 5.7e-16 of the record of the nominal data.
        # A steel pipe, a solid of aluminum and a copper sheath: radii in m, resistivity in Ω·m and
        # relative permeability.
        formula = LineCableModels.Engine.InternalImpedance.Formula(Val(:schelkunoff1934))
        cases = ((62.0e-3, 65.0e-3, 7.0e-8, 10.0), (0.0, 19.5e-3, 3.365e-8, 1.0),
            (20.5e-3, 21.5e-3, 1.7241e-8, 1.0))
        for (r_in, r_ex, ρ, μ_r) in cases, n in (1, 3, 20), f in (50.0, 1.0e5)
            s = im * (2π * f)
            surface = formula(n, ForwardDiff.Dual(r_in, 0.0), ForwardDiff.Dual(r_ex, 1.0),
                ForwardDiff.Dual(ρ, 0.0), ForwardDiff.Dual(μ_r, 0.0),
                complex(ForwardDiff.Dual(0.0, 0.0), ForwardDiff.Dual(2π * f, 0.0)))
            h = 1.0e-7 * r_ex
            difference = (formula(n, r_in, r_ex + h, ρ, μ_r, s).outer -
                          formula(n, r_in, r_ex - h, ρ, μ_r, s).outer) / 2h
            partial = complex(ForwardDiff.partials(real(surface.outer))[1],
                ForwardDiff.partials(imag(surface.outer))[1])
            @test isapprox(partial, difference; rtol = 1.0e-6)
            @test isapprox(nominal(surface.outer), formula(n, r_in, r_ex, ρ, μ_r, s).outer;
                rtol = 1.0e-14)
        end
    end
end

# A scaled Bessel function times the factor that removes its scaling is the unscaled
# function. That function has a derivative on the axis where the scaling has a cusp. Each
# partial derivative of the product is a sum of at most four products of numbers below 1.
# So it differs from the closed form by less than 16eps().
@testset "a scaled Bessel function times its factor keeps the derivative of the unscaled one" begin
    cases = (
        (z -> besselix(0, z) * exp(abs(real(z))), z -> besseli(1, z), (0.0, 1.2)),
        (z -> besselix(0, z) * exp(abs(real(z))), z -> besseli(1, z), (-0.0, 1.2)),
        (z -> besselix(0, z) * exp(abs(real(z))), z -> besseli(1, z), (0.3, 1.2)),
        (z -> besselix(1, z) * exp(abs(real(z))),
            z -> (besseli(0, z) + besseli(2, z)) / 2, (0.0, 1.2)),
        (z -> besselix(1, z) * exp(abs(real(z))),
            z -> (besseli(0, z) + besseli(2, z)) / 2, (-0.0, 1.2)),
        (z -> besseljx(0, z) * exp(abs(imag(z))), z -> -besselj(1, z), (1.2, 0.0)),
        (z -> besseljx(0, z) * exp(abs(imag(z))), z -> -besselj(1, z), (1.2, -0.0)),
        (z -> besseljx(0, z) * exp(abs(imag(z))), z -> -besselj(1, z), (1.2, -0.3)))
    @testset "at $(complex(x, y))" for (product, derivative, (x, y)) in cases
        exact = derivative(complex(x, y))
        a, b = measurement(x, 0.1), measurement(y, 0.1)
        w = product(complex(a, b))
        @test complex(Measurements.derivative(real(w), a), Measurements.derivative(imag(w), a)) ≈
              exact atol = 16eps()
        @test complex(Measurements.derivative(real(w), b), Measurements.derivative(imag(w), b)) ≈
              im * exact atol = 16eps()
        w = product(complex(Dual(x, 1.0, 0.0), Dual(y, 0.0, 1.0)))
        @test complex(partials(real(w), 1), partials(imag(w), 1)) ≈ exact atol = 16eps()
        @test complex(partials(real(w), 2), partials(imag(w), 2)) ≈ im * exact atol = 16eps()
    end
end
