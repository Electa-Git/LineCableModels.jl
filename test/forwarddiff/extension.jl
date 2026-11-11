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
