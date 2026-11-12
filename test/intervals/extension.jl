using LineCableModels: nominal
using IntervalArithmetic: interval, inf, sup, mid
using Arblib: Arblib, Arb, Acb
using SpecialFunctions: besselix, besselkx, besseljx, besselk, besseli, besselj

# Each function's definition in Arb, the reference of an enclosure at 256 bits.
const ARB_DEFINITIONS = (
    (besselix, (ν, w) -> besseli(Acb(ν), w) * exp(-abs(real(w)))),
    (besseljx, (ν, w) -> besselj(Acb(ν), w) * exp(-abs(imag(w)))),
    (besselkx, (ν, w) -> besselkx(Acb(ν), w)),
    (besselk, (ν, w) -> besselk(Acb(ν), w)))

# The enclosure e contains the ball of the definition at the point z, at 256 bits.
function encloses(e::Complex, definition, ν, z::Complex)
    value = definition(ν, Acb(z; prec = 256))
    return Arblib.contains(Arb(real(e); prec = 256), real(value)) &&
           Arblib.contains(Arb(imag(e); prec = 256), imag(value))
end

function interval_box(z, δ)
    return complex(interval(real(z) - δ, real(z) + δ), interval(imag(z) - δ, imag(z) + δ))
end

function enclosure_width(e::Complex)
    return max(sup(real(e)) - inf(real(e)), sup(imag(e)) - inf(imag(e)))
end

@testset "IntervalArithmetic extension" begin
    both = (0.7 + 0.4im, -0.7 + 0.4im, 0.3 - 1.1im, 2.5 + 2.5im)
    right = (0.7 + 0.4im, 0.3 - 1.1im, 2.5 + 2.5im)
    @testset "the enclosures of $f" for ((f, definition), points) in zip(ARB_DEFINITIONS,
        (both, both, right, right))
        for ν in (0, 1), z in points

            thin = interval_box(z, 1.0e-10)
            center = complex(mid(real(thin)), mid(imag(thin)))
            e = f(ν, thin)
            @test complex(mid(real(e)), mid(imag(e))) ≈ f(ν, center) rtol = 1.0e-12
            @test encloses(e, definition, ν, center)
            wide = f(ν, interval_box(z, 0.1))
            @test all(encloses(wide, definition, ν, complex(a, b))
            for a in range(real(z) - 0.1, real(z) + 0.1; length = 8),
            b in range(imag(z) - 0.1, imag(z) + 0.1; length = 8))
            point = f(ν, complex(interval(real(z)), interval(imag(z))))
            @test enclosure_width(point) < 1.0e-12 * abs(f(ν, z))
        end
    end

    @testset "the enclosures of besseljx on a real interval" begin
        for ν in (0, 1), x in (0.7, 2.0, 12.5)

            thin = interval(x - 1.0e-10, x + 1.0e-10)
            e = besseljx(ν, thin)
            @test mid(e) ≈ besseljx(ν, mid(thin)) rtol = 1.0e-12
            value = besselj(Arb(ν), Arb(mid(thin); prec = 256))
            @test Arblib.contains(Arb(e; prec = 256), value)
            wide = Arb(besseljx(ν, interval(x - 0.1, x + 0.1)); prec = 256)
            @test all(Arblib.contains(wide, besselj(Arb(ν), Arb(a; prec = 256)))
            for a in range(x - 0.1, x + 0.1; length = 64))
            point = besseljx(ν, interval(x))
            @test sup(point) - inf(point) < 1.0e-12 * abs(besseljx(ν, x))
        end
    end

    @testset "the distance of an interval point about the origin encloses its range" begin
        d = distance((interval(-1.0e-3, 1.0e-3), interval(0.0, 0.0)))
        @test inf(d) <= 0 && sup(d) >= 1.0e-3
    end

    @test nominal(interval(1, 2)) == 1.5
end
