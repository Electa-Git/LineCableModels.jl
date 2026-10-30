@testset "Engine spectral tools" begin
    Engine = LineCableModels.Engine
    integrate = Engine.integrate
    select = LineCableModels.select
    section = select(integrate => :quad, integrate => :quad)
    @test section == (method = Val(:quad), rtol = 1e-8, atol = 0.0, maxevals = 10^7)
    @test select(integrate => :quad, (method = :quad, rtol = 1e-10)) ==
          (method = Val(:quad), rtol = 1e-10, atol = 0.0, maxevals = 10^7)
    @test_throws ArgumentError select(integrate => :quad, (method = :simpson,))
    @test_throws ArgumentError select(integrate => :quad, (tolerance = 1e-6,))
    untargeted = select(integrate => :quad, (rtol = 0.0,))
    @test_throws ArgumentError LineCableModels.validate(untargeted, integrate, Val(:quad))
    buffers = Engine.initialize_buffers(Engine.SpectralIntegral, Val(:quad), Float64, (;))
    # On the real axis the path is the identity, and the test gives the breakpoints.
    decay = Engine.SpectralIntegral(λ -> exp(-λ), complex(1.0))
    points = buffers.integral.points
    append!(empty!(points), (0.5, 2.0))
    value, estimate = integrate(decay, section.method, section, buffers)
    @test isapprox(value, 1; rtol = 1e-8)
    @test estimate <= 1e-8 * abs(value)
    empty!(points)
    @test_logs (:warn, r"estimated error above the requested target") integrate(decay,
        section.method, select(integrate => :quad, (rtol = 1e-15, maxevals = 10)), buffers)
    push!(points, -1.0)
    @test_throws DomainError integrate(decay, section.method, section, buffers)

    # Along the path, the geometric factor of a buried and of an overhead pair integrates to
    # h / (h² + y²), the value on the real axis. The wavenumbers of air and of the earth are
    # those of 0.1 Ω·m at 1 kHz.
    wavenumbers = (complex(0.0, 2.0958450213811722e-5),
        complex(0.198691764763236, 0.19869176586860801))
    tight = select(integrate => :quad, (rtol = 1e-13,))
    for (height, separation) in ((2.0, 1.0), (11.0, 1.0))
        factor = Engine.SpectralIntegral(
            λ -> (exp(-(height + im * separation) * λ) +
                  exp(-(height - im * separation) * λ)) / 2,
            buffers, wavenumbers, (), (), height, separation)
        @test 0 < angle(factor.ζ) <= π / 6
        @test height * cos(angle(factor.ζ)) - separation * sin(angle(factor.ζ)) > 0
        @test issorted(points) && allunique(points) && first(points) >= 0 && 1.0 in points
        @test isapprox(first(integrate(factor, tight.method, tight, buffers)),
            height / (height^2 + separation^2); rtol = 1e-12)
    end
    # The path passes below a pole in the first quadrant. A zero pole is no pole.
    pole = complex(2.0, 0.5)
    @test angle(Engine.SpectralIntegral(
        identity, buffers, wavenumbers, (pole,), (), 2.0, 1.0).ζ) <=
          angle(pole) / 2
    path = Engine.SpectralIntegral(identity, buffers, wavenumbers, (), (), 2.0, 1.0).ζ
    written = copy(points)
    @test Engine.SpectralIntegral(
        identity, buffers, wavenumbers, (complex(0.0),), (), 2.0, 1.0).ζ == path
    @test points == written
    @test_throws DomainError Engine.SpectralIntegral(
        identity, buffers, wavenumbers, (), (), 0.0, 1.0)
    # After the first calls above, the constructor and the quadrature do not allocate.
    @test (@allocated LineCableModels.SpectralIntegral(
        identity, buffers, wavenumbers, (), (), 2.0, 1.0)) == 0
    @test (@allocated LineCableModels.integrate(decay, section.method, section, buffers)) ==
          0
end
