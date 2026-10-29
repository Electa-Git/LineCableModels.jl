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
    decay = Engine.SpectralIntegral(t -> complex(exp(-t)))
    value,
    estimate = integrate(decay, section.method, section, buffers; points = [0.5, 2.0])
    @test isapprox(value, 1; rtol = 1e-8)
    @test estimate <= 1e-8 * abs(value)
    @test_logs (:warn, r"estimated error above the requested target") integrate(decay,
        section.method, select(integrate => :quad, (rtol = 1e-15, maxevals = 10)), buffers;
        points = Float64[])
    @test_throws DomainError integrate(decay, section.method, section, buffers;
        points = [-1.0])

    # Along the ray, the geometric factor of a buried and of an overhead pair integrates to
    # h / (h² + y²), the value on the real axis. The wavenumbers of air and of the earth are
    # those of 0.1 Ω·m at 1 kHz.
    wavenumbers = (complex(0.0, 2.0958450213811722e-5),
        complex(0.198691764763236, 0.19869176586860801))
    tight = select(integrate => :quad, (rtol = 1e-13,))
    for (height, separation) in ((2.0, 1.0), (11.0, 1.0))
        ray = Engine.ray!(buffers, wavenumbers, (), (), height, separation)
        @test 0 < angle(ray) <= π / 6
        @test height * cos(angle(ray)) - separation * sin(angle(ray)) > 0
        points = buffers.integral.points
        @test issorted(points) && allunique(points) && first(points) >= 0 && 1.0 in points
        factor = Engine.SpectralIntegral(t -> ray *
                                              (exp(-(height + im * separation) * ray * t) +
                                               exp(-(height - im * separation) * ray * t)) /
                                              2)
        @test isapprox(first(integrate(factor, tight.method, tight, buffers; points)),
            height / (height^2 + separation^2); rtol = 1e-12)
    end
    # The ray passes below a pole in the first quadrant. A zero pole is no pole.
    pole = complex(2.0, 0.5)
    @test angle(Engine.ray!(buffers, wavenumbers, (pole,), (), 2.0, 1.0)) <= angle(pole) / 2
    ray = Engine.ray!(buffers, wavenumbers, (), (), 2.0, 1.0)
    points = copy(buffers.integral.points)
    @test Engine.ray!(buffers, wavenumbers, (complex(0.0),), (), 2.0, 1.0) == ray
    @test buffers.integral.points == points
    @test_throws DomainError Engine.ray!(buffers, wavenumbers, (), (), 0.0, 1.0)
    # After the first calls above, the ray and the quadrature do not allocate.
    @test (@allocated LineCableModels.ray!(buffers, wavenumbers, (), (), 2.0, 1.0)) == 0
    @test (@allocated LineCableModels.integrate(decay, section.method, section, buffers;
        points = buffers.integral.points)) == 0
end
