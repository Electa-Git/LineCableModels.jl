using LineCableModels: nominal
using SpecialFunctions: besselix, besselkx, besseljx, besselk

@testset "Measurements extension" begin
    @test Base.get_extension(LineCableModels, :LineCableModelsMeasurementsExt) isa Module
    @test nominal(measurement(1.0, 0.1)) === 1.0
    @test nominal(complex(1.0 ± 0.1, 2.0 ± 0.2)) === complex(1.0, 2.0)

    z = 0.8 + 0.6im
    σ = 0.01
    h = 1e-6
    for ν in (0, 1)
        I = besselix(ν, z)
        K = besselkx(ν, z)
        J = besseljx(ν, z)
        dI = ν == 0 ? besselix(1, z) : besselix(0, z) - besselix(1, z) / z
        dK = ν == 0 ? -besselkx(1, z) : -besselkx(0, z) - besselkx(1, z) / z
        dJ = ν == 0 ? -besseljx(1, z) : besseljx(0, z) - besseljx(1, z) / z
        dK_plain = ν == 0 ? -besselk(1, z) : -besselk(0, z) - besselk(1, z) / z
        for (f, dx, dy) in ((besselix, dI - sign(real(z)) * I, im * dI),
            (besselkx, K + dK, im * (K + dK)),
            (besseljx, dJ, im * dJ - sign(imag(z)) * J),
            (besselk, dK_plain, im * dK_plain))
            for (argument, direction, derivative) in
                ((complex(measurement(0.8, σ), 0.6), 1.0, dx),
                (complex(0.8, measurement(0.6, σ)), im, dy))
                result = f(ν, argument)
                @test nominal(result) === f(ν, z)
                difference = (f(ν, z + h * direction) - f(ν, z - h * direction)) / (2h)
                for (part, expected, finite) in
                    ((real(result), abs(real(derivative)) * σ, abs(real(difference)) * σ),
                    (imag(result), abs(imag(derivative)) * σ, abs(imag(difference)) * σ))
                    actual = Measurements.uncertainty(part)
                    @test actual ≈ expected rtol=1e-12
                    @test actual ≈ finite rtol=1e-7
                    @test length(Measurements.uncertainty_components(part)) == 1
                end
            end
        end
    end

    for (copper, expected) in
        (
        (Material(:conductor, 1.7241e-8 ± 1.7241e-10, 1, 1, 20, 0.00393),
            Measurement{Float64}),
        (material(kind = :conductor, rho = 1.7241f-8), Float64))
        design = build(CableDesign, "bare_wire",
            Terminal(:core, Region(:core_metal, Circle(0.0425), copper)))
        @test eltype(design) === expected
    end
    wire = build(CableDesign, "float32",
        Terminal(:core,
            Region(:core_metal, Circle(1.5f-3),
                material(kind = :conductor, rho = 1.7241f-8))))
    @test eltype(wire) === Float32

    reference = last(first(REFERENCES))
    copper = Material(:conductor, 1.7241e-8, 1, 0.999994, 20, 0.00393)
    wire = build(CableDesign, "bare_wire",
        Terminal(:core, Region(:core_metal, Circle(0.0425 ± 1e-4), copper)))
    system = build(CableSystem, [wire, wire],
        [Placement(0.0, -1.0), Placement(1.0, -1.0)];
        connections = [(core = 1,), (core = 2,)], id = "two_bare_wires", line_length = 1.0)
    frequencies = reference.frequencies[1:2]
    problem = LineParametersProblem(system, homogeneous(rho = 0.1 ± 0.01); frequencies)
    deterministic = LineParametersProblem(reference.system, reference.earth; frequencies)
    mapped = map(x -> Float64(nominal(x)), problem)
    @test mapped isa LineParametersProblem{Float64}
    @test mapped.temperature == deterministic.temperature
    @test mapped.frequencies == deterministic.frequencies
    @test mapped.earth.layers == deterministic.earth.layers
    for name in (:id, :line_length, :placements, :regions, :terminals)
        @test getfield(mapped.system, name) == getfield(deterministic.system, name)
    end
    @test length(mapped.system.designs) == length(deterministic.system.designs)
    for (actual, expected) in zip(mapped.system.designs, deterministic.system.designs)
        @test actual.tree.items == expected.tree.items
        for name in (:id, :regions, :outer, :terminals, :assemblies, :repeats, :lays,
            :compactions, :fills)
            @test getfield(actual, name) == getfield(expected, name)
        end
    end

    @test_throws MethodError CoaxialFormulation(; Γ = complex(0.0, 1e-3 ± 1e-5))
    reference = last(REFERENCES[2])
    problem = LineParametersProblem(reference.system, homogeneous(rho = 0.1 ± 0.01);
        frequencies = reference.frequencies[1:50:end])
    result = @inferred compute(problem, CoaxialFormulation())
    @test result isa LineParameters{Measurement{Float64}}
end
