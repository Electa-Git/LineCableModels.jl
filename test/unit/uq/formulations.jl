@testitem "UQ / computation options / validated constructors and concrete tuples" tags=[:unit] begin
    using Random
    inner = CableConstantsFormulation()
    normalize = LineCableModels.computation_options
    @test LineCableModels.ComputationOptions === NamedTuple

    for owner in (Combinatorial, LinearError, MonteCarlo)
        keyword = @inferred owner(inner; options=(retain_details=true,))
        positional = @inferred owner(inner, (retain_details=true,))
        @test keyword.options === positional.options
        @test isconcretetype(typeof(keyword.options))
        @test normalize(owner, keyword.options) === keyword.options
        @test_throws ArgumentError owner(inner; options=(unused=true,))
        @test_throws ArgumentError owner(inner, (unused=true,))
        @test_throws ArgumentError owner(inner, (retain_details=1,))
        # Explicit type parameters must not restore the generated raw constructor.
        @test_throws MethodError owner{typeof(inner), typeof((unused=true,))}(
            inner, (unused=true,))
    end
    @test_throws ArgumentError Combinatorial(Grid((1, 2)), (;))
    @test_throws ArgumentError Combinatorial(
        Gridspace{Int}(identity, (Grid((1, 2)),)), (;))

    plain = (trials=Int16(8), confidence=0.9f0, cdf_tol=0.1f0,
        distribution=:uniform, seed=Int16(42), return_samples=true,
        return_histograms=true, bins=Int16(2), retain_details=true,
        on_error=:retry, max_failures=Int16(5))
    normalized = @inferred normalize(MonteCarlo, plain)
    @test normalized.trials === 8
    @test normalized.seed === UInt64(42)
    @test normalized.confidence === Float64(0.9f0)
    @test normalized.cdf_tol === Float64(0.1f0)
    @test normalized.bins === 2
    @test normalized.max_failures === 5
    @test (@inferred MonteCarlo(inner; options=plain)).options === normalized
    @test (@inferred MonteCarlo(inner; plain...)).options === normalized
    @test (@inferred MonteCarlo(inner; trials=8, options=(seed=42,))).options ===
          (@inferred MonteCarlo(inner; options=(trials=8, seed=42))).options
    @test_throws ArgumentError MonteCarlo(inner; trials=8, options=(trials=8,))
    @test_throws ArgumentError MonteCarlo(inner; unused=true)
    @test_throws ArgumentError normalize(MonteCarlo, (on_error=:retry,))
    @test_throws ArgumentError normalize(MonteCarlo, (on_error=:ignore,))
    @test_throws ArgumentError normalize(MonteCarlo, (distribution=:unsupported,))

    for key in (:trials, :bins, :max_failures)
        for value in (true, 0, -1, 1.5, "2", big(typemax(Int)) + 1)
            options = NamedTuple{(key,)}((value,))
            @test_throws ArgumentError normalize(MonteCarlo, options)
        end
    end
    @test_throws ArgumentError normalize(MonteCarlo, (max_failures=nothing,))
    for key in (:confidence, :cdf_tol)
        for value in (true, 0, 1, -0.5, NaN, Inf, "0.5",
                big(1) - eps(BigFloat), big(2)^(-2000))
            @test_throws ArgumentError normalize(MonteCarlo, NamedTuple{(key,)}((value,)))
        end
    end
    for key in (:return_samples, :return_histograms, :retain_details)
        @test_throws ArgumentError normalize(MonteCarlo, NamedTuple{(key,)}((1,)))
    end
    for seed in (true, -1, 1.5, "42", big(typemax(UInt64)) + 1)
        @test_throws ArgumentError normalize(MonteCarlo, (; seed))
    end
    @test normalize(MonteCarlo, (seed=0,)).seed === UInt64(0)
    @test normalize(MonteCarlo, (seed=typemax(UInt64),)).seed === typemax(UInt64)
    automatic = @inferred MonteCarlo(inner)
    @test automatic.options.trials === nothing
    @test automatic.options.seed === nothing
    @test automatic.options.bins === nothing

    sampler = (rng, center, deviation) -> center + deviation * randn(rng)
    custom = @inferred MonteCarlo(inner; options=(distribution=sampler,))
    @test custom.options.distribution === sampler
    @test fieldtype(typeof(custom.options), :distribution) === typeof(sampler)
end
