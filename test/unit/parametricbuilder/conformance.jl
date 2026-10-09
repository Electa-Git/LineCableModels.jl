@testitem "ParametricBuilder / Gridspace / structural realization" tags=[:unit, :measurements] setup=[
    UseEngineSupport] begin
    using Random
    using Measurements
    import LineCableModels.ParametricBuilder as PB

    struct DuplicateValue end
    (::DuplicateValue)(value)=(value, value)

    duplicated=PB.Gridspace{Tuple}(
        DuplicateValue(),
        (PB.Grid(10.0, PB.AbsoluteError(0.5)),)
    )
    duplicated_point=first(PB.points(duplicated))
    propagated=PB.materialize(duplicated_point)
    @test propagated[1] === propagated[2]
    @test Measurements.cov(propagated[1], propagated[2]) ==
          Measurements.uncertainty(propagated[1])^2

    rng=MersenneTwister(42)
    @test all(1:32) do _
        draw=LineCableModels.realize(rng, duplicated_point, :normal)
        draw[1] === draw[2]
    end

    independent=PB.Gridspace{Tuple}(
        tuple,
        (
            PB.Grid(10.0, PB.AbsoluteError(0.5)),
            PB.Grid(10.0, PB.AbsoluteError(0.5))
        );
        combine = :zip
    )
    independent_values=PB.materialize(first(PB.points(independent)))
    @test independent_values[1] !== independent_values[2]
    @test iszero(Measurements.cov(independent_values...))

    reused=PB.Grid(10.0, PB.AbsoluteError(0.5))
    reused_values=PB.materialize(first(PB.points(PB.Gridspace{Tuple}(
        tuple,
        (reused, reused);
        combine = :zip
    ))))
    @test reused_values[1] !== reused_values[2]
    @test iszero(Measurements.cov(reused_values...))
end
