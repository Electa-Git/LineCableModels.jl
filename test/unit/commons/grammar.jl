@testitem "Commons / grammar / finite sources and containers" tags=[:unit, :commons] begin
    using Random
    const C = LineCableModels.Commons

    # A Grid declares a deterministic or an uncertainty-bearing finite source.
    @test C.Grid((1, 2)) isa C.DeterministicGrid
    @test C.DeterministicGrid([1, 2]).vals == (1, 2)
    @test C.DeterministicGrid(:one).vals == (:one,)
    @test C.Grid((1.0, 2.0), 5.0) isa C.RelativeGrid
    @test C.Grid(1.0, C.AbsoluteError(0.5)) isa C.AbsoluteGrid
    @test C.Grid((1.0, 2.0), 5.0) isa C.AbstractUncertainGrid
    @test C.has_uncertainty(C.Grid(1.0, C.AbsoluteError(0.5)))
    @test !C.has_uncertainty(C.Grid((1.0, 2.0)))
    @test C.sample_uncertainty(Random.Xoshiro(1), C.UncertainValue(2.0, 0.0), :normal) == 2.0

    # `has_grid` finds a source at any depth of tuples, vectors, named tuples and pairs.
    @test C.has_grid(C.Grid((1, 2)))
    @test C.has_grid((a = (b = C.Grid((1, 2)),), c = 1))
    @test C.has_grid(Any[1, C.Grid((1, 2))])
    @test C.has_grid(:core => C.Grid((1, 2)))
    @test !C.has_grid((a = 1, b = [1.0, 2.0], c = (d = :x,)))
    @test !C.has_grid(Dict(:core => C.Grid((1, 2))))
    @test @inferred(C.has_grid((verbosity = (default = 0,), trace = Val(false)))) === false

    # A container with a source inside becomes the space of its reconstructions.
    named = C.Gridspace((a = C.Grid((1, 2)), b = :fixed))
    @test named isa C.Gridspace{NamedTuple{(:a, :b)}}
    @test collect(named) == [(a = 1, b = :fixed), (a = 2, b = :fixed)]
    paired = C.Gridspace(:core => C.Grid((1, 2)))
    @test paired isa C.Gridspace{Pair}
    @test collect(paired) == [:core => 1, :core => 2]
    vector = C.Gridspace([1.0, C.Grid((2, 3))])
    @test collect(vector) == [[1.0, 2.0], [1.0, 3.0]]
    @test eltype(first(vector)) === Float64
    @test_throws MethodError C.Gridspace(1.0)
end

@testitem "Commons / grammar / a scalar or a space of it" tags=[:unit, :commons] begin
    const C = LineCableModels.Commons
    const RealOrGrid = C.ScalarOrGrid{Real}

    # A value of T, a Gridspace of T and a non-empty Grid of T are members. The values of a
    # Grid can mix subtypes of T.
    @test 2.0 isa RealOrGrid
    @test C.Grid((1.0, 2.0)) isa RealOrGrid
    @test C.Grid((1, 2.5)) isa RealOrGrid
    @test C.Gridspace{Float64}(float, (C.Grid((1, 2)),)) isa RealOrGrid
    @test Shell(C.Grid((1e-3, 2e-3))) isa C.ScalarOrGrid{Shell}

    # The empty Grid, a vector, a Gridpoint, an uncertain grid and a space of another type
    # are not.
    @test !(C.Grid(()) isa RealOrGrid)
    @test !([1.0, 2.0] isa RealOrGrid)
    @test !(C.Gridpoint{Float64}(float, (1,)) isa RealOrGrid)
    @test !(C.Grid((1.0, 2.0), 5.0) isa RealOrGrid)
    @test !(C.Grid((:a, :b)) isa RealOrGrid)
    @test !(C.Gridspace{Symbol}(Symbol, (C.Grid(("a", "b")),)) isa RealOrGrid)

    # One method receives the scalar and the space, and lifts the space.
    half(x::C.ScalarOrGrid{Real}) = C.parameterize(Float64, value -> value / 2, (x,))
    @test @inferred(half(2.0)) === 1.0
    @test collect(@inferred(half(C.Grid((2.0, 4.0))))) == [1.0, 2.0]
end

@testitem "Commons / grammar / lifting and composition" tags=[:unit, :commons] begin
    using Random
    const C = LineCableModels.Commons

    # Without a source, `parameterize` returns the scalar call itself.
    @test C.parameterize(Tuple, tuple, (1, (a = 2,))) === (1, (a = 2,))
    space = C.parameterize(Tuple, tuple, (C.Grid((1, 2)), (a = C.Grid((3, 4)),)))
    @test collect(space) == [(1, (a = 3,)), (2, (a = 3,)), (1, (a = 4,)), (2, (a = 4,))]

    # `combine` names a composition. Zip pairs the direct parameters of one call, at any
    # depth of its containers.
    @test C.Composition(:product) === C.Product()
    @test C.Composition(:zip) === C.Zip()
    @test C.Composition(C.Zip()) === C.Zip()
    @test_throws ArgumentError C.Composition(:outer)
    zipped = C.parameterize(Tuple, tuple, (C.Grid((1, 2)), (a = C.Grid((3, 4)),));
        combine = :zip)
    @test collect(zipped) == [(1, (a = 3,)), (2, (a = 4,))]
    literal(x, y) = C.Gridspace{Tuple}(tuple, (x, y); combine = :zip)
    @test @inferred(literal(C.Grid((1, 2)), C.Grid((3, 4)))) isa C.Gridspace{Tuple}

    # A nested space of another call multiplies with the parameters of this one, so a zip
    # that includes it is refused. Pairing complete objects is a Grid of them.
    nested = C.Gridspace{Symbol}(identity, (C.Grid((:a, :b)),))
    refused = try
        C.Gridspace{Tuple}(tuple, (nested, C.Grid((1, 2))); combine = :zip)
    catch error
        error
    end
    @test refused isa ArgumentError
    @test occursin("pass a Grid of them", sprint(showerror, refused))
    @test length(C.Gridspace{Tuple}(tuple, (nested, C.Grid((1, 2))))) == 4
    single = C.Gridspace{Symbol}(identity, (C.Grid((:a,)),))
    @test collect(C.Gridspace{Tuple}(tuple, (single, C.Grid((1, 2))); combine = :zip)) ==
          [(:a, 1), (:a, 2)]

    # A point is an unresolved argument tuple that materializes or realizes.
    point = C.points(space, 2)
    @test point isa C.Gridpoint{Tuple}
    @test C.materialize(point) == (2, (a = 3,))
    @test_throws BoundsError C.points(space, 5)
    uncertain = C.Gridspace{Tuple}(tuple, (C.Grid(1.0, C.AbsoluteError(0.0)),))
    drawn = first(C.points(uncertain))
    arguments = C.realize_arguments(Random.Xoshiro(1), drawn, :normal)
    @test C.realize(drawn, arguments) == (1.0,)
end

@testitem "Commons / grammar / families and options" tags=[:unit, :commons] begin
    const C = LineCableModels.Commons
    const U = LineCableModels.Units

    # The constructor of a type declared in Commons lifts through its family method.
    options = C.FormulationOptions(C.Grid(((a = 1,), (a = 2,))))
    @test options isa C.Gridspace{C.FormulationOptions}
    @test [option.data for option in options] == [(a = 1,), (a = 2,)]
    keywords = C.FormulationOptions(; a = C.Grid((1, 2)), b = :fixed)
    @test [option.data for option in keywords] == [(a = 1, b = :fixed), (a = 2, b = :fixed)]
    @test @inferred(C.FormulationOptions(; a = 1)) === C.FormulationOptions((a = 1,))
    # The positional named tuple lifts like the keywords.
    positional = C.FormulationOptions((kron_reduction = C.Grid((true, false)),))
    @test positional isa C.Gridspace{C.FormulationOptions}
    @test [option.data for option in positional] ==
          [(kron_reduction = true,), (kron_reduction = false,)]
    nested = C.FormulationOptions((integration = (method = C.Grid((:quad, :gauss)),),))
    @test [option.data.integration.method for option in nested] == [:quad, :gauss]
    @test @inferred(C.FormulationOptions((a = 1,))).data === (a = 1,)
    @test_throws MethodError C.FormulationOptions(1.0)

    # Computation options control one computation, so they do not vary.
    refused = try
        C.ComputationOptions(; trace = C.Grid((true, false)))
    catch error
        error
    end
    @test refused isa ArgumentError
    @test occursin("computation options do not vary", sprint(showerror, refused))
    @test C.ComputationOptions(; trace = false).data == (trace = false,)

    # The generics of the grammar, which modules extend for their own types.
    for generic in (C.add!, C.build, C.homogenize, C.constitutive, C.line_length,
            C.resistance, C.inductance, C.capacitance, C.description)
        @test parentmodule(generic) === C
    end
    @test C.description(nothing) == "none"
    @test C.description(missing) == "method unavailable"
    @test U.label(U.quantity(C.R)) == "Series resistance"
    @test U.label(U.quantity(C.L)) == "Series inductance"
    @test U.label(U.quantity(C.C)) == "Shunt capacitance"
    @test U.family(U.quantity(C.R)) === Val(:series)
    @test U.family(U.quantity(C.C)) === Val(:shunt)
end

@testitem "Commons / grammar / builders" tags=[:unit, :commons] begin
    const C = LineCableModels.Commons

    # A builder is a singleton under the grammar type, and its one call method lifts it
    # over all its arguments and keywords.
    @test formula isa C.FormulaBuilder
    @test C.FormulaBuilder <: C.Builder{C.FormulaDefinition}
    # A builder's target is its `Target`, unless its owner specifies one for its subject. A
    # space selects by the type of its points.
    @test C.target(formula) === C.FormulaDefinition
    @test C.target(formula, :carson1926) === C.FormulaDefinition
    @test C.target(formula, C.Grid((:carson1926, :deri1973))) === C.FormulaDefinition
    # A literal identifier infers at its call site. The identifier is part of the type.
    @test @inferred((() -> formula(:carson1926))()) isa C.FormulaDefinition{:carson1926}
    earths = formula(C.Grid((:carson1926, :deri1973)))
    @test earths isa C.Gridspace{C.FormulaDefinition}
    @test [C.formula_id(selection) for selection in earths] == [:carson1926, :deri1973]
    orders = formula(:default; order = C.Grid((:before, :after)))
    @test [NamedTuple(selection).order for selection in orders] == [:before, :after]
    paired = formula(C.Grid((:carson1926, :deri1973)); order = C.Grid((:before, :after)),
        combine = :zip)
    @test [(C.formula_id(s), NamedTuple(s).order) for s in paired] ==
          [(:carson1926, :before), (:deri1973, :after)]
end

@testitem "Commons / grammar / buffers extend the calculation record" tags=[:unit, :engine] begin
    const Engine=LineCableModels.Engine
    plan=Engine.ReductionPlan([1, 2, 2, 0]; reduce_bundle=true, kron_reduction=true,
        ideal_transposition=false)
    existing=(destination = zeros(ComplexF64, 2, 2),)
    buffers=LineCableModels.Commons.initialize_buffers(plan, ComplexF64, (;), (;), existing)
    @test keys(buffers) == (:destination, :reduction)
    @test buffers.destination === existing.destination
    @test keys(buffers.reduction) ==
          (:ordered, :potential, :factor, :coupling, :right_hand_side, :identity)
    @test all(matrix -> matrix isa Matrix{ComplexF64}, buffers.reduction)
    @test size(buffers.reduction.ordered) == (4, 4)
    @test size(buffers.reduction.potential) == (length(plan.keep), length(plan.keep))
    @test size(buffers.reduction.factor) == (length(plan.eliminate), length(plan.eliminate))
end
