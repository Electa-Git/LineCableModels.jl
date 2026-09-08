@testitem "Engine / numerical declarations are scoped by complete binding and actual provider" tags=[:unit] begin
    const E = LineCableModels.Engine
    const II = E.InternalImpedance
    const EI = E.EarthImpedance
    const FM = LineCableModels.FormulaMethod
    internal = FM(Val(:default), II.internal_impedance, Val(:outer))
    external = FM(Val(:default), EI.earth_impedance, Val(:self), Val(1), Val(1))
    @test computation_options(internal) == (;)
    @test computation_options(external).integration.method === :quad
    @test_throws ArgumentError computation_options(internal, (integration = (method = :quad,),))
    @test_throws ArgumentError computation_options(FM(Val(:Unspecified), II.internal_impedance, Val(:outer)))
    @test_throws ArgumentError computation_options(LineCableModelsCoaxial, (integration_method = :quad,))

    pair = E.EarthPair(1, 1, (1.0, 1.0), 0.0, (1, 1); radius = 0.01)
    replacement = (functor, pair, workspace) -> complex(3.0, 4.0)
    missing = EI.Formula(:default; hooks = (contribution = replacement,))
    @test_throws ArgumentError validate(missing, pair)
    @eval LineCableModels.computation_options(
        ::FM{:default, typeof(EI.earth_impedance)}, ::$(typeof(replacement))) = (;)
    admitted = validate(missing, pair)
    @test isempty(admitted.options)
    @test_throws ArgumentError validate(
        EI.Formula(:default;
            hooks = (contribution = replacement,),
            options = (integration = (method = :quad,),)),
        pair)
    @test_throws ArgumentError validate(missing,
        E.EarthPair(1, 2, (1.0, -1.0), 1.0, (1, 2)))
    rho = [Inf, 100.0]
    epsilon = 8.8541878128e-12 .* [1, 10]
    mu = fill(4pi*1e-7, 2)
    @test missing(rho, epsilon, mu, 100.0im, pair)() === 3.0 + 4.0im
    @test computation_options(external).integration.method === :quad
end

@testitem "Engine / public internal surfaces consume spectral options and retained resources" tags=[:unit] setup=[TestFixtures] begin
    using QuadGK
    const E=LineCableModels.Engine
    const II=E.InternalImpedance
    const FM=LineCableModels.FormulaMethod
    seen=Tuple[]
    outer=(functor,
        workspace)->begin
        push!(seen, (functor.options.integration.method, workspace))
        integral=E.SpectralIntegral(Val(:cosine), λ->complex(exp(-λ)),
            (height = 1.0, separation = 0.0), 1.0)
        E.integrate(functor.options.integration.method, integral,
            functor.options.integration.options, workspace)*1e-4
    end
    @eval LineCableModels.computation_options(
        ::FM{:default, typeof(II.internal_impedance), Tuple{Val{:outer}}},
        ::$(typeof(outer))) = (integration = (method = :quad, options = (;)),)
    @test_throws ArgumentError II.Formula(:default; options = (integration = (method = :quad,),))
    resources=(segments = alloc_segbuf(Float64, ComplexF64, Float64; size = 32),
        images = ComplexF64[], exponents = ComplexF64[])
    reference=II.surface_impedances(
        II.Formula(:default), 0.008, 0.01, 1.7241e-8, 1.0, 100.0im)
    for method in (:quad, :trapz, :cim)
        selected=II.Formula(formula(:default; hooks = (outer = outer,),
            options = (integration = (method = method,),)))
        values=II.surface_impedances(selected, 0.008, 0.01, 1.7241e-8, 1.0, 100.0im;
            workspace = resources)
        @test values.outer ≈ 5e-5 rtol=3e-6
        @test values.inner == reference.inner
        @test values.mutual == reference.mutual
        @test last(seen) === (Val(method), resources)
        @test isempty(selected.options.inner) && isempty(selected.options.mutual)
    end
    # The same override and numerical selection must reach the backend's internal terms.
    problem=TestFixtures.line_parameters_problem(frequencies = [50.0])
    results=map((:quad, :trapz, :cim)) do method
        empty!(seen)
        result=compute(problem,
            Formulation(internal_impedance = formula(:default;
                hooks = (outer = outer,), options = (integration = (method = method,),))))
        @test !isempty(seen)
        @test all(record -> record[1] === Val(method), seen)
        result
    end
    @test results[2].Z.values ≈ results[1].Z.values rtol=3e-6
    @test results[3].Z.values ≈ results[1].Z.values rtol=3e-6
end

@testitem "Engine / surface current basis reproduces concentric wall contributions" tags=[:unit] begin
    using LinearAlgebra
    const II = LineCableModels.Engine.InternalImpedance
    selected = II.Formula(:default)
    for frequency in (1e-4, 50.0, 1e5)
        coefficients = II.surface_impedances(selected, 0.008, 0.01, 1.7241e-8, 1.0,
            complex(0.0, 2pi*frequency))
        W = [coefficients.inner coefficients.mutual; coefficients.mutual coefficients.outer]
        # Terminal currents are (contained metal, enclosing wall). Surface
        # currents are (-contained, contained + wall), exactly the current map.
        B = [-1.0 1.0; 0.0 1.0]
        lifted = B * W * transpose(B)
        @test lifted[1, 1] ≈ coefficients.inner - 2coefficients.mutual + coefficients.outer
        @test lifted[1, 2] ≈ coefficients.outer - coefficients.mutual
        @test lifted[2, 2] == coefficients.outer
        if frequency == 1e-4
            resistance = 1.7241e-8 / (pi * (0.01^2 - 0.008^2))
            @test real(lifted[2, 2]) ≈ resistance rtol=1e-10
            @test abs(real(lifted[1, 1])) < 1e-10resistance
            @test abs(real(lifted[1, 2])) < 1e-10resistance
        end
    end
end

@testitem "Earth / artificial material values are distinct from source restrictions" tags=[:unit] setup=[FormulaContractModels] begin
    const M=FormulaContractModels
    material=M.EP.EarthMaterial(100.0, -10.0, 1.0)
    @test material.eps_r == -10
    @test M.EP.EarthLayer <: M.EP.AbstractEarthLayer <: M.EP.AbstractEarthModel
    @test M.EP.EarthMaterial <: M.EP.AbstractEarthMaterial <:
          LineCableModels.AbstractMaterial
    rho=[Inf, 100.0, 200.0]
    epsilon=8.8541878128e-12 .* [1, -10, -20]
    mu=fill(4pi*1e-7, 3)
    pair=M.E.EarthPair(1, 2, (-0.25, -1.5), 1.0, (2, 3))
    @test isfinite(M.selection(M.EI)(rho, epsilon, mu, 100.0im, pair;
        thickness = [Inf, 0.5, Inf])())
    @test_throws DomainError M.EI.Formula(:default)(
        rho[1:2], epsilon[1:2], mu[1:2], 100.0im,
        M.E.EarthPair(1, 2, (-0.25, -1.5), 1.0, (2, 2)))
end

@testitem "Engine / internal consumers request only their actual surface kinds" tags=[:unit] begin
    const II = LineCableModels.Engine.InternalImpedance
    const FM = LineCableModels.FormulaMethod
    II.internal_impedance(::Val{:ManufacturedOuter}, ::Val{:outer}, functor,
        workspace) = complex(functor.state.rho / (pi * functor.state.radius^2))
    LineCableModels.computation_options(::FM{
        :ManufacturedOuter, typeof(II.internal_impedance), Tuple{Val{:outer}}}) = (;)
    function (formula::II.Formula{:ManufacturedOuter})(r_in, r_ex, rho, mu_r, jω)
        state = (rho = rho, radius = r_ex)
        return II.Functor{
            :ManufacturedOuter, typeof(formula.binding), typeof(formula.hooks),
            typeof(state), typeof(formula.options)}(
            formula.binding, formula.hooks, state, formula.options)
    end
    binding = (outer = FM(Val(:ManufacturedOuter), II.internal_impedance, Val(:outer)),)
    controls = (outer = (;),)
    selected = II.Formula{
        :ManufacturedOuter, typeof(binding), NamedTuple{()}, NamedTuple{()},
        typeof(controls), Tuple{}}(binding, (;), (;), controls, ())
    @test validate(selected, (:outer,)) === selected
    @test_throws ArgumentError validate(selected, (:inner, :outer, :mutual))
    @test_throws ArgumentError II.surface_impedances(
        selected, 0.0, 0.01, 1.7e-8, 1.0, 100im)
    copper = Material(kind = :conductor, rho = 1.7e-8)
    insulation_material = Material(kind = :insulator, rho = Inf, eps_r = 2.3)
    design = build(CableDesign,
        "outer-only",
        terminal(:core,
            solid(copper, Disk(0.01)), insulation(insulation_material; t = 0.002)))
    system = build(LineCableSystem, design, Pose2(0.0, 5.0); connections = Dict(:core => 1))
    problem = LineParametersProblem(system; frequencies = [50.0], earth_props = homogeneous(rho = 100.0))
    result = compute(problem, Formulation(internal_impedance = selected))
    @test all(isfinite, result.Z.values)
    @test keys(details(result).formulations.numerical.internal_impedance) == (:outer,)
    # A configured inner override is an error when the geometry consumes only outer.
    inner = (functor, workspace) -> 1.0im
    @eval LineCableModels.computation_options(
        ::FM{:default, typeof(II.internal_impedance), Tuple{Val{:inner}}}, ::$(typeof(inner))) = (;)
    @test_throws ArgumentError compute(
        problem, Formulation(internal_impedance =
        formula(:default; hooks = (inner = inner,))))
end

@testitem "Engine / each metal prepares shared surface state once per frequency" tags=[:unit] setup=[TestFixtures] begin
    const II=LineCableModels.Engine.InternalImpedance
    const FM=LineCableModels.FormulaMethod
    const preparations=Tuple[]
    const evaluations=Tuple[]
    function II.internal_impedance(::Val{:ManufacturedSurfaces},
            kind::Union{Val{:inner}, Val{:outer}, Val{:mutual}}, functor, workspace)
        push!(evaluations, (functor.state.serial, kind))
        return kind===Val(:mutual) ? 0.5+0.1im : 2.0+1.0im
    end
    LineCableModels.computation_options(::FM{
        :ManufacturedSurfaces, typeof(II.internal_impedance)})=(;)
    function (formula::II.Formula{:ManufacturedSurfaces})(r_in, r_ex, rho, mu_r, jω)
        push!(preparations, (r_in, r_ex, rho, mu_r, jω))
        state=(serial = length(preparations),)
        return II.Functor{
            :ManufacturedSurfaces, typeof(formula.binding), typeof(formula.hooks),
            typeof(state), typeof(formula.options)}(
            formula.binding, formula.hooks, state, formula.options)
    end
    kinds=(:inner, :outer, :mutual)
    binding=NamedTuple{kinds}(map(
        kind->FM(Val(:ManufacturedSurfaces), II.internal_impedance, Val(kind)), kinds))
    controls=NamedTuple{kinds}(map(_->(;), kinds))
    selected=II.Formula{
        :ManufacturedSurfaces, typeof(binding), NamedTuple{()}, NamedTuple{()},
        typeof(controls), Tuple{}}(binding, (;), (;), controls, ())
    problem=TestFixtures.line_parameters_problem(frequencies = [50.0, 100.0])
    result=compute(problem, Formulation(internal_impedance = selected))
    expected=sum(length(design.terminal_order)
    for design in problem.system.designs)*length(problem.frequencies)
    @test length(preparations) == expected
    @test allunique(evaluations)
    @test count(record -> record[2] === Val(:outer), evaluations) == expected
    @test any(record -> record[2] === Val(:inner), evaluations)
    @test any(record -> record[2] === Val(:mutual), evaluations)
    @test all(isfinite, result.Z.values)
end
