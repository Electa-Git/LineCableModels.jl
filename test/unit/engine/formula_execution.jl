@testitem "Engine / explicit equation recipes are complete only for required cases" tags=[:unit] setup=[TestFixtures] begin
    const E = LineCableModels.Engine
    omitted = Formulation()
    absent = Formulation(earth_impedance=nothing, earth_admittance=nothing,
        internal_impedance=nothing)
    for slot in (:earth_impedance, :earth_admittance, :internal_impedance)
        @test formula_id(getproperty(absent.methods, slot)) ==
              formula_id(getproperty(omitted.methods, slot))
    end
    buried = TestFixtures.three_bare_wires_problem(heights=(-1.0,-1.0,-1.0), frequencies=[50.0])
    overhead = TestFixtures.three_bare_wires_problem(frequencies=[50.0])
    options = (reduce_bundle=false, kron_reduction=false, ideal_transposition=false)
    partial = Formulation(earth_impedance=(earth=formula(:unified),),
        earth_admittance=(earth=formula(:unified),),
        internal_impedance=(outer=formula(:default),); options)
    expected = compute(buried, Formulation(;options))
    actual = compute(buried, partial)
    @test Z(actual) == Z(expected)
    @test Y(actual) == Y(expected)
    @test_throws ArgumentError compute(overhead, partial)
    @test_throws ArgumentError compute(buried, Formulation(earth_impedance=(;)))
    @test_throws ArgumentError compute(buried,
        Formulation(earth_impedance=(earth=nothing,)))
    excess = Formulation(
        earth_impedance=(air=formula(:carson1926), earth=formula(:unified)),
        earth_admittance=(air=formula(:wise1948), earth=formula(:unified)); options)
    @test Z(compute(buried, excess)) == Z(expected)
    @test Y(compute(buried, excess)) == Y(expected)
    @test_throws ArgumentError Formulation(earth_impedance=(soil=formula(:unified),))
    @test_throws ArgumentError Formulation(earth_admittance=formula(:unknown))
    disabled = Formulation(earth_properties=nothing, temperature_dependence=nothing)
    @test disabled.methods.earth_properties === nothing
    @test disabled.methods.temperature_dependence === nothing
end

@testitem "Engine / registered author earth equations are explicit numerical stubs" tags=[:unit] begin
    const E = LineCableModels.Engine
    for (owner, operation) in ((E.EarthImpedance, E.EarthImpedance.earth_impedance),
                              (E.EarthAdmittance, E.EarthAdmittance.earth_potential_coefficient))
        @test formula_id(owner.Formula(:default)) === :unified
        for identifier in setdiff(owner.formulas(), (:default, :unified))
            selected = owner.Formula(identifier)
            @test formula_id(selected) === identifier
            @test !isempty(description(selected))
            # Stubs fail before reading a functor or executing any numerics.
            @test_throws ArgumentError operation(selected, Val(:mutual), Val(1), Val(1),
                nothing, nothing, nothing)
            @test_throws ArgumentError operation(selected, Val(:mutual), Val(2), Val(2),
                nothing, nothing, nothing)
        end
    end
end

@testitem "Engine / prescribed propagation is a unified formulation option" tags=[:unit] setup=[TestFixtures] begin
    const E = LineCableModels.Engine
    problem = TestFixtures.three_bare_wires_problem(frequencies=[50.0, 500.0])
    @test_throws MethodError LineParametersProblem(problem.system;
        earth_props=problem.earth_props, frequencies=problem.frequencies, Γ=[0im,0im])
    options = (reduce_bundle=false, kron_reduction=false, ideal_transposition=false)
    make(gamma) = Formulation(
        earth_impedance=formula(:unified; options=(Γ=gamma,)),
        earth_admittance=formula(:unified; options=(Γ=gamma,)); options)
    implicit = compute(problem, Formulation(;options))
    zero_gamma = compute(problem, make(0.0im))
    @test Z(zero_gamma) == Z(implicit)
    @test Y(zero_gamma) == Y(implicit)
    scalar = compute(problem, make(1e-4im))
    aligned = compute(problem, make(fill(1e-4im, 2)))
    @test Z(scalar) == Z(aligned)
    @test Y(scalar) == Y(aligned)
    @test_throws DimensionMismatch compute(problem, make([1e-4im]))
    @test_throws ArgumentError make(NaN)
    @test_throws ArgumentError make([0im,complex(Inf)])
    for invalid in (true, [], (0, 1), "0", [0, NaN])
        @test_throws ArgumentError make(invalid)
    end
    @test_throws ArgumentError Formulation(
        earth_impedance=formula(:unified; parameters=(Γ=1e-4im,)))
    method = LineCableModels.FormulaMethod(make(0).methods.earth_impedance,
        E.EarthImpedance.earth_impedance, Val(:self), Val(1), Val(1))
    @test formulation_options(method).data.Γ == 0
    @test formulation_options(method, FormulationOptions(Γ=2e-4im)).data.Γ == 2e-4im
    @test_throws ArgumentError formulation_options(method, FormulationOptions(Γ=NaN))
    @test_throws ArgumentError formulation_options(method, FormulationOptions(integration=1))
    using Measurements
    rho = measurement(0.1, 0.001)
    material_problem = TestFixtures.three_bare_wires_problem(; rho,
        radius=oftype(rho, 0.0425), frequencies=[50.0, 500.0])
    uncertain = compute(material_problem, make(1e-4im))
    @test eltype(Z(uncertain)) === Complex{Measurement{Float64}}
    @test any(>(0), uncertainty.(real.(Z(uncertain))))
    @test E.same_physical_state([rho], [rho])
    @test !E.same_physical_state([rho], [measurement(0.1, 0.001)])
    record = LineCableModels.ImportExport.serialize_value(aligned)
    restored = LineCableModels.ImportExport.deserialize_value(record)
    @test details(restored).data.formulations.methods.earth_impedance.options.Γ == fill(1e-4im, 2)
    prescribed = [1e-4im, 2e-4im]
    sweep = compute(problem, make(prescribed))
    # Run scalar samples in reverse call order. Prescriptions follow the
    # existing frequency coordinates; neither side is sorted independently.
    for index in reverse(eachindex(problem.frequencies))
        single = LineParametersProblem(problem.system; earth_props=problem.earth_props,
            temperature=problem.temperature, frequencies=[problem.frequencies[index]])
        result = compute(single, make(prescribed[index]))
        @test Z(result)[:, :, 1] == Z(sweep)[:, :, index]
        @test Y(result)[:, :, 1] == Y(sweep)[:, :, index]
    end
    setprecision(BigFloat,128) do
        mixed_precision = compute(problem, LineParametersFormulation[
            make(1e-4im), make(big"0.0001"*im)])
        @test eltype(Z(mixed_precision[1])) === ComplexF64
        @test eltype(Z(mixed_precision[2])) === Complex{BigFloat}
        @test Z(mixed_precision[2]) == Z(compute(problem, make(big"0.0001"*im)))
    end
end

@testitem "Engine / completed prescriptions are detached without replacing uncertainty sources" tags=[:unit] setup=[TestFixtures, FormulaFixtures] begin
    using Measurements
    E, IE = LineCableModels.Engine, LineCableModels.ImportExport
    problem = TestFixtures.three_bare_wires_problem(frequencies=[50., 500.])
    original = [1e-4+2e-4im, 3e-4-1e-4im]
    supplied = copy(original)
    selection = formula(:unified; options=(Γ=supplied,))
    formulation = Formulation(earth_impedance=selection, earth_admittance=selection)
    first_result = compute(problem, formulation)
    first_Z, first_Y = copy(Z(first_result)), copy(Y(first_result))
    record = details(first_result).data.formulations
    for family in (:earth_impedance, :earth_admittance)
        @test getproperty(record.methods, family).options.Γ == original
        @test getproperty(record.requested, family).options.Γ == original
        @test getproperty(record.methods, family).options.Γ !==
              getproperty(formulation.methods, family).options.data.Γ
    end
    for family in (:earth_impedance, :earth_admittance)
        getproperty(formulation.methods, family).options.data.Γ[2] = 7e-4im
    end
    supplied[2] = 7e-4im
    second_result = compute(problem, formulation)
    @test Z(first_result) == first_Z && Y(first_result) == first_Y
    @test Z(second_result) != first_Z || Y(second_result) != first_Y
    restored = IE.deserialize_value(IE.serialize_value(first_result))
    for retained in (first_result, restored), family in (:earth_impedance, :earth_admittance)
        data = details(retained).data.formulations
        @test getproperty(data.methods, family).options.Γ == original
        @test getproperty(data.requested, family).options.Γ == original
        @test Z(retained) == first_Z && Y(retained) == first_Y
    end
    rho = measurement(1.0, 0.01)
    custom = FormulaFixtures.selection(E.EarthImpedance; layers=2:2, scale=rho)
    retained = NamedTuple(Formulation(earth_impedance=custom))
    @test retained.methods.earth_impedance.parameters.scale === rho
    @test uncertainty(retained.methods.earth_impedance.parameters.scale-rho) == 0
end

@testitem "Engine / prescribed options govern coupled sharing without entering physical UQ" tags=[:unit] setup=[TestFixtures] begin
    using Measurements
    E = LineCableModels.Engine
    prescribed = [1e-4+2e-4im, 3e-4-1e-4im]
    problem = TestFixtures.three_bare_wires_problem(frequencies=[50., 500.])
    select(gamma; controls=(;)) = formula(:unified;
        options=(Γ=gamma, integration=(method=:quad, options=controls)))
    form(z, p) = Formulation(earth_impedance=z, earth_admittance=p;
        options=(reduce_bundle=false, kron_reduction=false, ideal_transposition=false))
    joined = compute(problem, form(select(prescribed), select(copy(prescribed)));
        options=(trace=true,))
    changed = prescribed .* 2
    separate = compute(problem, form(select(prescribed), select(changed)); options=(trace=true,))
    other = compute(problem, form(select(changed), select(changed)); options=(trace=true,))
    @test length(details(separate).data.trace.integrals) ==
          2length(details(joined).data.trace.integrals)
    @test Z(separate) == Z(joined)
    @test Y(separate) == Y(other)
    controls = compute(problem, form(select(prescribed), select(prescribed; controls=(rtol=1e-9,)));
        options=(trace=true,))
    @test length(details(controls).data.trace.integrals) ==
          2length(details(joined).data.trace.integrals)
    selected = form(select(prescribed), select(prescribed))
    space = Gridspace{LineParametersProblem}(
        rho -> TestFixtures.three_bare_wires_problem(; rho,
            radius=oftype(rho, 0.0425), frequencies=[50., 500.]),
        (Grid(0.1, AbsoluteError(0.001)),))
    sampled = compute(ParametricProblem(space),
        MonteCarlo(selected; trials=3, seed=314, distribution=:normal,
            return_samples=true, retain_details=true))
    @test sampled.trial_counts == [3]
    @test all(isempty, details(sampled).data.failures)
    for trial in only(details(sampled).data.trials), family in (:earth_impedance, :earth_admittance)
        @test getproperty(trial.data.formulations.methods, family).options.Γ == prescribed
    end
    linear = compute(ParametricProblem(space), LinearError(selected))
    @test any(>(0), uncertainty.(real.(Z(first(linear)))))
    @test selected.methods.earth_impedance.options.data.Γ == prescribed
end

@testitem "Engine / numerical declarations follow selected types and indexed equations" tags=[:unit] setup=[FormulaFixtures] begin
    const E=LineCableModels.Engine
    const II=E.InternalImpedance
    const EI=E.EarthImpedance
    const FM=LineCableModels.FormulaMethod
    const M=FormulaFixtures
    internal=FM(II.Formula(:default), II.internal_impedance, Val(:outer))
    external=FM(EI.Formula(:default), EI.earth_impedance, Val(:self), Val(1), Val(1))
    @test formulation_options(internal)==FormulationOptions()
    @test formulation_options(external).data.integration.method === :quad
    @test_throws ArgumentError formulation_options(
        internal, FormulationOptions(integration = (method = :quad,)))
    @test_throws ArgumentError computation_options(
        LineCableModelsCoaxial, ComputationOptions((integration_method = :quad,)))
    # The custom selection uses existing admission and equation generics. It is
    # not a changed implementation of the built-in's claimed scientific identity.
    custom=M.selection(EI; layers = 2:2)
    pair=E.EarthPair(1, 2, (-1.0, -1.0), 1.0, (2, 2))
    bound=validate(custom, pair)
    @test bound.equation.selection === custom
    @test isempty(bound.options.data)
    @test_throws ArgumentError validate(
        M.selection(EI; layers = 2:2,
            options = (integration = (method = :quad,),)), pair)
    rho=[Inf, 100.0]
    epsilon=8.8541878128e-12 .* [1, 10]
    mu=fill(4pi*1e-7, 2)
    value=custom(rho, epsilon, mu, 100.0im, pair; thickness = [Inf, Inf])()
    @test isfinite(value)
    @test formula_id(custom) !== formula_id(external.selection)
    @test formulation_options(external).data.integration.method === :quad
end

@testitem "Engine / required indexed consumers alone initialize formula storage" tags=[:unit] setup=[TestFixtures, FormulaFixtures] begin
    const E=LineCableModels.Engine
    const M=FormulaFixtures
    buried=TestFixtures.three_bare_wires_problem(heights=(-1.,-1.,-1.), frequencies=[50.])
    active=M.selection(E.EarthImpedance;layers=2:2)
    unused=M.selection(E.EarthImpedance;layers=3:3,
        options=(integration=(method=:quad,),))
    potential=M.selection(E.EarthAdmittance;layers=2:2)
    selected=Formulation(earth_impedance=(air=unused,earth=active),earth_admittance=potential)
    first_result=compute(buried,selected)
    @test isempty(unused.initialized)
    required, numerical=only(active.initialized)
    @test required == fill((2,2),9)
    @test isempty(numerical.segments)
    saved=copy(Z(first_result))
    second_result=compute(buried,selected)
    @test Z(first_result) == saved == Z(second_result)
    @test length(active.initialized)==2
    @test active.initialized[2][2].segments !== numerical.segments
    @test isempty(unused.initialized)
    # The same concrete formula needs numerical storage for its (1,1) equation.
    overhead=TestFixtures.three_bare_wires_problem(frequencies=[50.])
    compute(overhead,Formulation(earth_impedance=active,earth_admittance=potential))
    required, numerical=last(active.initialized)
    @test required == fill((1,1),9)
    @test !isempty(numerical.segments)
end

@testitem "Engine / formula initialization preserves another owner's buffer identity" tags=[:unit] setup=[FormulaFixtures] begin
    const E=LineCableModels.Engine
    buffers=(destination = zeros(ComplexF64, 2, 2),)
    @test E.initialize_buffers((nothing,), Float64, (;), (;), buffers) === buffers
    @test_throws ArgumentError E.initialize_buffers(
        (FormulaFixtures.BufferReplacement(),),
        Float64, (;), (;), buffers)
end

@testitem "Engine / public internal surfaces consume spectral options and initialized workspace" tags=[:unit] setup=[
    TestFixtures, FormulaFixtures] begin
    using QuadGK
    const E=LineCableModels.Engine
    const II=E.InternalImpedance
    const M=FormulaFixtures
    @test_throws ArgumentError II.Formula(:default; options = (integration = (method = :quad,),))
    standalone_workspace=(buffers = (quadrature = E.integration_workspace(Float64, ComplexF64),),)
    base=II.Formula(:default)
    args=(0.008, 0.01, 1.7241e-8, 1.0, 100.0im)
    reference=II.surface_impedances(base, args...)
    for method in (:quad,)
        selected=M.SpectralSurface(method)
        surfaces=(inner = base, outer = selected, transfer = base)
        values=II.surface_impedances(surfaces, args...; workspace = standalone_workspace)
        @test values.outer ≈ 5e-5 rtol=3e-6
        @test values.inner == reference.inner
        @test values.transfer == reference.transfer
        @test last(selected.seen) === (Val(method), standalone_workspace)
    end
    problem=TestFixtures.line_parameters_problem(frequencies = [50.0])
    results=map((:quad,)) do method
        selected=M.SpectralSurface(method)
        result=compute(problem,
            Formulation(internal_impedance =
            (inner = base, outer = selected, transfer = base)))
        @test !isempty(selected.seen)
        @test first(selected.seen)[1] === :initialize
        evaluations=filter(record->record[1]===Val(method), selected.seen)
        @test !isempty(evaluations)
        @test all(record -> record[2] === first(evaluations)[2], evaluations)
        @test first(evaluations)[2] isa E.LineParametersWorkspace
        @test !isempty(first(evaluations)[2].buffers.quadrature.segments)
        result
    end
    @test all(isfinite, only(results).Z.values)
    @test all(isfinite, only(results).Y.values)

    selected=M.SpectralSurface()
    local_problem=CableConstantsProblem(TestFixtures.coaxial_design())
    local_formulation=CableConstantsFormulation(
        internal_impedance = (inner = base, outer = selected, transfer = base))
    local_result=compute(local_problem, local_formulation)
    @test first(selected.seen)[1] === :initialize
    local_evaluations=filter(record->record[1]===Val(:quad), selected.seen)
    @test !isempty(local_evaluations)
    @test first(local_evaluations)[2] isa E.CableConstantsWorkspace
    @test all(record->record[2] === first(local_evaluations)[2], local_evaluations)
    @test all(isfinite, resistance(local_result))
end

@testitem "Engine / surface current basis reproduces concentric wall contributions" tags=[:unit] begin
    using LinearAlgebra
    const II = LineCableModels.Engine.InternalImpedance
    selected = II.Formula(:default)
    for frequency in (1e-4, 50.0, 1e5)
        coefficients = II.surface_impedances(selected, 0.008, 0.01, 1.7241e-8, 1.0,
            complex(0.0, 2pi*frequency))
        W = [coefficients.inner coefficients.transfer;
             coefficients.transfer coefficients.outer]
        # Terminal currents are (contained metal, enclosing wall). Surface
        # currents are (-contained, contained + wall), exactly the current map.
        B = [-1.0 1.0; 0.0 1.0]
        lifted = B * W * transpose(B)
        @test lifted[1, 1] ≈
              coefficients.inner - 2coefficients.transfer + coefficients.outer
        @test lifted[1, 2] ≈ coefficients.outer - coefficients.transfer
        @test lifted[2, 2] == coefficients.outer
        if frequency == 1e-4
            resistance = 1.7241e-8 / (pi * (0.01^2 - 0.008^2))
            @test real(lifted[2, 2]) ≈ resistance rtol=1e-10
            @test abs(real(lifted[1, 1])) < 1e-10resistance
            @test abs(real(lifted[1, 2])) < 1e-10resistance
        end
    end
end

@testitem "Earth / artificial material values are distinct from source restrictions" tags=[:unit] setup=[FormulaFixtures] begin
    const M=FormulaFixtures
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

@testitem "Engine / internal consumers request only their actual surface kinds" tags=[:unit] setup=[FormulaFixtures] begin
    const II=LineCableModels.Engine.InternalImpedance
    selected=FormulaFixtures.SurfaceLaw(kinds = (:outer,))
    @test keys(II.surface_impedances(
        selected, 0.0, 0.01, 1.7e-8, 1.0, 100im)) == (:outer,)
    @test_throws ArgumentError II.surface_impedances(
        selected, 0.008, 0.01, 1.7e-8, 1.0, 100im)
    empty!(selected.preparations)
    empty!(selected.evaluations)
    copper=Material(kind = :conductor, rho = 1.7e-8)
    dielectric=Material(kind = :insulator, rho = Inf, eps_r = 2.3)
    design=build(CableDesign, "outer-only",
        terminal(:core, solid(copper, Disk(0.01)), insulation(dielectric; t = 0.002)))
    system=build(LineCableSystem, design, Pose2(0.0, 5.0); connections = Dict(:core=>1))
    problem=LineParametersProblem(system; frequencies = [50.0], earth_props = homogeneous(rho = 100.0))
    result=compute(problem, Formulation(internal_impedance = selected))
    @test all(isfinite, result.Z.values)
    @test keys(details(result).data.formulations.methods.internal_impedance.options)==(:outer,)
    @test length(selected.preparations)==1
    @test only(selected.evaluations)[2] === Val(:outer)
end

@testitem "Engine / each metal prepares shared surface state once per frequency" tags=[:unit] setup=[
    TestFixtures, FormulaFixtures] begin
    selected=FormulaFixtures.SurfaceLaw()
    problem=TestFixtures.line_parameters_problem(frequencies = [50.0, 100.0])
    result=compute(problem, Formulation(internal_impedance = selected))
    expected=sum(length(design.terminal_order)
    for design in problem.system.designs)*length(problem.frequencies)
    @test length(selected.preparations)==expected
    @test allunique(selected.evaluations)
    @test count(record -> record[2] === Val(:outer), selected.evaluations)==expected
    @test any(record -> record[2] === Val(:inner), selected.evaluations)
    @test any(record -> record[2] === Val(:transfer), selected.evaluations)
    @test all(isfinite, result.Z.values)
    empty!(selected.preparations)
    empty!(selected.evaluations)
    composite=compute(problem,
        Formulation(internal_impedance =
        (inner = selected, outer = selected, transfer = selected)))
    @test length(selected.preparations)==expected
    @test allunique(selected.evaluations)
    @test Z(composite)==Z(result) && Y(composite)==Y(result)
end

@testitem "Engine / first tubular primitive requests all surfaces without extra assembly terms" tags=[:unit] setup=[FormulaFixtures] begin
    M = FormulaFixtures
    copper = Material(kind=:conductor, rho=1.7e-8)
    dielectric = Material(kind=:insulator, rho=Inf, eps_r=2.3)
    design = build(CableDesign, "single-hollow-wall",
        Group(:wall, Region(:metal, Annulus(0.008, 0.01), copper)),
        Region(:insulation, Annulus(0.01, 0.012), dielectric))
    problem = CableConstantsProblem(design)
    first = M.SurfaceLaw()
    result = compute(problem, CableConstantsFormulation(internal_impedance=first))
    @test length(first.preparations) == 1
    @test Set(last.(first.evaluations)) == Set((Val(:inner), Val(:outer), Val(:transfer)))
    @test only(first.preparations)[1] == 0.008
    changed = M.SurfaceLaw(coefficients=(inner=100+20im, outer=2+1im, transfer=40+30im))
    other = compute(problem, CableConstantsFormulation(internal_impedance=changed))
    @test other.R == result.R
    @test other.L == result.L
    @test_throws ArgumentError compute(problem,
        CableConstantsFormulation(internal_impedance=(outer=M.SurfaceLaw(kinds=(:outer,)),)))
end

@testitem "Engine / unused tubular surface declarations allocate and evaluate nothing" tags=[:unit] setup=[FormulaFixtures] begin
    M = FormulaFixtures
    copper = Material(kind=:conductor, rho=1.7e-8)
    dielectric = Material(kind=:insulator, rho=Inf, eps_r=2.3)
    design = build(CableDesign, "solid-with-excess-recipe",
        Group(:core, Region(:metal, Disk(0.01), copper)),
        Region(:insulation, Annulus(0.01, 0.012), dielectric))
    unused = M.SpectralSurface()
    outer = M.SurfaceLaw(kinds=(:outer,))
    selected = (inner=unused, outer=outer, transfer=unused)
    local_result = compute(CableConstantsProblem(design),
        CableConstantsFormulation(internal_impedance=selected))
    @test isempty(unused.seen)
    system = build(LineCableSystem, design, Pose2(0.0, -1.0); connections=(core=1,))
    problem = LineParametersProblem(system; earth_props=homogeneous(rho=100.), frequencies=[50.])
    result = compute(problem, Formulation(internal_impedance=selected))
    @test isempty(unused.seen)
    @test all(isfinite, result.Z)
    @test all(isfinite, local_result.R)
end
