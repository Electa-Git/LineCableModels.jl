@testitem "Engine / indexed formula declarations and immutable source domains" tags=[:unit, :engine] begin
    const E = LineCableModels.Engine
    const EI = E.EarthImpedance
    const EA = E.EarthAdmittance
    const Expression = LineCableModels.Expression
    using LineCableModels.Commons: bindings
    air = E.EarthPair(1, 2, (10.0, 12.0), 1.0, (1, 1))
    soil = E.EarthPair(1, 2, (-1.0, -2.0), 1.0, (2, 2))
    mixed = E.EarthPair(1, 2, (10.0, -2.0), 1.0, (1, 2))
    self = E.EarthPair(1, 1, (-1.0, -1.0), 0.0, (2, 2); radius = 0.01)
    @test validate(self) === self
    @test self.radius == 0.01 && iszero(self.separation)
    @test_throws DomainError validate(E.EarthPair(1, 1, (-1.0, -1.0), 0.0, (2, 2)))
    @test_throws ArgumentError validate(E.EarthPair(
        1, 1, (-1.0, -1.0), 0.01, (2, 2); radius = 0.01))
    @test_throws ArgumentError validate(E.EarthPair(
        1, 1, (-1.0, -2.0), 0.0, (2, 2); radius = 0.01))
    @test_throws ArgumentError validate(E.EarthPair(1, 2, (1.0, 2.0), 1.0, (2, 2)))
    # The signatures of a formula's methods are its only declaration of the earth layers
    # it handles. A missing expression fails before the frequency loop.
    earth = homogeneous(rho = 100.0)
    for owner in (EI, EA)
        selected = owner.Formula(:default)
        @test only(bindings(selected, (air,))).expression.arguments[2:3] == (Val(1), Val(1))
        @test only(bindings(selected, (soil,))).expression.arguments[2:3] == (Val(2), Val(2))
        @test only(bindings(selected, (self,))).kind === :self
        @test only(bindings(selected, (soil,))).kind === :mutual
        @test only(bindings(selected, (mixed,))).kind === :mutual
        @test_throws ArgumentError bindings(selected, (E.EarthPair(
            1, 2, (-1.0, -2.0), 1.0, (2, 3)),))
        for pair in (air, soil, mixed, self)
            expression = Expression(selected, pair)
            @test validate(expression, earth) === expression
        end
        author=owner.Formula(:xue2018)
        @test_throws "the earth model has 2 layers and formula :xue2018 is defined up to layer 2; it has no expression for a mutual interaction from layer 1 to layer 2" validate(
            Expression(author, mixed), earth)
        @test only(bindings(selected, (air,))).expression isa Expression
        @test only(bindings(selected, (air,))).expression.selection === selected
    end
    for id in (:ametani2009, :lucca1994)
        selected = EI.Formula(id)
        @test only(bindings(selected, (mixed,))).kind === :mutual
        @test validate(Expression(selected, mixed), earth) isa Expression
        @test_throws ArgumentError validate(Expression(selected, air), earth)
        @test_throws ArgumentError validate(Expression(selected, soil), earth)
        @test_throws ArgumentError validate(Expression(selected, self), earth)
    end
    vertical=E.EarthPair(1, 2, (-1.0, -2.0), 0.0, (2, 2))
    @test only(bindings(EI.Formula(:saad1996), (vertical,))).kind === :mutual
    @test only(bindings(EI.Formula(:saad1996), (self,))).kind === :self
    @test only(bindings(EI.Formula(:wedepohl1973), (vertical,))).kind === :mutual
    @test only(bindings(EI.Formula(:wedepohl1973), (self,))).kind === :self
    @test_throws ArgumentError validate(Expression(EI.Formula(:pollaczek1926), air), earth)
    @test_throws "formula :carson1926 is defined up to layer 1" validate(
        Expression(EI.Formula(:carson1926), soil), earth)
end

@testitem "Engine / the plan reads an earth formula's media from its expression signatures" tags=[:unit, :engine] setup=[FormulaFixtures] begin
    const M=FormulaFixtures
    const E=M.E
    const EP=M.EP
    const EH=M.EH
    material=Material(kind = :conductor, rho = 1.7241e-8)
    design=build(CableDesign, "media-decision",
        Stack(Group(:phase, Region(:core, Disk(0.01), material))))
    two=homogeneous(rho = 100.0, eps_r = 10.0)
    three=build(EP.EarthModel, (
        EP.EarthLayer(100.0, 10.0, 1.0, 0.5), EP.EarthLayer(200.0, 20.0, 1.0)))
    four=build(EP.EarthModel, (EP.EarthLayer(100.0, 10.0, 1.0, 0.5),
        EP.EarthLayer(200.0, 20.0, 1.0, 0.5), EP.EarthLayer(300.0, 30.0, 1.0)))
    execution=computation_options(LineCableModelsCoaxial, ComputationOptions((;)))
    function plan(heights, earth, formulation)
        system=build(LineCableSystem, fill(design, length(heights)),
            [(Float64(i), h) for (i, h) in enumerate(heights)];
            connections = [(phase = i,) for i in eachindex(heights)])
        problem=LineParametersProblem(system; earth_props = earth, frequencies = [50.0])
        blueprints=E.CableBlueprint{Float64}[E.flatten(LineCableModelsCoaxial(), d, Float64)
                                             for d in problem.system.designs]
        workspace=E.LineParametersWorkspace(problem, formulation, execution, blueprints)
        return workspace.plan.earth_calculations, workspace.buffers.earth_materials
    end
    overhead, buried, deep=(10.0, 12.0), (-0.25, -0.3), (-0.25, -1.5)
    # A formula whose expressions stop at layer 2 sees a two-layer earth as it is.
    calculations, materials=plan(overhead, two, Formulation())
    @test only(calculations).earth isa EP.EarthModel
    @test size(only(materials).rho, 1) == 2 && only(materials).thickness === nothing
    # On more layers it consumes the `:default` reduction unless the formula specifies one.
    calculations, materials=plan(overhead, three, Formulation())
    @test only(calculations).earth isa EH.AfterFD
    @test formula_id(EH.rule(only(calculations).earth)) === :bottommost
    @test size(only(materials).rho, 1) == 2 && only(materials).thickness === nothing
    explicit=formula(:default; equivalent_earth = formula(:default; order = :before))
    calculations, _=plan(overhead, two,
        Formulation(earth_impedance = explicit, earth_admittance = explicit))
    @test only(calculations).earth isa EH.BeforeFD
    # A formula with an expression above layer 2 consumes the layered earth.
    layered=Formulation(earth_impedance = M.selection(M.EI), earth_admittance = M.selection(M.EA))
    calculations, materials=plan(deep, three, layered)
    @test all(calculation -> calculation.earth isa EP.EarthModel, calculations)
    @test all(buffer -> size(buffer.rho, 1) == 3 && buffer.thickness == [Inf, 0.5, Inf],
        materials)
    _, materials=plan(overhead, two, layered)
    @test all(buffer -> size(buffer.rho, 1) == 2 && buffer.thickness === nothing, materials)
    # Beyond its highest layer, a layered formula has no expression.
    @test_throws "the earth model has 4 layers and formula :LayerImpedance is defined up to layer 3" plan(
        (-0.25, -0.75, -1.5), four, layered)
    # The plan rejects a missing expression or an unadmitted reduction before the frequency loop.
    @test_throws "the earth model has 2 layers and formula :saad1996 is defined up to layer 2; it has no expression for a self interaction from layer 1 to layer 1" plan(
        overhead, two, Formulation(earth_impedance = formula(:saad1996)))
    @test_throws "does not admit equivalent-earth reduction :bottommost" plan(
        buried, three, Formulation(earth_impedance = formula(:saad1996)))
end

@testitem "Engine / a layered earth formula is never reduced without an explicit reduction" tags=[:unit, :engine] setup=[FormulaFixtures] begin
    const M=FormulaFixtures
    const EP=M.EP
    const EH=M.EH
    material=Material(kind = :conductor, rho = 1.7241e-8)
    design=build(CableDesign, "layered-safety",
        Stack(Group(:phase, Region(:core, Disk(0.01), material))))
    # Conductors in air, in the first earth layer and in the second.
    system=build(LineCableSystem, fill(design, 3), [(1.0, 2.0), (2.0, -0.25), (3.0, -1.5)];
        connections = [(phase = i,) for i in 1:3])
    three=build(EP.EarthModel, (
        EP.EarthLayer(100.0, 10.0, 1.0, 0.5), EP.EarthLayer(200.0, 20.0, 1.0)))
    problem=LineParametersProblem(system; earth_props = three, frequencies = [50.0])
    recorded(result, slot)=getproperty(details(result).data.formulations.methods, slot).equivalent_earth
    requested(result, slot)=get(getproperty(details(result).data.formulations.requested, slot),
        :equivalent_earth, nothing)
    potential=M.selection(M.EA)
    # The layered fixture sees all three media, as do a formula with typed runtime arguments
    # and a formula with generic layers. The record of each keeps `equivalent_earth` at
    # `nothing`.
    typed, generic=M.TypedLayerImpedance(), M.GenericLayerImpedance()
    empty!(M.calls)
    for impedance in (M.selection(M.EI), typed, generic)
        result=compute(problem, Formulation(earth_impedance = impedance, earth_admittance = potential))
        @test all(isfinite, result.Z.values) && all(isfinite, result.Y.values)
        @test recorded(result, :earth_impedance) === nothing
        @test recorded(result, :earth_admittance) === nothing
    end
    @test !isempty(M.calls) && all(call -> length(call[7]) == 3, M.calls)
    @test !isempty(typed.media) && all(==(3), typed.media)
    @test !isempty(generic.media) && all(==(3), generic.media)
    # A homogeneous formula on the same earth consumes the `:default` reduction and records
    # it as an explicit one would. The requested selection stays as declared.
    unspecified=compute(problem, Formulation())
    explicit=formula(:default; equivalent_earth = formula(:default))
    given=compute(problem, Formulation(earth_impedance = explicit, earth_admittance = explicit))
    @test unspecified.Z.values == given.Z.values && unspecified.Y.values == given.Y.values
    for slot in (:earth_impedance, :earth_admittance)
        @test recorded(unspecified, slot).order === :after
        @test recorded(unspecified, slot).rule.identifier === :bottommost
        @test recorded(unspecified, slot) == recorded(given, slot)
        @test requested(unspecified, slot) === nothing
    end
    # An explicit reduction on a layered formula applies, and the record shows it.
    empty!(M.calls)
    reduced=M.selection(M.EI; equivalent_earth = EH.AfterFD(:default))
    result=compute(problem, Formulation(earth_impedance = reduced, earth_admittance = potential))
    @test recorded(result, :earth_impedance).rule.identifier === :bottommost
    @test recorded(result, :earth_admittance) === nothing
    impedance_calls=filter(call -> call[1] === :EarthImpedance, M.calls)
    @test !isempty(impedance_calls) && all(call -> length(call[7]) == 2, impedance_calls)
end

@testitem "Engine / common earth functors retain evaluated material values" tags=[:unit, :engine] begin
    const E=LineCableModels.Engine
    μ0, ε0=4pi*1e-7, 8.8541878128e-12
    rho, epsilon, mu=[Inf, 100.0], [ε0, 10ε0], [μ0, μ0]
    jω=complex(0.0, 2pi*50)
    pair=E.EarthPair(1, 2, (-1.0, -2.0), 0.75, (2, 2))
    for owner in (E.EarthImpedance, E.EarthAdmittance)
        selected=owner.Formula(:default)
        functor=selected(rho, epsilon, mu, jω, pair)
        @test functor.binding.expression.selection === selected
        @test functor.state.mu === mu
        @test functor.state.epsilon === epsilon
        @test functor.state.sigma == Tuple(inv.(rho))
        @test_throws ArgumentError owner.Formula(:default; parameters = (unknown = 1,))
        @test_throws DimensionMismatch selected(
            [Inf, 100.0, 999.0], [ε0, 10ε0, 20ε0], [μ0, μ0, μ0], jω, pair)
    end
    @test LineCableModels.ComputationOptions === LineCableModels.Commons.ComputationOptions
    @test !(ComputationOptions() isa NamedTuple)
end

@testitem "Engine / earth field formulas consume material values without constitutive dispatch" tags=[:unit, :engine] begin
    const E = LineCableModels.Engine
    rho = [1e8, 100.0]
    epsilon = 8.8541878128e-12 .* [1.5, 10.0]
    mu = 4pi * 1e-7 .* [1.25, 2.5]
    s = 100pi * im
    for (owner, ids, heights, layers) in (
        (E.EarthImpedance, (:carson1926, :gary1976, :wise1934), (1.0, 1.5), (1, 1)),
        (E.EarthImpedance,
            (:pollaczek1926, :saad1996, :wedepohl1973, :xue2018, :unified),
            (-1.0, -1.5), (2, 2)),
        (E.EarthImpedance, (:ametani2009, :lucca1994), (1.0, -1.5), (1, 2)),
        (E.EarthAdmittance, (:wise1948,), (1.0, 1.5), (1, 1)),
        (E.EarthAdmittance, (:pollaczek1926, :xue2018, :unified), (-1.0, -1.5), (2, 2)))
        pair = E.EarthPair(1, 2, heights, 0.4, layers)
        for id in ids
            selected = owner.Formula(id)
            functor = selected(rho, epsilon, mu, s, pair)
            @test functor.state.rho === rho
            @test functor.state.mu === mu
            @test functor.state.epsilon === epsilon
            @test functor.state.sigma == Tuple(inv.(rho))
        end
    end
    @test mu == 4pi * 1e-7 .* [1.25, 2.5]
    @test epsilon == 8.8541878128e-12 .* [1.5, 10.0]
end

@testitem "Engine / indexed finite-layer methods determine admission and ordered assembly" tags=[:unit, :parametric, :slow] setup=[FormulaFixtures] begin
    using LinearAlgebra
    const M=FormulaFixtures
    const E=M.E
    const EI=M.EI
    const EA=M.EA
    const EP=M.EP
    using LineCableModels.Commons: bindings
    empty!(M.calls)
    inventories=(EI.formulas(EI.Formula), EA.formulas(EA.Formula))
    selected=M.selection(EI; options = (integration = (method = :quad, options = (;)),))
    heights=(2.0, -0.25, -1.5)
    pairs=[E.EarthPair(t, s, (heights[s], heights[t]), s==t ? 0.0 : 1.0,
               (s, t); radius = s==t ? 0.01 : nothing) for s in 1:3 for t in 1:3]
    bound=bindings(selected, pairs)
    @test isempty(M.calls) # Structural preflight does not execute kernels.
    @test count(case -> haskey(case.options.data, :integration), bound) == 1
    @test bound[1].options.data.integration.method === Val(:quad)
    @test_throws ArgumentError bindings(selected, pairs[2:end]) # Unconsumed integral controls.
    absent=E.EarthPair(1, 2, (-2.0, -3.0), 1.0, (3, 4))
    four=build(EP.EarthModel, (EP.EarthLayer(100.0, 10.0, 1.0, 0.5),
        EP.EarthLayer(200.0, 20.0, 1.0, 0.5), EP.EarthLayer(300.0, 30.0, 1.0)))
    message="the earth model has 4 layers and formula :LayerImpedance is defined up to layer 3; it has no expression for a mutual interaction from layer 3 to layer 4"
    @test_throws message E.validate(LineCableModels.Expression(selected, absent), four)
    @test isempty(M.calls)

    material=Material(kind = :conductor, rho = 1.7241e-8)
    design=build(CableDesign, "indexed-fixture",
        Stack(Group(:phase, Region(:core, Disk(0.01), material))))
    system=build(LineCableSystem, fill(design, 3),
        [(Float64(i), h) for (i, h) in enumerate(heights)];
        connections = [(phase = i,) for i in 1:3])
    earth=build(EP.EarthModel, (
        EP.EarthLayer(100.0, 10.0, 1.0, 0.5), EP.EarthLayer(200.0, 20.0, 1.0)))
    problem=LineParametersProblem(system; earth_props = earth, frequencies = [50.0, 500.0])
    fd=M.DispersiveEarth(exponent = 0)
    fd_calls=fd.seen
    formulation=Formulation(earth_impedance = selected, earth_admittance = M.selection(EA),
        earth_properties = fd,
        options = (ideal_transposition = false,))
    result=compute(problem, formulation; options = (trace = true,))
    @test length(fd_calls) == 4
    execution=computation_options(LineCableModelsCoaxial, ComputationOptions((;)))
    T=eltype(problem)
    blueprints=E.CableBlueprint{T}[E.flatten(LineCableModelsCoaxial(), d, T)
                                   for d in problem.system.designs]
    workspace=E.LineParametersWorkspace(problem, formulation, execution, blueprints)
    @test length(fd_calls) == 4 # Allocation did not reevaluate material laws.
    @test workspace.buffers.quadrature !== nothing
    @test length(M.calls) == 36
    for (_, _, row, column, layers, geometry, rho, physical) in M.calls
        @test layers == (column, row)
        @test geometry == (heights[column], heights[row])
        @test rho == [Inf, 50.0, 100.0]
        @test physical.layers == layers && physical.heights == geometry
    end
    for k in eachindex(problem.frequencies), t in 1:3, s in 1:3
        # Declared manufactured rule, evaluated independently of the probe method.
        coefficient=11s+17t+3t+5s+problem.frequencies[k]/100+(s==t ? 101 : 0)
        @test details(result).data.trace.Zg[t, s, k] ≈ coefficient*(1e-4+1e-3im) rtol=3e-6
        @test details(result).data.trace.Pg[t, s, k] ≈ coefficient*1e9
    end
    @test result.Z.values[1, 2, 1] ≈ (11*2+17+3+5*2+0.5)*(1e-4+1e-3im)
    @test result.Z.values[2, 1, 1] ≈ (11+17*2+3*2+5+0.5)*(1e-4+1e-3im)
    @test result.Y.values[:, :, 1] ≈ 2pi * 50im * inv(details(result).data.trace.P[:, :, 1])
    # A full-layer impedance and a separately reduced homogeneous potential
    # consume different material inventories in the same solve.
    buried=build(LineCableSystem, fill(design, 3),
        [(1.0, -0.1), (2.0, -0.25), (3.0, -1.5)];
        connections = [(phase = i,) for i in 1:3])
    buried_problem=LineParametersProblem(buried; earth_props = earth, frequencies = [50.0])
    for order in (:before, :after)
        empty!(M.calls)
        empty!(fd_calls)
        hybrid=Formulation(earth_impedance = M.selection(EI),
            earth_admittance = formula(:unified; equivalent_earth = formula(:default; order)),
            earth_properties = fd,
            options = (ideal_transposition = false,))
        value=compute(buried_problem, hybrid)
        @test all(isfinite, value.Z.values) && all(isfinite, value.Y.values)
        @test length(M.calls) == 9
        @test all(record -> length(record[7]) == 3, M.calls)
        @test Set(record[5] for record in M.calls) == Set(((2, 2), (2, 3), (3, 2), (3, 3)))
        @test length(fd_calls) == (order === :before ? 11 : 2)
        @test details(value).data.formulations.methods.earth_impedance.equivalent_earth ===
              nothing
        @test details(value).data.formulations.methods.earth_admittance.equivalent_earth.order ===
              order
    end
    @test (EI.formulas(EI.Formula), EA.formulas(EA.Formula)) === inventories
end

@testitem "Engine / unrelated coupled formula calculates in the main workspace without quadrature" tags=[:unit, :parametric] setup=[
    TestFixtures, FormulaFixtures] begin
    const E=LineCableModels.Engine
    const M=FormulaFixtures
    problem=TestFixtures.three_bare_wires_problem(heights = (-1.0, -1.0, -1.0),
        frequencies = [50.0, 50.0, 500.0])
    selected=Formulation(earth_impedance = M.CoupledImpedance(),
        earth_admittance = M.selection(M.EA),
        options = (
            reduce_bundle = false, kron_reduction = false, ideal_transposition = false))
    execution=computation_options(LineCableModelsCoaxial, ComputationOptions(trace = true))
    blueprints=[E.flatten(LineCableModelsCoaxial(), d, Float64)
                for d in problem.system.designs]
    workspace=E.LineParametersWorkspace(problem, selected, execution, blueprints)
    @test workspace.buffers.coupled isa Matrix{ComplexF64}
    @test !haskey(workspace.buffers, :quadrature)
    @test isempty(selected.methods.earth_impedance.solves)
    @test (@inferred E._solve!(workspace, selected)) === workspace
    @test length(selected.methods.earth_impedance.solves)==3
    @test isempty(workspace.trace.integrals)
    for k in 1:3, p in 1:3, q in 1:3
        @test workspace.trace.Zg[p, q, k]≈2π * im * problem.frequencies[k] * 1e-6 *
                                            (6+p+2q+(p==q ? 3 : 0))
    end
    E._solve!(workspace, selected)
    @test length(selected.methods.earth_impedance.solves)==6
    other=E.LineParametersWorkspace(problem, selected, execution, blueprints)
    @test other.buffers.coupled !== workspace.buffers.coupled
    @test compute(problem, selected).Z.values≈workspace.buffers.Zout
end

@testitem "Engine / one workspace overwrites its buffers independently of layout" tags=[:unit, :parametric] setup=[TestFixtures] begin
    const E=LineCableModels.Engine
    selection=Formulation(options = (reduce_bundle = false, kron_reduction = false,
        ideal_transposition = false))
    execution=computation_options(LineCableModelsCoaxial, ComputationOptions(trace = true))
    layouts=((1.0, 1.0, 1.0), (-1.0, -1.0, -1.0), (1.0, -1.0, -1.0),
        (-1.0, 1.0, -1.0), (-1.0, -1.0, 1.0))
    problems=map(
        heights->TestFixtures.three_bare_wires_problem(;
            heights, frequencies = [50.0, 50.0, 500.0]),
        layouts)
    workspaces=map(problems) do problem
        blueprints=only(E.flatten(LineCableModelsCoaxial(), problem.system.designs,
            Float64, [selection]))
        E.LineParametersWorkspace(problem, selection, execution, blueprints)
    end
    @test all(w -> typeof(w) === typeof(first(workspaces)), workspaces)
    two_wire_problem=TestFixtures.line_parameters_problem(TestFixtures.two_wire_system();
        frequencies = [50.0])
    two_wire_blueprints=only(E.flatten(LineCableModelsCoaxial(),
        two_wire_problem.system.designs, Float64, [selection]))
    two_wire_workspace=E.LineParametersWorkspace(two_wire_problem, selection,
        execution, two_wire_blueprints)
    @test typeof(two_wire_workspace) === typeof(first(workspaces))
    for (w, problem) in zip(workspaces, problems)
        @test length(w.buffers.earth_materials) == 1
        calculation=only(w.plan.earth_calculations)
        @test calculation.impedance_indices == calculation.potential_indices
        # Poison only numerical buffers. Identity matrices and geometry are inputs.
        for array in (w.buffers.Zearth, w.buffers.Pearth, w.buffers.Zprimitive,
            w.buffers.Pprimitive, w.buffers.rho_cond, w.buffers.dielectric_admittivity,
            w.buffers.earth.evaluated...)
            fill!(array, NaN)
        end
        for array in (
            w.buffers.axial_field, w.buffers.source_potential, w.buffers.current_map,
            w.buffers.enclosed_impedance, w.buffers.enclosed_potential, w.buffers.current_factor,
            w.buffers.current_rhs, w.buffers.radial_argument, w.buffers.source_logscale,
            w.buffers.circumference_average, w.buffers.radial_current, w.buffers.earth_spectrum...)
            fill!(array, NaN)
        end
        for materials in w.buffers.earth_materials
            foreach(a->fill!(a, NaN), (materials.rho, materials.epsilon, materials.mu))
        end
        E._solve!(w, selection)
        first_result=E._finish(w, problem, selection, execution.data.output_basis)
        Z, Y=copy(first_result.Z.values), copy(first_result.Y.values)
        @test all(isfinite, Z) && all(isfinite, Y)
        @test Z[:, :, 1] == Z[:, :, 2] && Y[:, :, 1] == Y[:, :, 2]
        # Γ=0: one Z and one voltage integral per ordered interaction. A second
        # solve for the potential consumer would double this observed work.
        @test length(w.trace.integrals) == 2*3^2*3
        fill!(w.buffers.Zearth, NaN)
        fill!(w.buffers.Pearth, NaN)
        E._solve!(w, selection)
        @test w.buffers.Zout == Z && w.buffers.Yout == Y
        @test first_result.Z.values == Z && first_result.Y.values == Y
        @test first_result.Z.values !== w.buffers.Zout
        @test details(first_result).data.trace.Zg !== w.trace.Zg
        @test details(first_result).data.trace.integrals !== w.trace.integrals
        @test length(w.trace.integrals) == 2*3^2*3 # Trace belongs to this solve only.
        # Public scans remain sorted. Revisit the actual material-earth stages
        # in reverse index order to expose any prior-frequency readiness state.
        for frequency in (3, 2, 1)
            E.materials!(w, selection, frequency)
            E.earth!(w, frequency)
            @test w.buffers.Zearth == w.trace.Zg[:, :, frequency]
            @test w.buffers.Pearth == w.trace.Pg[:, :, frequency]
        end
    end
    first, second=workspaces[1:2]
    for field in (:Zearth, :Pearth, :Zprimitive, :Pprimitive)
        @test getproperty(first.buffers, field) !== getproperty(second.buffers, field)
    end
    @test first.buffers.axial_field !== second.buffers.axial_field
    @test first.buffers.quadrature.segments !== second.buffers.quadrature.segments
end

@testitem "Engine / coupled material admission precedes integration and publication" tags=[:unit, :parametric] setup=[TestFixtures] begin
    const E=LineCableModels.Engine
    problem=TestFixtures.three_bare_wires_problem(frequencies = [50.0])
    selected=Formulation()
    execution=computation_options(LineCableModelsCoaxial, ComputationOptions(trace = true))
    blueprints=only(E.flatten(LineCableModelsCoaxial(), problem.system.designs, Float64, [selected]))
    w=E.LineParametersWorkspace(problem, selected, execution, blueprints)
    E.materials!(w, selected)
    E.materials!(w, selected, 1)
    fill!(w.buffers.axial_field, NaN)
    materials=only(w.buffers.earth_materials)
    materials.rho[2, 1]=-1
    @test_throws DomainError E.earth!(w, 1)
    @test isempty(w.trace.integrals)
    @test all(isnan, w.buffers.axial_field)
end

@testitem "Engine / distinct coupled controls perform separate complete calculations" tags=[:unit, :parametric, :slow] setup=[TestFixtures] begin
    problem=TestFixtures.three_bare_wires_problem(heights = (1.0, -1.0, -1.0), frequencies = [50.0])
    z=formula(:unified; options = (integration = (options = (rtol = 1e-8,),),))
    p=formula(:unified; options = (integration = (options = (rtol = 1e-9,),),))
    selected=Formulation(earth_impedance = z, earth_admittance = p,
        options = (
            reduce_bundle = false, kron_reduction = false, ideal_transposition = false))
    actual=compute(problem, selected; options = (trace = true,))
    trace=details(actual).data.trace
    @test length(trace.integrals) == 2*2*3^2
    for (definition, field) in ((z, :Z), (p, :Y))
        reference=compute(problem,
            Formulation(earth_impedance = definition,
                earth_admittance = definition, options = (reduce_bundle = false,
                    kron_reduction = false, ideal_transposition = false)))
        @test getproperty(actual, field).values ≈ getproperty(reference, field).values rtol=1e-12
    end
end

@testitem "Engine / fixed material evaluation follows allocation and repeats per solve" tags=[:unit, :parametric] setup=[
    TestFixtures, FormulaFixtures] begin
    const E=LineCableModels.Engine
    const M=FormulaFixtures
    td, fd=M.ScaledResistivity(), M.DispersiveEarth(exponent = 0)
    problem=TestFixtures.three_bare_wires_problem(frequencies = [50.0, 500.0])
    selected=Formulation(temperature_dependence = td, earth_properties = fd)
    execution=computation_options(LineCableModelsCoaxial, ComputationOptions())
    blueprints=only(E.flatten(LineCableModelsCoaxial(), problem.system.designs, Float64, [selected]))
    w=E.LineParametersWorkspace(problem, selected, execution, blueprints)
    @test isempty(td.seen) && isempty(fd.seen)
    E._solve!(w, selected)
    @test length(td.seen) == length(w.input.cable.conductor_materials)
    @test length(fd.seen) == length(problem.frequencies)
    @test all(record -> last(record) === w, td.seen)
    @test all(record -> last(record) === w, fd.seen)
    E._solve!(w, selected)
    @test length(td.seen) == 2length(w.input.cable.conductor_materials)
    @test length(fd.seen) == 2length(problem.frequencies)
end

@testitem "Engine / material laws and numerical options do not change physical routing" tags=[:unit, :parametric] setup=[
    TestFixtures, FormulaFixtures] begin
    const E=LineCableModels.Engine
    base=TestFixtures.line_parameters_problem(frequencies = [50.0, 1e5])
    fd=FormulaFixtures.DispersiveEarth(exponent = 0)
    @test_throws ArgumentError compute(base,
        Formulation(earth_impedance = :carson1926, earth_properties = fd))
    # Unsupported equation fails when invoked, after legitimate material stages.
    empty!(fd.seen)
    selected=Formulation(earth_properties = fd)
    changed=compute(base, selected)
    material=base.earth_props.layers[2]
    equivalent=LineParametersProblem(base.system; temperature = base.temperature,
        earth_props = homogeneous(rho = material.rho/2, eps_r = material.eps_r, mu_r = material.mu_r),
        frequencies = base.frequencies)
    expected=compute(equivalent)
    @test changed.Z.values ≈ expected.Z.values rtol=1e-10
    @test changed.Y.values ≈ expected.Y.values rtol=1e-10
    @test Set(record[2] for record in fd.seen)==Set(base.frequencies)
    @test details(changed).data.formulations.methods.earth_properties.identifier ===
          :DispersiveEarth
    reference=compute(base)
    for rtol in (1e-7, 1e-9)
        result=compute(base,
            Formulation(
                earth_impedance = formula(
                    :default; options = (integration = (
                        method = :quad, options = (; rtol)),)),
                earth_admittance = formula(
                    :default; options = (integration = (
                        method = :quad, options = (; rtol)),))))
        @test result.Z.values ≈ reference.Z.values rtol=3e-6
        @test result.Y.values ≈ reference.Y.values rtol=3e-6
        controls=details(result).data.formulations.methods.earth_impedance.options.integration
        @test controls.method === :quad
        @test controls.options.rtol === rtol
    end
end

@testitem "Engine / consumer earth state preserves scalar types and uncertainty" tags=[:unit, :measurements] setup=[FormulaFixtures] begin
    using Measurements
    const E=LineCableModels.Engine
    for T in (Float32, Float64, BigFloat)
        rho=T[Inf, 100]
        epsilon=T(8.8541878128e-12) .* T[1, 10]
        mu=fill(T(4pi*1e-7), 2)
        pair=E.EarthPair(1, 2, (-one(T), -T(2)), T(0.75), (2, 2))
        s=complex(zero(T), T(2)*T(pi)*T(50))
        for owner in (E.EarthImpedance, E.EarthAdmittance), method in (:quad,)

            selected=FormulaFixtures.selection(owner; scale = one(T))
            functor=selected(rho, epsilon, mu, s, pair)
            @test functor.state.jω isa Complex{T}
            result=functor()
            @test result isa Complex{T}
            @test isfinite(result)
        end
    end
    rho=measurement.([Inf, 100.0], [0.0, 1.0])
    epsilon=measurement.(8.8541878128e-12 .* [1, 10], 0.0)
    mu=measurement.(fill(4pi*1e-7, 2), 0.0)
    s=complex(measurement(0.0, 0.0), measurement(2pi*50, 0.0))
    pair=E.EarthPair(1, 2, (-1.0, -2.0), 0.75, (2, 2))
    for owner in (E.EarthImpedance, E.EarthAdmittance)
        value=FormulaFixtures.selection(owner; scale = rho[2])(
            rho, epsilon, mu, s, pair)()
        @test value isa Complex{Measurement{Float64}}
        @test isfinite(value)
        @test uncertainty(real(value)) > 0
        @test uncertainty(abs(value)) > 0
    end
end
