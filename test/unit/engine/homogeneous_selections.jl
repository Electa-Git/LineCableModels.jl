@testitem "Engine / homogeneous selections retain each indexed formula through assembly" tags=[:unit] setup=[FormulaContractModels] begin
    const E=LineCableModels.Engine
    # Manufactured potential coefficients isolate the selection-routing check.
    potential=FormulaContractModels.selection(E.EarthAdmittance; layers = 2:2)
    material=Material(kind = :conductor, rho = 1.7241e-8)
    dielectric=Material(kind = :insulator, rho = 1e14, eps_r = 2.3)
    design=build(CableDesign,
        "selection-contract",
        Stack(Group(:phase,
                Region(:core, Disk(0.01), material)),
            Region(:insulation, Annulus(0.01, 0.012), dielectric)))
    positions=[(0.0, 10.0), (1.0, 12.0), (2.0, -1.0), (3.0, -2.0)]
    system=build(LineCableSystem, fill(design, 4), positions;
        connections = [(phase = i,) for i in 1:4])
    problem=LineParametersProblem(system; earth_props = homogeneous(rho = 100.0), frequencies = [
        50.0, 500.0])
    choices=(
        air = formula(:Carson1926; options = (integration = (options = (rtol = 1e-6,),),)),
        earth = formula(:Pollaczek1926; options = (integration = (options = (rtol = 1e-8,),),)),
        mixed = formula(:Lucca1994))
    selected=Formulation(earth_impedance = choices, earth_admittance = potential,
        options = (ideal_transposition = false,))
    result=compute(problem, selected; options = (trace = true,))
    @test all(isfinite, result.Z.values) && all(isfinite, result.Y.values)
    @test result.Z.values ≈ permutedims(result.Z.values, (2, 1, 3))
    # The two same-medium subproblems recover the corresponding assembled blocks.
    for (indices, leaf) in ((1:2, choices.air), (3:4, choices.earth))
        subsystem=build(LineCableSystem, fill(design, 2), positions[indices];
            connections = [(phase = i,) for i in 1:2])
        subproblem=LineParametersProblem(subsystem; earth_props = problem.earth_props,
            frequencies = problem.frequencies)
        subresult=compute(subproblem,
            Formulation(earth_impedance = leaf,
                earth_admittance = potential, options = (ideal_transposition = false,));
            options = (trace = true,))
        @test result.Z.values[indices, indices, :] ≈ subresult.Z.values
        @test details(result).trace.Pg[indices, indices, :] ≈ details(subresult).trace.Pg
    end
    provenance=details(result).formulations
    @test provenance.requested.earth_impedance == NamedTuple(selected).requested.earth_impedance
    @test provenance.effective.earth_impedance ==
          map(record -> record.identifier,provenance.requested.earth_impedance)
    records=provenance.numerical.earth_impedance
    @test Set((record.formula, record.source, record.target) for record in records) ==
          Set(((:Carson1926, 1, 1), (:Pollaczek1926, 2, 2),
        (:Lucca1994, 1, 2), (:Lucca1994, 2, 1)))
    @test all(record -> isempty(record.options), filter(r -> r.formula === :Lucca1994, records))
    @test all(record -> record.options.integration.options.rtol == 1e-6,
        filter(r -> r.formula === :Carson1926, records))
    @test all(record -> record.options.integration.options.rtol == 1e-8,
        filter(r -> r.formula === :Pollaczek1926, records))
    @test_throws ArgumentError compute(problem, Formulation(earth_impedance = formula(:Lucca1994)))
    hybrid=compute(problem,
        Formulation(earth_impedance = merge(choices,
                (mixed = formula(:default),)),
            earth_admittance = potential,
            options = (ideal_transposition = false,));
        options = (trace = true,))
    complete=compute(problem,
        Formulation(earth_impedance = :default,
            earth_admittance = potential, options = (ideal_transposition = false,));
        options = (trace = true,))
    for p in 1:4, q in 1:4
        # A selected default block still solves the complete four-wire system.
        reference=sign(last(positions[p]))==sign(last(positions[q])) ? result : complete
        @test details(hybrid).trace.Zg[p, q, :]≈details(reference).trace.Zg[p, q, :] rtol=1e-10
    end

    # Different configurations of the same source must not be merged by identity tag.
    calls=Tuple[]
    air_hook=(functor, pair,
        workspace)->begin
        push!(calls, (:air, pair.layers, pair.row, pair.column))
        1.0+2.0im
    end
    earth_hook=(functor, pair,
        workspace)->begin
        push!(calls, (:earth, pair.layers, pair.row, pair.column))
        3.0+4.0im
    end
    mixed_hook=(functor, pair,
        workspace)->begin
        push!(calls, (:mixed, pair.layers, pair.row, pair.column))
        complex(pair.layers[1], pair.layers[2])
    end
    for (id,
        hook) in ((:default, air_hook), (:default, earth_hook), (:Lucca1994, mixed_hook))
        @eval LineCableModels.computation_options(
            ::LineCableModels.FormulaMethod{$(QuoteNode(id)),
                typeof(E.EarthImpedance.earth_impedance)},
            ::$(typeof(hook))) = (;)
    end
    modified=compute(problem,
        Formulation(
            earth_impedance = (
                air = formula(:default; hooks = (contribution = air_hook,)),
                earth = formula(:default; hooks = (contribution = earth_hook,)),
                mixed = formula(:Lucca1994; hooks = (contribution = mixed_hook,))),
            earth_admittance = potential, options = (ideal_transposition = false,));
        options = (trace = true,))
    @test length(calls) == 32
    @test all(call -> call[2] == (1, 1), filter(call -> call[1] === :air, calls))
    @test all(call -> call[2] == (2, 2), filter(call -> call[1] === :earth, calls))
    @test details(modified).trace.Zg[3, 1, 1] == 1 + 2im
    @test details(modified).trace.Zg[1, 3, 1] == 2 + 1im
    @test details(modified).formulations.modified.earth_impedance ==
          (air = true, earth = true, mixed = true)

    # Y has independent selections and hooks under exactly the same grammar.
    potential_calls=Symbol[]
    potential_hook=name->(
        functor, pair, workspace)->begin
        push!(potential_calls, name)
        coefficient=pair.row==pair.column ? 100 : 10
        coefficient+pair.layers[1]+0.1*pair.layers[2]
    end
    potential_choices=map((:air, :earth, :mixed)) do name
        hook=potential_hook(name)
        @eval LineCableModels.computation_options(
            ::LineCableModels.FormulaMethod{:ContractLayers,
                typeof(E.EarthAdmittance.earth_potential_coefficient)},
            ::$(typeof(hook))) = (;)
        FormulaContractModels.selection(E.EarthAdmittance;
            layers = 2:2, hooks = (contribution = hook,))
    end
    independent_y=compute(problem,
        Formulation(earth_impedance = choices,
            earth_admittance = NamedTuple{(:air, :earth, :mixed)}(potential_choices),
            options = (ideal_transposition = false,));
        options = (trace = true,))
    @test count(==(:air), potential_calls) == 8
    @test count(==(:earth), potential_calls) == 8
    @test count(==(:mixed), potential_calls) == 16
    @test details(independent_y).trace.Pg[3, 1, 1] == 11.2
    @test details(independent_y).trace.Pg[1, 3, 1] == 12.1
    @test independent_y.Z.values == result.Z.values
    @test details(independent_y).formulations.modified.earth_admittance ==
          (air = true, earth = true, mixed = true)
end

@testitem "Engine / scalar and homogeneous shorthand preserve numerical and model contracts" tags=[:unit] setup=[TestFixtures] begin
    const E=LineCableModels.Engine
    problem=TestFixtures.line_parameters_problem(frequencies = [50.0, 500.0])
    scalar=compute(problem, Formulation())
    same=(air = formula(:default), earth = formula(:default), mixed = formula(:default))
    shorthand=compute(problem, Formulation(earth_impedance = same, earth_admittance = same))
    @test shorthand.Z.values == scalar.Z.values
    @test shorthand.Y.values == scalar.Y.values
    @test_throws ArgumentError Formulation(earth_impedance = (air = formula(:default),))
    @test_throws ArgumentError Formulation(earth_admittance = merge(same, (other = formula(:default),)))
    reordered=Formulation(earth_impedance = (
        mixed = same.mixed, earth = same.earth, air = same.air))
    @test keys(reordered.methods.earth_impedance) === (:air, :earth, :mixed)
    @test keys(reordered.definitions.earth_impedance) === (:air, :earth, :mixed)
    @test_throws ArgumentError E.Formulation(reordered.methods.earth_impedance, Val(2), Val(3))
    model=build(E.EarthModel,
        (LineCableModels.Earth.EarthLayer(100.0, 10.0, 1.0, 0.5),
            LineCableModels.Earth.EarthLayer(200.0, 20.0, 1.0)))
    layered=LineParametersProblem(problem.system; earth_props = model, frequencies = [50.0])
    @test_throws ArgumentError compute(layered, Formulation(earth_impedance = same))
    air_system=build(LineCableSystem, problem.system.designs,
        [Pose2(i, 10.0) for i in eachindex(problem.system.designs)];
        connections = [Dict(:core=>i, :sheath=>0, :jacket=>0)
                       for i in eachindex(problem.system.designs)])
    overhead_layered=LineParametersProblem(air_system; earth_props = model, frequencies = [50.0])
    @test_throws ArgumentError compute(overhead_layered,
        Formulation(earth_impedance = same, earth_admittance = same))
    # Scalar EHEM remains explicit and can consume the full physical soil inventory.
    reduced=compute(layered,
        Formulation(
            earth_impedance = formula(:default; equivalent_earth = formula(:default)),
            earth_admittance = formula(:default; equivalent_earth = formula(:default))))
    @test all(isfinite, reduced.Z.values) && all(isfinite, reduced.Y.values)
end
