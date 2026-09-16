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
        air = formula(:carson1926; options = (integration = (options = (rtol = 1e-6,),),)),
        earth = formula(:pollaczek1926; options = (integration = (options = (rtol = 1e-8,),),)),
        mixed = formula(:lucca1994))
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
        # The manufactured potential deliberately identifies local row/column
        # indices. Rebuilding a subproblem renumbers them; it cannot preserve the
        # corresponding numerical block as a physical potential formula would.
        layer = last(positions[first(indices)]) > 0 ? 1 : 2
        for (row, p) in enumerate(indices), (column, q) in enumerate(indices),
                (k, frequency) in enumerate(problem.frequencies)
            coefficient = (i, j) -> 1e9 * (11layer + 17layer + 3i + 5j +
                frequency/100 + (i == j ? 101 : 0))
            @test details(result).data.trace.Pg[p, q, k] ≈ coefficient(p, q)
            @test details(subresult).data.trace.Pg[row, column, k] ≈ coefficient(row, column)
        end
    end
    provenance=details(result).data.formulations
    @test provenance.requested.earth_impedance == NamedTuple(selected).requested.earth_impedance
    @test provenance.effective.earth_impedance ==
          map(record -> record.identifier,provenance.requested.earth_impedance)
    records=provenance.numerical.earth_impedance
    @test Set((record.formula, record.source, record.target) for record in records) ==
          Set(((:carson1926, 1, 1), (:pollaczek1926, 2, 2),
        (:lucca1994, 1, 2), (:lucca1994, 2, 1)))
    @test all(record -> isempty(record.options), filter(r -> r.formula === :lucca1994, records))
    @test all(record -> record.options.integration.options.rtol == 1e-6,
        filter(r -> r.formula === :carson1926, records))
    @test all(record -> record.options.integration.options.rtol == 1e-8,
        filter(r -> r.formula === :pollaczek1926, records))
    @test_throws ArgumentError compute(problem, Formulation(earth_impedance = formula(:lucca1994)))
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
        @test details(hybrid).data.trace.Zg[p, q, :]≈details(reference).data.trace.Zg[p, q, :] rtol=1e-10
    end

    # Different parameterizations of one native type remain distinct selections.
    M=FormulaContractModels
    empty!(M.calls)
    native_choices=(air=M.selection(E.EarthImpedance;layers=2:2,scale=1.0),
        earth=M.selection(E.EarthImpedance;layers=2:2,scale=2.0),
        mixed=M.selection(E.EarthImpedance;layers=2:2,scale=3.0))
    custom=compute(problem,Formulation(earth_impedance=native_choices,
        earth_admittance=potential,options=(ideal_transposition=false,));options=(trace=true,))
    impedance_calls=filter(record->record[1] === :EarthImpedance,M.calls)
    @test length(impedance_calls)==32
    for k in eachindex(problem.frequencies), row in 1:4, column in 1:4
        s=positions[column][2]>0 ? 1 : 2
        t=positions[row][2]>0 ? 1 : 2
        scale=s==t ? Float64(s) : 3.0
        coefficient=11s+17t+3row+5column+problem.frequencies[k]/100+(row==column ? 101 : 0)
        @test details(custom).data.trace.Zg[row,column,k] ≈ scale*coefficient*(1e-4+1e-3im)
    end
    @test details(custom).data.formulations.effective.earth_impedance==
        (air=:LayerImpedance,earth=:LayerImpedance,mixed=:LayerImpedance)
    @test details(custom).data.formulations.requested.earth_impedance.earth.parameters.scale==2.0

    # Potential coefficients use exactly the same air/earth/mixed grammar.
    empty!(M.calls)
    potential_choices=(air=M.selection(E.EarthAdmittance;layers=2:2,scale=1.0),
        earth=M.selection(E.EarthAdmittance;layers=2:2,scale=2.0),
        mixed=M.selection(E.EarthAdmittance;layers=2:2,scale=3.0))
    independent_y=compute(problem,Formulation(earth_impedance=choices,
        earth_admittance=potential_choices,options=(ideal_transposition=false,));options=(trace=true,))
    potential_calls=filter(record->record[1] === :EarthAdmittance,M.calls)
    @test count(record->record[5]==(1,1),potential_calls)==8
    @test count(record->record[5]==(2,2),potential_calls)==8
    @test count(record->record[5][1]!=record[5][2],potential_calls)==16
    for row in 1:4, column in 1:4
        s=positions[column][2]>0 ? 1 : 2
        t=positions[row][2]>0 ? 1 : 2
        scale=s==t ? Float64(s) : 3.0
        @test details(independent_y).data.trace.Pg[row,column,:] ≈ scale .* details(result).data.trace.Pg[row,column,:]
    end
    @test independent_y.Z.values==result.Z.values
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
        connections = [Dict(:core=>i, :sheath=>0)
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
