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
    @test map(value -> value.identifier, details(custom).data.formulations.methods.earth_impedance)==
        (air=:LayerImpedance,earth=:LayerImpedance,mixed=:LayerImpedance)
    @test details(custom).data.formulations.requested.earth_impedance.earth.parameters.scale==2.0

    # Partial Unified publication must still use the complete physical system.
    # The other entries come from the explicit manufactured equation above.
    hybrid=compute(problem, Formulation(
        earth_impedance=merge(native_choices, (mixed=formula(:default),)),
        earth_admittance=potential, options=(ideal_transposition=false,));
        options=(trace=true,))
    complete=compute(problem, Formulation(earth_impedance=:default,
        earth_admittance=potential, options=(ideal_transposition=false,));
        options=(trace=true,))
    for p in 1:4, q in 1:4
        reference=sign(last(positions[p]))==sign(last(positions[q])) ? custom : complete
        @test details(hybrid).data.trace.Zg[p,q,:] ≈
              details(reference).data.trace.Zg[p,q,:] rtol=1e-10
    end

    # Potential coefficients use exactly the same air/earth/mixed grammar.
    empty!(M.calls)
    potential_choices=(air=M.selection(E.EarthAdmittance;layers=2:2,scale=1.0),
        earth=M.selection(E.EarthAdmittance;layers=2:2,scale=2.0),
        mixed=M.selection(E.EarthAdmittance;layers=2:2,scale=3.0))
    independent_y=compute(problem,Formulation(earth_impedance=native_choices,
        earth_admittance=potential_choices,options=(ideal_transposition=false,));options=(trace=true,))
    potential_calls=filter(record->record[1] === :EarthAdmittance,M.calls)
    @test count(record->record[5]==(1,1),potential_calls)==8
    @test count(record->record[5]==(2,2),potential_calls)==8
    @test count(record->record[5][1]!=record[5][2],potential_calls)==16
    for row in 1:4, column in 1:4
        s=positions[column][2]>0 ? 1 : 2
        t=positions[row][2]>0 ? 1 : 2
        scale=s==t ? Float64(s) : 3.0
        @test details(independent_y).data.trace.Pg[row,column,:] ≈ scale .* details(custom).data.trace.Pg[row,column,:]
    end
    @test independent_y.Z.values==custom.Z.values
end

@testitem "Engine / scalar and homogeneous shorthand preserve numerical and model contracts" tags=[:unit] setup=[TestFixtures] begin
    const E=LineCableModels.Engine
    problem=TestFixtures.line_parameters_problem(frequencies = [50.0, 500.0])
    scalar=compute(problem, Formulation())
    same=(air = formula(:default), earth = formula(:default), mixed = formula(:default))
    shorthand=compute(problem, Formulation(earth_impedance = same, earth_admittance = same))
    @test shorthand.Z.values == scalar.Z.values
    @test shorthand.Y.values == scalar.Y.values
    partial=Formulation(earth_impedance=(air=formula(:default),))
    @test keys(partial.methods.earth_impedance) === (:air,)
    @test_throws ArgumentError compute(problem, partial)
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
    @test_throws DimensionMismatch compute(overhead_layered,
        Formulation(earth_impedance = same, earth_admittance = same))
    # Scalar EHEM remains explicit and can consume the full physical soil inventory.
    reduced=compute(layered,
        Formulation(
            earth_impedance = formula(:default; equivalent_earth = formula(:default)),
            earth_admittance = formula(:default; equivalent_earth = formula(:default))))
    @test all(isfinite, reduced.Z.values) && all(isfinite, reduced.Y.values)
end
