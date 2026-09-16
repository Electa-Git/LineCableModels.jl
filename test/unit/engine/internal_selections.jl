@testitem "Engine / internal selections preserve native transfer dispatch and scalar assembly" tags=[:unit] setup=[TestFixtures,FormulaContractModels] begin
    const II=LineCableModels.Engine.InternalImpedance
    const M=FormulaContractModels
    const IO=LineCableModels.ImportExport
    args=(0.008,0.01,1.7241e-8,1.,100.0im)
    scalar=II.Formula(:default)
    same=Formulation(II.Formula,(transfer=:default,inner=:default,outer=:default))
    @test keys(same)==(:inner,:outer,:transfer)
    reference=@inferred II.surface_impedances(scalar,Val((:outer,:transfer,:inner)),args...)
    @test (@inferred II.surface_impedances(same,Val((:outer,:transfer,:inner)),args...))==reference
    @test II.surface_impedances(same,args...)==II.surface_impedances(scalar,args...)
    coefficients=(inner=11+12im,outer=21+22im,transfer=31+32im)
    customized=(inner=M.SurfaceLaw(kinds=(:inner,);coefficients),
        outer=M.SurfaceLaw(kinds=(:outer,);coefficients),
        transfer=M.SurfaceLaw(kinds=(:transfer,);coefficients))
    selections=Formulation(II.Formula,customized)
    @test selections === customized
    @test II.surface_impedances(selections,args...)==coefficients
    @test all(length(leaf.evaluations)==1 for leaf in selections)
    empty!(selections.transfer.evaluations)
    changed=II.surface_impedances(merge(same,(transfer=selections.transfer,)),args...)
    @test changed.inner==reference.inner && changed.outer==reference.outer
    @test changed.transfer==31+32im
    @test only(selections.transfer.evaluations)[2] === Val(:transfer)
    @test_throws ArgumentError Formulation(II.Formula,(inner=:default,outer=:default,mutual=:default))
    @test_throws ArgumentError Formulation(II.Formula,(outer=:default,))
    @test_throws ArgumentError validate(merge(same,(outer=selections.inner,)),(:outer,))
    @test_throws ArgumentError validate(selections,(:outer,))
    @test validate(same,(:outer,)) === same
    @test keys(II.surface_impedances(same,Val((:outer,)),args...))==(:outer,)

    struct ObservedSurfaces{F,P,O,C} <: LineCableModels.Engine.InternalImpedanceFormulation
        base::F
        parameters::P
        options::O
        configured_options::C
        observations::Vector{Tuple}
    end
    observed=ObservedSurfaces(scalar,(;),scalar.options,(),Tuple[])
    function (leaf::ObservedSurfaces)(args...)
        prepared=leaf.base(args...)
        II.Functor(leaf,prepared.state,leaf.options)
    end
    function II.internal_impedance(leaf::ObservedSurfaces,kind::Union{Val{:inner},Val{:outer},Val{:transfer}},
            functor,workspace)
        push!(leaf.observations,(functor.state,workspace))
        II.internal_impedance(leaf.base,kind,functor,workspace)
    end
    args32=(0.003f0,0.005f0,2f-8,1f0,20Float32(pi)*im)
    workspace=Ref(:surface_workspace)
    delegated=II.surface_impedances(observed,args32...;workspace)
    @test delegated === II.surface_impedances(scalar,args32...)
    @test all(value->value isa ComplexF32,delegated)
    @test length(observed.observations)==3
    @test all(record->record[2] === workspace,observed.observations)
    @test all(record->record[1] === first(observed.observations)[1],observed.observations)
    state=first(observed.observations)[1]
    @test (state.r_in,state.r_ex,state.rho_c,state.mur_c,state.jω) === args32

    alternative=M.SurfaceLaw(kinds=(:transfer,),coefficients=(transfer=41+42im,))
    distinct=merge(same,(transfer=alternative,))
    @test validate(distinct,(:inner,:outer,:transfer)) === distinct
    values=II.surface_impedances(distinct,args...)
    @test values.inner==reference.inner && values.outer==reference.outer
    @test values.transfer==41+42im && length(alternative.preparations)==1
    @test_throws ArgumentError validate(merge(same,(inner=alternative,)),(:inner,:outer,:transfer))

    problem=TestFixtures.line_parameters_problem(frequencies=[50.,500.])
    original=compute(problem,Formulation())
    composed=compute(problem,Formulation(internal_impedance=same))
    @test Z(composed)==Z(original) && Y(composed)==Y(original)
    @test details(composed).data.formulations.effective.internal_impedance==
        (inner=:schelkunoff1934,outer=:schelkunoff1934,transfer=:schelkunoff1934)
    cable_problem=CableConstantsProblem(first(problem.system.designs);frequency=50.)
    @test compute(cable_problem,CableConstantsFormulation())==
        compute(cable_problem,CableConstantsFormulation(internal_impedance=same))
    custom=compute(problem,Formulation(internal_impedance=distinct))
    @test Z(custom)!=Z(original) && Y(custom)==Y(original)
    @test details(custom).data.formulations.effective.internal_impedance.transfer === :SurfaceLaw
    candidates=Formulation(internal_impedance=Grid((formula(:default),customized));combine=:zip)
    @test length(candidates)==2
    @test collect(candidates)[2].methods.internal_impedance.transfer === customized.transfer
    native=Formulation(internal_impedance=same)
    for source in (native,MonteCarlo(native),LinearError(native))
        saved=IO.deserialize_value(Val(:formulation),NamedTuple(source))
        @test description([source];quantity=R)==description([saved];quantity=R)
        @test occursin("internal Z(transfer)",only(description([source];quantity=R)))
        @test !occursin("internal Z",only(description([source];quantity=B)))
    end
end
