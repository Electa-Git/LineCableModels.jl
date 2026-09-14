@testitem "Engine / internal selections preserve transfer providers and scalar assembly" tags=[:unit] setup=[TestFixtures] begin
    const II=LineCableModels.Engine.InternalImpedance
    const FM=LineCableModels.FormulaMethod
    const IO=LineCableModels.ImportExport
    calls=Symbol[]
    inner=(functor,workspace)->(push!(calls,:inner);11+12im)
    outer=(functor,workspace)->(push!(calls,:outer);21+22im)
    transfer=(functor,workspace)->(push!(calls,:transfer);31+32im)
    for (kind,provider) in ((:inner,inner),(:outer,outer),(:transfer,transfer))
        @eval LineCableModels.computation_options(
            ::FM{:default,typeof(II.internal_impedance),Tuple{Val{$(QuoteNode(kind))}}},
            ::$(typeof(provider)))=(;)
    end
    args=(0.008,0.01,1.7241e-8,1.,100.0im)
    scalar=II.Formula(:default)
    same=Formulation(II.Formula,(transfer=:default,inner=:default,outer=:default))
    @test keys(same)==(:inner,:outer,:transfer)
    reference=@inferred II.surface_impedances(scalar,Val((:outer,:transfer,:inner)),args...)
    @test (@inferred II.surface_impedances(same,Val((:outer,:transfer,:inner)),args...))==reference
    @test II.surface_impedances(same,args...)==II.surface_impedances(scalar,args...)
    customized=(inner=formula(:default;hooks=(inner=inner,)),
        outer=formula(:default;hooks=(outer=outer,)),
        transfer=formula(:default;hooks=(transfer=transfer,)))
    selections=Formulation(II.Formula,customized)
    @test II.surface_impedances(selections,args...)==(inner=11+12im,outer=21+22im,transfer=31+32im)
    @test calls==[:inner,:outer,:transfer]
    empty!(calls)
    modified=II.surface_impedances(merge(same,(transfer=selections.transfer,)),args...)
    @test modified.inner==reference.inner && modified.outer==reference.outer
    @test modified.transfer==31+32im && calls==[:transfer]
    @test_throws ArgumentError Formulation(II.Formula,(inner=:default,outer=:default,mutual=:default))
    @test_throws ArgumentError Formulation(II.Formula,(outer=:default,))
    @test_throws ArgumentError II.Formula(:default;hooks=(mutual=transfer,))
    @test_throws ArgumentError validate(merge(same,(outer=selections.inner,)),(:outer,))
    @test_throws ArgumentError validate(selections,(:outer,))
    @test validate(same,(:outer,))===same
    @test keys(II.surface_impedances(same,Val((:outer,)),args...))==(:outer,)

    # A different ID must dispatch to its own equation and prepare only that
    # provider. Hooks on a shared default are not a substitute for this gate.
    preparations=Ref(0)
    II.internal_impedance(::Val{:TestTransfer},::Val{:transfer},functor,workspace)=41+42im
    LineCableModels.computation_options(::FM{:TestTransfer,typeof(II.internal_impedance),Tuple{Val{:transfer}}})=(;)
    @eval function (leaf::II.Formula{:TestTransfer})(r_in,r_ex,rho,mu_r,jω)
        $preparations[]+=1
        II.Functor{:TestTransfer,typeof(leaf.binding),typeof(leaf.hooks),Nothing,typeof(leaf.options)}(
            leaf.binding,leaf.hooks,nothing,leaf.options)
    end
    binding=(transfer=FM(Val(:TestTransfer),II.internal_impedance,Val(:transfer)),)
    controls=(transfer=(;),)
    alternative=II.Formula{:TestTransfer,typeof(binding),NamedTuple{()},NamedTuple{()},typeof(controls),Tuple{}}(
        binding,(;),(;),controls,())
    distinct=merge(same,(transfer=alternative,))
    @test validate(distinct,(:inner,:outer,:transfer))===distinct
    coefficients=II.surface_impedances(distinct,args...)
    @test coefficients.inner==reference.inner && coefficients.outer==reference.outer
    @test coefficients.transfer==41+42im && preparations[]==1
    @test_throws ArgumentError validate(merge(same,(inner=alternative,)),(:inner,:outer,:transfer))

    # New named selections must reach both calculation entry points; equal
    # surface formulas must not alter their current-basis assembly or Y.
    problem=TestFixtures.line_parameters_problem(frequencies=[50.,500.])
    original=compute(problem,Formulation())
    composed=compute(problem,Formulation(internal_impedance=same))
    @test Z(composed)==Z(original) && Y(composed)==Y(original)
    @test details(composed).formulations.effective.internal_impedance==
        (inner=:default,outer=:default,transfer=:default)
    cable_problem=CableConstantsProblem(first(problem.system.designs);frequency=50.)
    a=compute(cable_problem,CableConstantsFormulation())
    b=compute(cable_problem,CableConstantsFormulation(internal_impedance=same))
    @test a==b
    changed=compute(problem,Formulation(internal_impedance=distinct))
    @test Z(changed)!=Z(original) && Y(changed)==Y(original)
    @test details(changed).formulations.effective.internal_impedance.transfer===:TestTransfer
    candidates=Formulation(internal_impedance=Grid((formula(:default),customized));combine=:zip)
    @test length(candidates)==2
    @test collect(candidates)[2].methods.internal_impedance.transfer.hooks.transfer===transfer

    native=Formulation(internal_impedance=customized)
    empty!(calls)
    for source in (native,MonteCarlo(native),LinearError(native))
        saved=IO.deserialize_value(Val(:formulation),NamedTuple(source))
        @test description([source];quantity=R)==description([saved];quantity=R)
        @test occursin("internal Z(transfer)",only(description([source];quantity=R)))
        @test !occursin("internal Z",only(description([source];quantity=B)))
    end
    @test isempty(calls)
    # Historical internal names are interpreted only in the saved boundary;
    # the neighboring earth mutual vocabulary must never be rewritten.
    declaration=NamedTuple(native)
    old=merge(declaration,(requested=merge(declaration.requested,
        (internal_impedance=(identifier=:default,hooks=(mutual=transfer,)),)),
        methods=merge(declaration.methods,(internal_impedance=(identifier=:default,),))))
    saved=IO.deserialize_value(Val(:formulation),old)
    internal=only(value for (scope,value) in pairs(saved...) if last(scope)==(:internal_impedance,))
    @test formulation_options(internal).hooks==(transfer=transfer,)
    @test old.requested.internal_impedance.hooks==(mutual=transfer,)
    @test isempty(calls)
end
