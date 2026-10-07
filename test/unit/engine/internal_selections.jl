@testitem "Engine / one internal formula serves its surfaces through native dispatch" tags=[:unit, :importexport, :slow] setup=[TestFixtures,FormulaFixtures] begin
    const II=LineCableModels.Engine.InternalImpedance
    const M=FormulaFixtures
    args=(0.008,0.01,1.7241e-8,1.,100.0im)
    scalar=II.Formula(:default)
    reference=@inferred II.surface_impedances(scalar,Val((:outer,:transfer,:inner)),args...)
    @test II.surface_impedances(scalar,args...)==reference[(:inner,:outer,:transfer)]
    # A slot takes one formula. Its surface impedances are that formula's expressions.
    @test_throws MethodError Formulation(internal_impedance=(inner=:default,outer=:default,transfer=:default))
    coefficients=(inner=11+12im,outer=21+22im,transfer=31+32im)
    custom=M.SurfaceLaw(;coefficients)
    @test II.surface_impedances(custom,args...)==coefficients
    @test length(custom.state_inputs)==1 && length(custom.evaluations)==3
    outer_only=M.SurfaceLaw(kinds=(:outer,);coefficients)
    @test (@inferred II.surface_impedances(outer_only,Val((:outer,)),args...)).outer==coefficients.outer
    @test_throws MethodError II.surface_impedances(outer_only,args...)

    struct ObservedSurfaces{F,P,O} <: LineCableModels.Engine.InternalImpedanceFormulation
        base::F
        parameters::P
        options::O
        observations::Vector{Tuple}
    end
    observed=ObservedSurfaces(scalar,(;),scalar.options,Tuple[])
    function (leaf::ObservedSurfaces)(args...)
        surface_functor=leaf.base(args...)
        II.Functor(leaf,surface_functor.state,leaf.options)
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

    problem=TestFixtures.line_parameters_problem(frequencies=[50.,500.])
    original=compute(problem,Formulation())
    explicit=compute(problem,Formulation(internal_impedance=:default))
    @test Z(explicit)==Z(original) && Y(explicit)==Y(original)
    cable_problem=CableConstantsProblem(first(problem.system.designs);frequency=50.)
    local_default=compute(cable_problem,CableConstantsFormulation())
    local_explicit=compute(cable_problem,CableConstantsFormulation(internal_impedance=:default))
    @test local_default.R==local_explicit.R && local_default.L==local_explicit.L
    changed=compute(problem,Formulation(internal_impedance=M.SurfaceLaw(;coefficients)))
    @test Z(changed)!=Z(original) && Y(changed)==Y(original)
    @test details(changed).data.formulations.methods.internal_impedance.identifier === :SurfaceLaw
    formulations=Formulation(internal_impedance=Grid((formula(:default),custom));combine=:zip)
    @test length(formulations)==2
    @test collect(formulations)[2].methods.internal_impedance === custom
    native=Formulation(internal_impedance=:default)
    for source in (native,MonteCarlo(native),LinearError(native))
        @test occursin("internal Z",only(description([source];quantity=R)))
        @test !occursin("internal Z",only(description([source];quantity=B)))
    end
end
