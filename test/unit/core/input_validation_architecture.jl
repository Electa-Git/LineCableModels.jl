@testitem "Core / InputValidation / owned inputs validate without mutation" tags=[:unit] setup=[TestFixtures] begin
    using RequiredInterfaces, Serialization
    const IV=LineCableModels.InputValidation
    const E=LineCableModels.Engine
    @test parentmodule(validate) === IV
    @test !applicable(validate, (answer=42,))
    for root in (AbstractMaterial, AbstractEarthModel, AbstractProblemDefinition)
        @test RequiredInterfaces.isInterface(root)
        @test Tuple(RequiredInterfaces.functions(RequiredInterfaces.getInterface(root))) == (validate,)
    end
    design=TestFixtures.coaxial_design()
    system=TestFixtures.three_phase_system()
    problem=TestFixtures.line_parameters_problem(system)
    values=(TestFixtures.conductor_material(), homogeneous(rho=100.0),
        layer(rho=200.0), design, system, problem, CableConstantsProblem(design),
        TestFixtures.two_conductor_results())
    for input in values
        owner=parentmodule(typeof(input))
        @test which(validate, Tuple{typeof(input)}).module === owner
        bytes(x)=let buffer=IOBuffer(); Serialization.serialize(buffer,x); take!(buffer); end
        before=bytes(input)
        @test applicable(validate,input)
        @test validate(input) === input
        @test bytes(input)==before
    end
    # Construction and compute must reject damaged current inputs. The public
    # materialization path for new problem implementors is tested in behavior.jl.
    problem.frequencies[1]=-1
    @test_throws DomainError validate(problem)
    @test_throws DomainError compute(problem,Formulation())
    bad=TestFixtures.coaxial_design()
    pop!(bad.terminal_map)
    snapshot=copy(bad.terminal_map)
    @test_throws DimensionMismatch validate(bad)
    @test bad.terminal_map == snapshot
end
