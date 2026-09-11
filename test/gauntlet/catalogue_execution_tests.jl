@testitem "Gauntlet / catalogue declarations and explicit frequencies" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport.Gauntlet
    using LineCableModels
    root=joinpath(pkgdir(LineCableModels),"gauntlet","benchmarks")
    for (directory,_,names) in walkdir(root), name in names
        endswith(name,".jl") || continue
        id=Symbol(first(splitext(name)))
        benchmark=benchmark_definition(id;frequencies=[.1,1.,7.,100.])
        @test benchmark.id === id
        @test benchmark.reference.problem === benchmark.candidate.problem
        problem=benchmark.reference.problem
        if problem isa LineParametersProblem
            @test problem.frequencies == [.1,1.,7.,100.]
            if benchmark.reference.formulation isa PSCAD.PSCADFormulation
                @test_throws ArgumentError validate(problem,benchmark.reference.formulation)
            elseif benchmark.reference.formulation isa LineCableModelsFEM
                @test problem.frequencies == [.1,1.,7.,100.]
            else
                @test frequencies(compute(problem,benchmark.reference.formulation)) == problem.frequencies
            end
        else
            @test first(problem.space).frequencies == [.1,1.,7.,100.]
        end
    end
    fem = benchmark_definition(:benchmark_two_bare_wires_fem;
        frequencies=[50.0], reference_options=(frequency_workers=1, mesh_policy=:reuse,))
    @test !hasproperty(fem.reference.formulation, :execution)
    @test fem.reference.options isa ComputationOptions
    execution = computation_options(LineCableModelsFEM, fem.reference.options)
    @test execution.frequency_workers == 1
    @test execution.mesh_policy === :reuse
    @test execution.trace === Val(true)
    @test execution.keep_run_directory

end
