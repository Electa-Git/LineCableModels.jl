@testitem "Gauntlet / catalogue declarations and explicit frequencies" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport.Gauntlet
    using LineCableModels
    root=joinpath(pkgdir(LineCableModels),"gauntlet","benchmarks")
    for (directory,_,names) in walkdir(root), name in names
        endswith(name,".jl") || continue
        id=Symbol(first(splitext(name)))
        println("Catalogue declaration: ",id)
        flush(stdout)
        benchmark=benchmark_definition(id;frequencies=[.1,1.,7.,100.])
        @test benchmark.id === id
        @test benchmark.reference.problem === benchmark.candidate.problem
        problem=benchmark.reference.problem
        if problem isa LineParametersProblem
            @test problem.frequencies == [.1,1.,7.,100.]
            if benchmark.reference.formulation isa PSCAD.PSCADFormulation
                @test_throws ArgumentError validate(problem,benchmark.reference.formulation)
            end
        else
            @test first(problem.space).frequencies == [.1,1.,7.,100.]
        end
    end
    # Catalogue transport is exhaustive above. One current analytical calculation
    # checks the execution path; solving every cable adds no frequency contract.
    basic = benchmark_definition(:benchmark_two_bare_wires_fem;
        frequencies=[.1,1.,7.,100.]).candidate
    @test frequencies(compute(basic.problem,first(basic.formulation);options=basic.options)) ==
        [.1,1.,7.,100.]
    fem = benchmark_definition(:benchmark_two_bare_wires_fem;
        frequencies=[50.0], reference_options=(frequency_workers=1, mesh_policy=:reuse,))
    @test !hasproperty(fem.reference.formulation, :execution)
    @test fem.reference.options isa ComputationOptions
    execution = computation_options(LineCableModelsFEM, fem.reference.options)
    @test execution.data.frequency_workers == 1
    @test execution.data.mesh_policy === :reuse
    @test execution.data.trace === Val(true)
    @test execution.data.keep_run_directory

end
