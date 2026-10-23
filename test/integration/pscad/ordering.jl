@testitem "PSCAD / matrix coordinates and complete native reuse keys" tags=[:integration, :measurements, :slow] setup=[PSCADStandIn] begin
    using LineCableModels.Engine: completed_formulation
    using TOML, SHA, Measurements
    const P = LineCableModels.PSCAD
    const S = PSCADStandIn
    # The stand-in's matrices are deliberately nonsymmetric protocol matrices, not
    # scientific references.
    copper = Material(:conductor,1.72e-8,1.0)
    insulator = Material(:insulator,1e14,2.3)
    coax = build(CableDesign,"coax", terminal(:core,core(copper;r=0.004),insulation(insulator;t=0.002)),
        terminal(:sheath,sheath(copper;t=0.001),insulation(insulator;t=0.001)))
    wire = build(CableDesign,"wire",terminal(:core,core(copper;r=0.003),insulation(insulator;t=0.002)))
    function problem(connections)
        system = build(LineCableSystem,[coax,wire],[Pose2(0,-1),Pose2(1,-1)];
            connections, line_length=42.0, system_id="ordering")
        LineParametersProblem(system; earth_props=homogeneous(rho=100.0),
            frequencies=10.0 .^ range(-1,5;length=101))
    end
    permuted = problem([Dict(:core=>3,:sheath=>1),Dict(:core=>2)])
    natural = problem([Dict(:core=>1,:sheath=>2),Dict(:core=>3)])
    mktempdir() do directory
        application = S.install(directory)
        launches() = S.compiles(application)
        remote = S.remote(directory)
        selected = Formulation(:pscad)
        first_result = compute(permuted,selected;options=(;remote))
        @test Z(first_result)[1,2,1] == 231
        @test Z(first_result)[2,1,1] == 321
        @test Y(first_result)[1,3,1] ≈ 211e-6
        @test first_result.terminals == ["cable:1:sheath","cable:2:core","cable:1:core"]
        source = details(first_result).data.execution.source_run
        raw_hash = bytes2hex(open(sha256,joinpath(source,"outputs","result_zm.out")))
        reused = compute(natural,selected;options=(;remote,resume_run_directory=source))
        @test launches() == 1
        @test Z(reused)[1,2,1] == 121
        @test Y(reused)[3,1,1] ≈ 311e-6
        @test bytes2hex(open(sha256,joinpath(source,"outputs","result_zm.out"))) == raw_hash
        total = compute(permuted,selected;options=(;remote,resume_run_directory=source))
        @test Z(total) .* line_length(total) == 42 .* Z(first_result)
        @test Y(total) .* line_length(total) == 42 .* Y(first_result)
        repeated = compute(permuted,selected;options=(;remote))
        @test launches() == 2
        @test !details(repeated).data.execution.reused
        @test details(repeated).data.execution.source_run != source
        @test Z(repeated) == Z(first_result)
        callbacks = Int[]
        batch = compute(permuted,Grid((selected,Formulation(:pscad;earth_impedance=:pollaczek1926)));
            options=(;remote,on_result=(p,i,r)->push!(callbacks,i)))
        @test launches() == 3
        @test callbacks == [1,2]
        @test Z(batch[1]) == Z(batch[2]) && Z(batch[1]) !== Z(batch[2])
        @test Y(batch[1]) == Y(batch[2]) && Y(batch[1]) !== Y(batch[2])
        @test details(batch[1]).data.gridpoint.source_id == details(batch[2]).data.gridpoint.source_id
        @test [details(result).data.gridpoint.formulation_index for result in batch] == [1,2]
        @test completed_formulation(details(batch[2])).formulations.formulas.earth_impedance.identifier === :pollaczek1926
        near = [Formulation(:pscad;options=(base_frequency=f,)) for f in (50.000001,50.000002)]
        project_inputs = [P._pscad_inputs(permuted,f,P._pscad_blueprints(permuted.system)) for f in near]
        @test project_inputs[1].project == project_inputs[2].project
        distinct = compute(permuted,Grid(Tuple(near));options=(;remote))
        @test launches() == 5
        @test !any(r->details(r).data.execution.reused,distinct)
        @test [details(r).data.base_frequency for r in distinct] == [50.000001,50.000002]
        @test_throws r"callback failed" compute(permuted,selected;
            options=(;remote,on_result=(args...)->error("callback failed")))
        @test launches() == 6
        grid = Formulation(:pscad;earth_impedance=Grid((:default,:pollaczek1926)))
        traversed = compute(permuted,grid;options=(;remote))
        @test length(traversed) == 2
        @test Z(first(traversed)) == Z(first_result)
        @test [details(result).data.gridpoint.formulation_index for result in traversed] == [1,2]
        # The builder records the scalar type of each realization that it builds, the
        # deterministic point and the Monte Carlo nominal and trials included.
        built_types = Type[]
        function realized_problem(radius)
            T = typeof(radius)
            push!(built_types, T)
            design = build(CableDesign,"realization",terminal(:core,
                core(Material(:conductor,convert(T,1.72e-8),1.0);r=radius),
                insulation(insulator;t=0.002)))
            system = build(LineCableSystem,[design],[Pose2(zero(T),-one(T))];
                connections=[Dict(:core=>1)])
            LineParametersProblem(system;earth_props=homogeneous(rho=100.0),
                frequencies=natural.frequencies)
        end
        space = Gridspace{LineParametersProblem}(realized_problem,
            (Grid((0.004,),AbsoluteError(0.0001)),))
        realized_types = Type[]
        # The backend's options are those of the problem, and the callable is an option of
        # the run, in the parametric call.
        parametric = ParametricProblem(space,ComputationOptions((;remote)))
        callback = (on_result=(p,i,r)->push!(realized_types,eltype(p)),)
        before = launches()
        deterministic = Gridspace{LineParametersProblem}(realized_problem,(Grid((0.004,)),))
        nominal = compute(ParametricProblem(deterministic,parametric.options),Combinatorial(selected);
            options=callback)
        @test only(nominal) isa LineParameters
        sampled = compute(parametric,MonteCarlo(selected;trials=2,seed=71,return_samples=true);
            options=callback)
        @test sampled isa MonteCarloResult{<:LineParameters}
        @test launches() == before + 3
        # The deterministic point, the nominal and the two trials build with plain numbers.
        @test length(built_types) == 4
        @test all(T->!LineCableModels.Engine.has_uncertainty_type(T),built_types)
        # `on_result` runs once per point: the deterministic point and the Monte Carlo
        # point, with its nominal problem.
        @test length(realized_types) == 2
        @test all(T->!LineCableModels.Engine.has_uncertainty_type(T),realized_types)
        @test_throws r"PSCAD does not support Measurement" compute(parametric,Combinatorial(selected))
        @test_throws r"PSCAD does not support Measurement" compute(parametric,LinearError(selected))
        @test launches() == before + 3
        oversized = build(LineCableSystem,[coax,wire],[Pose2(0,-1),Pose2(1,-1)];
            connections=[Dict(:core=>1,:sheath=>2),Dict(:core=>3)],line_length=floatmax(Float64))
        overflowing = LineParametersProblem(oversized;earth_props=natural.earth_props,
            frequencies=natural.frequencies)
        overflowed = compute(overflowing,selected;options=(;remote))
        # The line's totals overflow, and the result per meter does not.
        @test !all(isfinite, Z(overflowed) .* line_length(overflowed))
        @test all(isfinite, Z(overflowed))
        # A completion from the former protocol never qualifies for new reuse.
        old = joinpath(directory,"old-protocol"); mkpath(old)
        open(io->TOML.print(io,Dict("schema_version"=>4)),joinpath(old,"complete.toml"),"w")
        @test_throws r"fresh version-5" compute(permuted,selected;options=(;remote,resume_run_directory=old))
        before = launches()
        latest = compute(permuted,selected;options=(;remote,resume_run_directory=:latest))
        @test details(latest).data.execution.reused
        @test details(latest).data.execution.source_run != old
        @test launches() == before
        record_path = joinpath(source,"complete.toml")
        corrupted = TOML.parsefile(record_path)
        corrupted["request_sha256"] = repeat("0",64)
        open(io->TOML.print(io,corrupted),record_path,"w")
        @test_throws r"input integrity check failed" compute(permuted,selected;
            options=(;remote,resume_run_directory=source))
    end
end
