@testitem "PSCAD / matrix coordinates and complete native reuse keys" tags=[:integration] begin
    using TOML, SHA, Measurements
    const P = LineCableModels.PSCAD
    const launches = Ref(0)
    const identity = Dict("schema"=>"1", "version"=>"5.1.0", "pscad_sha256"=>repeat("1",64),
        "line_constants_sha256"=>repeat("2",64), "master_library_sha256"=>repeat("3",64))
    # Deliberately nonsymmetric protocol matrices, not scientific references.
    function P.remote_command(::Val{:ordering_fixture}, ::P.RemoteConfig, command::AbstractString)
        if occursin("-SharedCase", command)
            launches[] += 1
            root = replace(only(match(r"-SharedCase '([^']+)'", command).captures), '\\'=>'/')
            script = """
                using TOML
                root = $(repr(root))
                input = TOML.parsefile(joinpath(root, "computation.toml"))
                output = joinpath(root, "outputs")
                n = input["matrix_size"][1]
                identity = TOML.parse($(repr(sprint(TOML.print, identity))))
                get(input, "expected_solver", identity) == identity || error("unexpected identity")
                for (name, scale) in (("zm",1.0), ("zp",0.0), ("ym",1e-6), ("yp",0.0))
                    open(joinpath(output,"result_"*name*".out"),"w") do io
                        println(io,"LOG10(FN) FN matrix")
                        for (k,f) in enumerate(input["frequencies"])
                            println(io,join((log10(f), f, (scale*(100r+10c+k) for r in 1:n for c in 1:n)...),' '))
                        end
                    end
                end
                readback = Dict(name=>Dict(field=>control["readback"] for (field,control) in fields)
                    for (name,fields) in input["native_settings"])
                open(io->TOML.print(io,readback),joinpath(output,"native-settings.toml"),"w")
                open(io->TOML.print(io,identity),joinpath(output,"solver.toml"),"w")
                write(joinpath(output,"timing.txt"),"1.25")
                write(joinpath(output,"pscad-console.txt"),"protocol fixture")
                """
        else
            script = "print(" * repr(sprint(TOML.print, identity)) * ")"
        end
        `$(Base.julia_cmd()) --startup-file=no --project=@stdlib -e $script`
    end
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
        remote = P.RemoteConfig("fixture",directory,"scratch","julia","python";
            local_root=directory,transport=:ordering_fixture)
        selected = Formulation(:pscad)
        first_result = compute(permuted,selected;options=(;remote))
        @test Z(first_result)[1,2,1] == 231
        @test Z(first_result)[2,1,1] == 321
        @test Y(first_result)[1,3,1] ≈ 211e-6
        @test details(first_result).data.coordinates == ["cable:1:sheath","cable:2:core","cable:1:core"]
        source = details(first_result).data.execution.source_run
        raw_hash = bytes2hex(open(sha256,joinpath(source,"outputs","result_zm.out")))
        reused = compute(natural,selected;options=(;remote,resume_run_directory=source))
        @test launches[] == 1
        @test Z(reused)[1,2,1] == 121
        @test Y(reused)[3,1,1] ≈ 311e-6
        @test bytes2hex(open(sha256,joinpath(source,"outputs","result_zm.out"))) == raw_hash
        total = compute(permuted,selected;options=(;remote,resume_run_directory=source,output_basis=:total))
        @test Z(total) == 42 .* Z(first_result)
        @test Y(total) == 42 .* Y(first_result)
        repeated = compute(permuted,selected;options=(;remote))
        @test launches[] == 2
        @test !details(repeated).data.execution.reused
        @test details(repeated).data.execution.source_run != source
        @test Z(repeated) == Z(first_result)
        callbacks = Int[]
        batch = compute(permuted,[selected,Formulation(:pscad;earth_impedance=:pollaczek1926)];
            options=(;remote,on_result=(p,i,r)->push!(callbacks,i)))
        @test launches[] == 3
        @test callbacks == [1,2]
        @test Z(batch[1]) == Z(batch[2]) && Z(batch[1]) !== Z(batch[2])
        @test Y(batch[1]) == Y(batch[2]) && Y(batch[1]) !== Y(batch[2])
        @test details(batch[1]).data.gridpoint.source_id == details(batch[2]).data.gridpoint.source_id
        @test [details(result).data.gridpoint.formulation_index for result in batch] == [1,2]
        @test details(batch[2]).data.formulations.requested.earth_impedance.identifier === :pollaczek1926
        near = [Formulation(:pscad;options=(base_frequency=f,)) for f in (50.000001,50.000002)]
        preparations = [P._prepare_pscad(permuted,f,P._pscad_blueprints(permuted.system)) for f in near]
        @test preparations[1].project == preparations[2].project
        distinct = compute(permuted,near;options=(;remote))
        @test launches[] == 5
        @test !any(r->details(r).data.execution.reused,distinct)
        @test [details(r).data.base_frequency for r in distinct] == [50.000001,50.000002]
        @test_throws r"callback failed" compute(permuted,selected;
            options=(;remote,on_result=(args...)->error("callback failed")))
        @test launches[] == 6
        grid = Formulation(:pscad;earth_impedance=Grid((:default,:pollaczek1926)))
        traversed = compute(permuted,grid;options=(;remote))
        @test length(traversed) == 2
        @test Z(first(traversed)) == Z(first_result)
        @test [details(result).data.gridpoint.formulation_index for result in traversed] == [1,2]
        function realized_problem(radius)
            T = typeof(radius)
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
        parametric = ParametricProblem(space,ComputationOptions((;remote,
            on_result=(p,i,r)->push!(realized_types,eltype(p)))))
        before = launches[]
        deterministic = Gridspace{LineParametersProblem}(realized_problem,(Grid((0.004,)),))
        nominal = compute(ParametricProblem(deterministic,parametric.options),Combinatorial(selected))
        @test only(nominal) isa LineParameters
        sampled = compute(parametric,MonteCarlo(selected;trials=2,seed=71,return_samples=true))
        @test sampled isa MonteCarloResult{<:LineParameters}
        @test launches[] == before + 3
        @test length(realized_types) == 3
        @test all(T->!LineCableModels.Engine.has_uncertainty_type(T),realized_types)
        @test_throws r"PSCAD does not support Measurement" compute(parametric,Combinatorial(selected))
        @test_throws r"PSCAD does not support Measurement" compute(parametric,LinearError(selected))
        @test launches[] == before + 3
        oversized = build(LineCableSystem,[coax,wire],[Pose2(0,-1),Pose2(1,-1)];
            connections=[Dict(:core=>1,:sheath=>2),Dict(:core=>3)],line_length=floatmax(Float64))
        overflowing = LineParametersProblem(oversized;earth_props=natural.earth_props,
            frequencies=natural.frequencies)
        overflow_callbacks = Int[]
        @test_throws DomainError compute(overflowing,selected;options=(;remote,output_basis=:total,
            on_result=(p,i,r)->push!(overflow_callbacks,i)))
        @test isempty(overflow_callbacks)
        # A completion from the former protocol never qualifies for new reuse.
        old = joinpath(directory,"old-protocol"); mkpath(old)
        open(io->TOML.print(io,Dict("schema_version"=>3)),joinpath(old,"complete.toml"),"w")
        @test_throws r"fresh version-4" compute(permuted,selected;options=(;remote,resume_run_directory=old))
        before = launches[]
        latest = compute(permuted,selected;options=(;remote,resume_run_directory=:latest))
        @test details(latest).data.execution.reused
        @test details(latest).data.execution.source_run != old
        @test launches[] == before
        record_path = joinpath(source,"complete.toml")
        corrupted = TOML.parsefile(record_path)
        corrupted["request_sha256"] = repeat("0",64)
        open(io->TOML.print(io,corrupted),record_path,"w")
        @test_throws r"input integrity check failed" compute(permuted,selected;
            options=(;remote,resume_run_directory=source))
    end
end
