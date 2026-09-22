@testitem "PSCAD / real station repair acceptance" tags=[:pscad_native] begin
    using TOML, SHA, CairoMakie, LinearAlgebra
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    const P = LineCableModels.PSCAD
    config_path = get(ENV, "LINECABLEMODELS_PSCAD_CONFIG", "")
    isempty(config_path) && error("Set LINECABLEMODELS_PSCAD_CONFIG to your station TOML filename")
    remote = P.RemoteConfig(config_path)
    copper = Material(:conductor, 1.72e-8, 1.0, 1.0, 20.0, 0.004)
    dielectric = Material(:insulator, 1e8, 2.3; tan_delta=0.03)
    designs = [build(CableDesign,"native-wire-$i",terminal(:core,
        core(copper;r=0.003 + i*0.0005),insulation(dielectric;t=0.002))) for i in 1:4]
    function native_problem(order)
        system = build(LineCableSystem,designs,[Pose2(0,2),Pose2(1,3),Pose2(2,-1),Pose2(3,-2)];
            connections=[Dict(:core=>i) for i in order],line_length=42.0,system_id="native-repair")
        LineParametersProblem(system;earth_props=homogeneous(rho=100.0,eps_r=10.0),
            temperature=60.0,frequencies=10.0 .^ range(log10(50.0),log10(1000.0);length=101))
    end
    problem = native_problem([1,2,3,4])
    options = (remote=remote,verbosity=(default=0,PSCAD=1))
    selected = Formulation(:pscad)
    first_timing = @timed compute(problem,selected;options)
    first_result = first_timing.value
    next_timing = @timed compute(problem,selected;options)
    next_result = next_timing.value
    @test !details(first_result).data.execution.reused
    @test !details(next_result).data.execution.reused
    @test details(first_result).data.execution.source_run != details(next_result).data.execution.source_run
    @test Z(first_result) == Z(next_result)
    @test Y(first_result) == Y(next_result)
    @test details(first_result).data.formulations.methods.earth_admittance.identifier === :ideal
    @test details(first_result).data.formulations.methods.internal_impedance.identifier === :wedepohl1973
    @test map(value -> value.identifier, details(first_result).data.formulations.methods.earth_impedance) ==
        (air=:carson1926, earth=:pollaczek1926, mixed=:lucca1994)
    # Independent electrostatic reference: air images plus each wire's insulation.
    # Buried and mixed external coefficients are zero. Matrix output has 8-digit
    # magnitudes, so the relative comparison allows their printed rounding.
    epsilon0 = 8.8541878128e-12
    radii = [0.003 + i*0.0005 for i in 1:4]
    exterior = radii .+ 0.002
    heights = [2.0,3.0,-1.0,-2.0]
    potential = Matrix(Diagonal(log.(exterior ./ radii) ./ (2pi * epsilon0 * 2.3)))
    for i in 1:2, j in 1:2
        potential[i,j] += i == j ? log(2heights[i]/exterior[i])/(2pi*epsilon0) :
            log(hypot(i-j, heights[i]+heights[j])/hypot(i-j, heights[i]-heights[j]))/(2pi*epsilon0)
    end
    for k in (1,51,101)
        expected_y = (2pi*problem.frequencies[k]*im) .* inv(potential)
        # Direct integration has a measured aerial conductance. Compare its
        # capacitance here and retain the full complex output below.
        @test imag.(Y(first_result)[:,:,k]) ≈ imag.(expected_y) rtol=1e-7
        @test iszero(Y(first_result)[3,4,k])
        @test iszero(Y(first_result)[1,3,k])
    end
    source = details(first_result).data.execution.source_run
    source_hash = bytes2hex(open(sha256,joinpath(source,"outputs","result_zm.out")))
    reuse_timing = @timed compute(native_problem([3,1,4,2]),selected;
        options=(;options...,resume_run_directory=source))
    reordered = reuse_timing.value
    @test details(reordered).data.execution.reused
    @test Z(reordered) == Z(first_result)[[2,4,1,3],[2,4,1,3],:]
    @test Y(reordered) == Y(first_result)[[2,4,1,3],[2,4,1,3],:]
    @test bytes2hex(open(sha256,joinpath(source,"outputs","result_zm.out"))) == source_hash
    callbacks = Int[]
    selections = [selected,
        Formulation(:pscad;earth_impedance=(air=:carson1926,earth=:pollaczek1926,mixed=:lucca1994)),
        Formulation(:pscad;earth_impedance=(air=:gary1976,earth=:wedepohl1973,mixed=:ametani2009)),
        Formulation(:pscad;earth_impedance=(air=:gary1976,earth=:saad1996,mixed=:lucca1994),
            insulation_admittance=:lossy)]
    batch = compute(problem,selections;options=(;options...,resume_run_directory=:latest,
        on_result=(p,i,r)->push!(callbacks,i)))
    # The first two selections share native settings. Distinct settings require
    # separate native outputs, which :latest may reuse from an earlier test run.
    @test callbacks == collect(1:length(selections))
    @test Z(batch[1]) == Z(batch[2])
    @test Z(batch[1]) !== Z(batch[2])
    # The Gary/Wedepohl case agrees with strict ideal images. Direct integration
    # departs from that law; the adapter must retain the measured conductance.
    for k in (1,51,101)
        expected_y = (2pi*problem.frequencies[k]*im) .* inv(potential)
        @test Y(batch[3])[:,:,k] ≈ expected_y rtol=1e-7
        @test iszero(Y(batch[3])[3,4,k])
        @test iszero(Y(batch[3])[1,3,k])
    end
    @test real(Y(first_result)[1,1,end]) > 0
    @test maximum(abs, real.(Y(batch[3]))) <= eps(Float64)*maximum(abs,Y(batch[3]))
    for result in (first_result,next_result,reordered,batch...)
        @test result isa LineParameters
        @test size(Z(result)) == size(Y(result)) == (4,4,101)
        @test all(isfinite,Z(result)) && all(isfinite,Y(result))
        @test frequencies(result) == problem.frequencies
        @test details(result).data.execution.solver_identity["version"] == "5.1.0"
    end
    @test any(loss.requested > 0 for values in details(last(batch)).data.dielectric_losses for loss in values)
    total = compute(problem,selected;options=(;options...,resume_run_directory=source,output_basis=:total))
    @test Z(total) == 42 .* Z(first_result)
    @test Y(total) == 42 .* Y(first_result)
    observed = ObservedResult(first_result)
    artifact = report(BenchmarkTableDefinition(),(reference=first_result,candidate=last(batch)))
    @test !isempty(artifact.tables.formulations.label)
    pages = LineCableModels.plot(artifact,(Z,);backend=:cairo,display_plot=false,
        controls=false,open_export=false)
    @test !isempty(pages)
    for (index,page) in enumerate(pages)
        CairoMakie.save(joinpath(remote.local_root,"native-$(basename(source))-$index.png"),page.figure)
    end
    evidence = Dict("transport"=>string(remote.transport), "source_run"=>source,
        "first_compute_seconds"=>first_timing.time,"next_compute_seconds"=>next_timing.time,
        "reuse_seconds"=>reuse_timing.time,
        "first_compile_seconds"=>details(first_result).data.execution.elapsed_seconds,
        "next_compile_seconds"=>details(next_result).data.execution.elapsed_seconds,
        "runs"=>[details(result).data.execution.source_run for result in batch],
        "ideal_deviation"=>Dict(
            "frequency_hz"=>problem.frequencies[end],
            "direct_y11_real_s_per_m"=>real(Y(first_result)[1,1,end]),
            "direct_y11_imag_s_per_m"=>imag(Y(first_result)[1,1,end]),
            "direct_y11_phase_degrees"=>rad2deg(angle(Y(first_result)[1,1,end])),
            "gary_y11_phase_degrees"=>rad2deg(angle(Y(batch[3])[1,1,end])),
            "direct_max_aerial_conductance_s_per_m"=>maximum(abs,real.(Y(first_result)[1:2,1:2,:]))))
    evidence_path = joinpath(remote.local_root,"acceptance-$(basename(source)).toml")
    open(io->TOML.print(io,evidence;sorted=true),evidence_path,"w")
    println("Native PSCAD acceptance evidence: ",evidence_path)
    println("First compute / next compute / reuse seconds: ",
        (first_timing.time,next_timing.time,reuse_timing.time))
end

@testitem "PSCAD / real station expected identity precedes compilation" tags=[:pscad_native] begin
    const P = LineCableModels.PSCAD
    remote = P.RemoteConfig(ENV["LINECABLEMODELS_PSCAD_CONFIG"])
    expected = P.identify(remote)
    expected["line_constants_sha256"] = repeat("0",64)
    metal = Material(:conductor,1.72e-8,1.0)
    dielectric = Material(:insulator,1e14,2.3)
    design = build(CableDesign,"identity-guard",terminal(:core,
        core(metal;r=0.004),insulation(dielectric;t=0.002)))
    system = build(LineCableSystem,[design],[Pose2(0,-1)];connections=[Dict(:core=>1)])
    problem = LineParametersProblem(system;earth_props=homogeneous(rho=100.0),
        frequencies=10.0 .^ range(0,3;length=101))
    work_root = mktempdir(mkpath(remote.local_root);prefix="identity-rejection-",cleanup=false)
    @test_throws r"does not match the expected solver identity" compute(problem,Formulation(:pscad);
        options=(;remote,work_root,solver_identity=expected))
    run_directory = only(readdir(work_root;join=true))
    @test !isfile(joinpath(run_directory,"complete.toml"))
    console = read(joinpath(run_directory,"outputs","pscad-console.txt"),String)
    @test count("Launching PSCAD 5.1.0",console) == 1
    @test !occursin("Starting PSCAD line-constants calculation",console)
    @test !isfile(joinpath(run_directory,"outputs","result_zm.out"))
    println("Native identity rejection evidence: ",run_directory)
end
