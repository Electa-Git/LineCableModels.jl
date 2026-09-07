@testitem "Gauntlet / UQ campaign preflight and scientific settings" tags=[:gauntlet_toolkit] begin
    using LineCableModels
    include(joinpath(pkgdir(LineCableModels), "test", "gauntlet", "runner.jl"))
    owner = GauntletSupport
    inner = Formulation()
    settings = Dict("seed"=>"1234", "trials"=>4)
    linear = owner.campaign_propagation(:linear_error, :coaxial, inner, false)
    sampled = owner.campaign_propagation(:monte_carlo, :coaxial, inner, settings)
    @test linear isa LinearError
    @test sampled isa MonteCarlo
    @test sampled.seed == 1234
    @test sampled.trials == 4
    @test sampled.options.on_error === :fail
    @test owner.campaign_propagation(:deterministic, :fem, inner, false) === inner
    @test_throws ArgumentError owner.campaign_propagation(:monte_carlo, :fem, inner, settings)
    @test_throws ArgumentError owner.campaign_propagation(:linear_error, :fem, inner, false)
    @test_throws ArgumentError owner.campaign_propagation(:linear_error, :pscad, inner, false)
    @test_throws ArgumentError owner.campaign_propagation(:unknown, :coaxial, inner, false)
    @test owner.campaign_propagation(:monte_carlo, :pscad, inner, settings) isa MonteCarlo
    @test owner.campaign_implementation(linear).selection_sha256 !=
        owner.campaign_implementation(sampled).selection_sha256
    changed = owner.campaign_propagation(:monte_carlo, :coaxial, inner, merge(settings, Dict("seed"=>"1235")))
    @test owner.campaign_implementation(changed).selection_sha256 !=
        owner.campaign_implementation(sampled).selection_sha256
    paths = getproperty.(owner.campaign_implementation(sampled).blobs, :path)
    @test "src/grid.jl" in paths
    @test "src/uq/montecarlo/compute.jl" in paths
    @test "src/uq/linearerror.jl" ∉ paths
    for formulation in (linear, sampled)
        implementation_paths = getproperty.(owner.campaign_implementation(formulation).blobs, :path)
        # These observations determine the RLCG moments recorded by UQ. Plotting
        # and rebinned publication products do not determine those moments.
        @test "src/engine/lineparameters/lineparameters.jl" in implementation_paths
        @test "src/uq/publication.jl" ∉ implementation_paths
        @test "src/uq/observations.jl" ∉ implementation_paths
        @test "ext/LineCableModelsMakieExt/montecarlo.jl" ∉ implementation_paths
    end
    vector = [1.0, 2.0, 3.0, 4.0]
    @test owner.semantic_sha256(reshape(vector, 2, 2)) != owner.semantic_sha256(reshape(vector, 1, 4))
    @test owner.semantic_sha256(reshape(vector, 2, 2)) != owner.semantic_sha256(vector)
    mktempdir() do root
        destination = joinpath(root, "forbidden")
        variation = owner.RelativeStandardUncertainty(1.0; tags=(:geometry, :cable_layer))
        @test_throws ArgumentError owner.run_campaign(destination, [:two_bare_wires];
            backends=(:fem,), catalogue=false, propagation=(:monte_carlo,),
            uncertainty=variation, seed=1234, trials=4)
        @test !ispath(destination)
        @test_throws ArgumentError owner.run_campaign(destination, [:two_bare_wires];
            backends=(:coaxial,), catalogue=false, propagation=(:linear_error,))
        @test !ispath(destination)
        @test_throws ArgumentError owner.run_campaign(destination, [:two_bare_wires];
            backends=(:coaxial,), catalogue=false, propagation=(:monte_carlo,), uncertainty=variation)
        @test !ispath(destination)
    end
end

@testitem "Gauntlet / UQ campaign results, exact replay and partial resumption" tags=[:gauntlet] begin
    using LineCableModels, JLD2, SHA
    include(joinpath(pkgdir(LineCableModels), "test", "gauntlet", "runner.jl"))
    owner = GauntletSupport
    mktempdir() do temporary
        root = joinpath(temporary, "uq")
        variation = owner.RelativeStandardUncertainty(1.0; tags=(:geometry, :cable_layer))
        choices = (internal_impedance=Grid((:default, :Schelkunoff1934)),)
        @test owner.run_campaign(root, [:two_bare_wires]; backends=(:coaxial,),
            catalogue=false, choices, propagation=(:deterministic, :linear_error, :monte_carlo),
            uncertainty=variation, trials=4, seed=1234)
        @test all(row -> row.completed == row.requested == 2, owner.campaign_status(root))
        _, plan = owner.campaign_plan(root)
        @test plan["uncertainty"]["percent"] == 1.0
        @test plan["uncertainty"]["tags"] == ["geometry", "cable_layer"]
        @test plan["monte_carlo"] == Dict("trials"=>4, "seed"=>"1234")
        model = owner.campaign_models([:two_bare_wires], variation)[(:two_bare_wires, true)]
        @test model.problem isa Gridspace
        @test has_uncertainty(model.problem)
        @test owner.campaign_input(model, :linear_error) == owner.campaign_input(model, :monte_carlo)
        for method in (:linear_error, :monte_carlo)
            job = only(filter(job -> job["propagation"] == string(method), plan["jobs"]))
            records = [JLD2.load(joinpath(root, job["id"], lpad(index, 4, '0') * ".jld2")) for index in 1:2]
            @test all(record -> record["kind"] === :gauntlet_moments, records)
            @test records[1]["moments"] == records[2]["moments"]
            @test records[1]["numerical_reference_approval"] === :unreviewed
            restored = LineCableModels.ImportExport.deserialize_value(records[1]["problem"])
            @test owner.numerical_input_sha256(restored) == owner.numerical_input_sha256(model.nominal_problem)
            @test records[1]["formulation"].definitions.internal_impedance === :default
            @test length(records[1]["frequencies"]) == 101
            @test records[1]["correlation"].rule === :parameter_identity
            inner = owner.campaign_formulation(:coaxial, first(job["selections"]), :default)
            direct = compute(ParametricProblem(model.problem),
                owner.campaign_propagation(method, :coaxial, inner, plan["monte_carlo"]))
            @test records[1]["moments"] == NamedTuple(owner.extract_moments(direct, model.port_order))
            if method === :monte_carlo
                @test records[1]["sampling"].root_seed == 1234
                @test records[1]["sampling"].trial_counts == [4]
                @test records[1]["sampling"].point_seeds == [1234]
                @test length(only(records[1]["computation_details"].trials)) == 4
            else
                @test records[1]["sampling"] === nothing
            end
        end
        files = [joinpath(root, job["id"], lpad(index, 4, '0') * ".jld2")
            for job in plan["jobs"] for index in 1:2]
        before = [(sha256(read(file)), stat(file).mtime) for file in files]
        @test owner.resume_campaign(root)
        @test [(sha256(read(file)), stat(file).mtime) for file in files] == before

        partial = joinpath(temporary, "partial")
        partial_plan = deepcopy(plan)
        partial_plan["jobs"] = filter(job -> job["propagation"] == "monte_carlo", partial_plan["jobs"])
        partial_plan["jobs"][1]["selections"][2]["earth_impedance"] = "DirectNumericalIntegration"
        owner.write_campaign_state(joinpath(partial, "campaign.toml"), partial_plan)
        @test !owner.resume_campaign(partial)
        @test only(owner.campaign_status(partial)).completed == 1
        saved = joinpath(partial, partial_plan["jobs"][1]["id"], "0001.jld2")
        saved_bytes, saved_mtime = read(saved), stat(saved).mtime
        @test !owner.resume_campaign(partial)
        @test read(saved) == saved_bytes
        @test stat(saved).mtime == saved_mtime
        @test !isfile(joinpath(dirname(saved), "0002.jld2"))

        plan["monte_carlo"]["seed"] = "1235"
        owner.write_campaign_state(joinpath(root, "campaign.toml"), plan)
        @test !owner.resume_campaign(root)
        @test [(sha256(read(file)), stat(file).mtime) for file in files] == before
        @test only(filter(row -> endswith(row.id, "monte_carlo"), owner.campaign_status(root))).state == "failed"
    end
end
