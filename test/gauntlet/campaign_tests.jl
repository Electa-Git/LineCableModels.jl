@testitem "Gauntlet / manual campaign checkpoints and resume" tags=[:gauntlet] setup=[GauntletSupport] begin
    using SHA
    using JLD2
    using LineCableModels
    using .GauntletSupport
    owner=GauntletSupport
    mktempdir() do temporary
        root=joinpath(temporary, "campaign")
        @test owner.run_campaign(root, [:two_bare_wires]; backends = (:coaxial,), catalogue = false)
        rows=owner.campaign_status(root)
        @test length(rows) == 1
        @test rows[1].state == "complete"
        @test rows[1].completed == rows[1].requested == 1
        artifact=joinpath(root, "two_bare_wires_coaxial", "0001.jld2")
        bytes=read(artifact)
        modified=stat(artifact).mtime
        document=JLD2.load(artifact)
        @test document["status"] === :complete
        @test document["numerical_reference_approval"] === :unreviewed
        @test document["basis"] === :pul
        @test document["domain"] === :PhaseDomain
        @test length(document["frequencies"]) == 101
        @test first(document["frequencies"]) >= 0.1
        @test last(document["frequencies"]) ≈ 1.0e6
        @test document["selection"]["earth_impedance"] == "default"
        restored=LineCableModels.ImportExport.deserialize_value(document["problem"])
        @test restored isa LineParametersProblem
        @test owner.numerical_input_sha256(restored) == document["input_sha256"]
        declared=document["formulation"]
        replay=compute(restored, Formulation(; declared.definitions..., options = declared.options))
        @test replay.Z.values == document["Z"]
        @test replay.Y.values == document["Y"]
        @test restored.frequencies == document["frequencies"]
        @test document["elapsed_at_completion_seconds"] > 0
        @test !haskey(document, "batch_execution_seconds")
        @test isfile(artifact * ".sha256")
        @test owner.resume_campaign(root)
        @test read(artifact) == bytes
        @test length(readdir(joinpath(dirname(artifact), "attempts"))) == 1

        # A registered formula with an incompatible indexed domain can fail
        # after an earlier result from the same Gridspace batch. Repeated
        # failure must not overwrite it.
        partial=joinpath(temporary, "partial")
        mkpath(partial)
        _, plan=owner.campaign_plan(root)
        push!(plan["jobs"][1]["selections"],
            Dict("id"=>"unsupported_coaxial",
                "earth_impedance"=>"Carson1926", "earth_admittance"=>"default"))
        owner.write_campaign_state(joinpath(partial, "campaign.toml"), plan)
        @test !owner.resume_campaign(partial)
        row=only(owner.campaign_status(partial))
        @test row.state == "failed"
        @test row.completed == 1
        @test row.requested == 2
        saved=joinpath(partial, "two_bare_wires_coaxial", "0001.jld2")
        saved_bytes=read(saved)
        saved_mtime=stat(saved).mtime
        @test !isfile(joinpath(dirname(saved), "0002.jld2"))
        @test !owner.resume_campaign(partial)
        @test read(saved) == saved_bytes
        @test stat(saved).mtime == saved_mtime
        @test length(readdir(joinpath(dirname(saved), "attempts"))) == 2
        @test stat(artifact).mtime == modified
        state_path=joinpath(dirname(artifact), "state.toml")
        owner.write_campaign_state(state_path, Dict("state"=>"running", "completed"=>0))
        # A leftover lock file after a dead process does not block resumption.
        @test owner.campaign_status(root)[1].state == "interrupted"
        open(joinpath(root, "execution.lock"), "a+") do lock
            @test ccall(:flock, Cint, (Cint, Cint), Base.fd(lock), 6) == 0
            @test owner.campaign_status(root)[1].state == "running"
            @test_throws ArgumentError owner.resume_campaign(root)
        end
        @test owner.resume_campaign(root)
        @test read(artifact) == bytes
        # Changed declarations cannot replace completed numerical evidence.
        _, plan=owner.campaign_plan(root)
        plan["jobs"][1]["input_sha256"]="changed"
        owner.write_campaign_state(joinpath(root, "campaign.toml"), plan)
        @test !owner.resume_campaign(root)
        @test owner.campaign_status(root)[1].state == "failed"
        @test read(artifact) == bytes
        @test length(readdir(joinpath(dirname(artifact), "attempts"))) == 2
    end
end

@testitem "Gauntlet / explicit campaign frequency range survives checkpoint and resume" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using JLD2
    using LineCableModels
    using .GauntletSupport
    owner=GauntletSupport
    id=:cable_132kv_630mm2_flathor
    for invalid in ((0.01, 1e6), (50.0, 50.0), (1e6, 0.1), (0.1, Inf),
        (0.1,), [nothing, 1e6], :all)
        @test_throws ArgumentError owner.campaign_models([id], nothing; frequency_range = invalid)
    end
    selected=owner.campaign_models([id], nothing; frequency_range = (0.1, 1e6))[(
        id, false)]
    @test first(selected.nominal_problem.frequencies) == 1.0
    @test selected.problem.frequencies == owner._loggrid(0.1, 1e6, 101)
    uncertain=owner.campaign_models([:two_bare_wires],
        owner.RelativeStandardUncertainty(1.0; tags = (:geometry, :cable_layer));
        frequency_range = (0.1, 1e5))
    @test first(uncertain[(:two_bare_wires, true)].problem).frequencies ==
          owner._loggrid(0.1, 1e5, 101)
    mktempdir() do temporary
        root=joinpath(temporary, "common-band")
        @test owner.run_campaign(root, [id]; backends = (:coaxial,), catalogue = false,
            frequency_range = (0.1, 1e6))
        _, plan=owner.campaign_plan(root)
        @test plan["frequency_range"] == [0.1, 1e6]
        path=joinpath(root, "$(id)_coaxial", "0001.jld2")
        bytes=read(path)
        document=JLD2.load(path)
        @test document["frequencies"] == selected.problem.frequencies
        restored=LineCableModels.ImportExport.deserialize_value(document["problem"])
        @test restored.frequencies == document["frequencies"]
        @test owner.resume_campaign(root)
        @test read(path) == bytes
        @test first(selected.nominal_problem.frequencies) == 1.0
    end
end

@testitem "Gauntlet / unsupported all-bare PSCAD campaign never invokes a solver" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using LineCableModels
    using .GauntletSupport
    owner=GauntletSupport
    mktempdir() do temporary
        root=joinpath(temporary, "all-bare")
        @test owner.run_campaign(root, [:two_bare_wires]; backends = (:pscad,), catalogue = false)
        state=only(owner.campaign_status(root))
        @test state.state == "inapplicable"
        @test state.requested == state.completed == 0
        @test state.skipped == 1
        @test occursin("at least one insulated cable", state.message)
        @test owner.resume_campaign(root)
        @test !isfile(joinpath(root, "two_bare_wires_pscad", "0001.jld2"))
    end
end

@testitem "Gauntlet / manual selections use the complete formulation grammar" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using LineCableModels
    using .GauntletSupport
    owner=GauntletSupport
    choices=owner.parse_selections([
        "--select", "insulation_admittance=default,Ametani2004",
        "--select", "semicon_admittance=default,Ametani2004"])
    @test choices.insulation_admittance isa AbstractGrid
    paired=owner.campaign_selections(nothing, :coaxial, false; choices, combine = :zip)
    product=owner.campaign_selections(nothing, :coaxial, false; choices)
    @test length(paired.selections) == 2
    @test length(product.selections) == 4
    @test isempty(product.skipped)
    normalized=owner.reference_case(:cable_320kv_armoured_dc_bipole)
    @test normalized.problem.frequencies != normalized.nominal_problem.frequencies
    @test owner.campaign_input(normalized, :deterministic) ==
          owner.numerical_input_sha256(normalized.problem)
    @test owner.campaign_input(normalized, :deterministic) !=
          owner.numerical_input_sha256(normalized.nominal_problem)
    @test [(value.insulation_admittance, value.semicon_admittance)
           for value in paired.selections] ==
          [(:default, :default), (:Ametani2004, :Ametani2004)]
    @test Set((value.insulation_admittance, value.semicon_admittance)
    for value in product.selections) ==
          Set(Iterators.product((:default, :Ametani2004), (:default, :Ametani2004)))

    # This iterates the public contract, not a frozen count of slots or formulas.
    # Every backend must receive each requested slot even when its implementation
    # subsequently reports a fixed equation, an absence, or unsupported physics.
    withenv("LINECABLEMODELS_GETDP"=>"/unused/planning-only/getdp") do
        for name in keys(Formulation().definitions)
            axis=NamedTuple{(name,)}((Grid((:default, :default)),))
            plan=owner.campaign_selections(nothing, :coaxial, false; choices = axis)
            @test length(plan.selections) == 2
            @test all(value -> getproperty(value, name) === :default, plan.selections)
            record=Dict(string(key)=>string(value)
            for (key, value) in pairs(first(plan.selections)))
            for backend in (:coaxial, :fem, :pscad)
                requested=owner.campaign_formulation(backend, record, :Ametani2004)
                @test getproperty(requested.definitions, name) === :default
                @test requested.definitions.insulation_admittance === :default
            end
        end
    end
    alternatives=owner.parse_selections(["--select", "earth_impedance=default,Saad1996",
        "--select", "earth_properties=default,default",
        "--select", "pipe_impedance=default"])
    selected=owner.campaign_selections(
        nothing, :coaxial, false; choices = alternatives, combine = :zip)
    @test length(selected.selections) == 2
    @test selected.selections[2].earth_properties === :default
    @test all(value -> value.pipe_impedance === :default, selected.selections)
    override=owner.campaign_selections(nothing, :coaxial, false;
        choices = (insulation_admittance = :default,), dielectric = :Ametani2004)
    @test only(override.selections).insulation_admittance === :default
    @test only(override.selections).semicon_admittance === :Ametani2004

    @test_throws ArgumentError owner.campaign_selections(nothing, :coaxial, true; choices)
    @test_throws ArgumentError owner.campaign_selections(nothing, :coaxial, false;
        choices = (not_a_formula = Grid(:default),))
    @test_throws ArgumentError owner.campaign_selections(nothing, :coaxial, false;
        choices = (pipe_impedance = Grid(:InventedPipeAuthor2099),))
    @test_throws ArgumentError owner.campaign_selections(nothing, :coaxial, false;
        choices, combine = :invented)
    for arguments in (["--select"], ["--select", "earth_impedance"],
        ["--select", "earth_impedance=default,"],
        ["--select", "earth_impedance=default",
            "--select", "earth_impedance=Pollaczek1926"])
        @test_throws ArgumentError owner.parse_selections(arguments)
    end
end

@testitem "Gauntlet / normalized frequency snapshots replay the computed problem" tags=[:gauntlet] setup=[GauntletSupport] begin
    using LineCableModels, JLD2, SHA
    using .GauntletSupport
    owner=GauntletSupport
    mktempdir() do root
        directory=joinpath(root, "normalized")
        id=:cable_320kv_armoured_dc_bipole
        @test owner.run_campaign(directory, [id]; backends = (:coaxial,), catalogue = false)
        path=joinpath(directory, "$(id)_coaxial", "0001.jld2")
        bytes=read(path)
        record=JLD2.load(path)
        restored=LineCableModels.ImportExport.deserialize_value(record["problem"])
        @test restored.frequencies == record["frequencies"]
        @test first(restored.frequencies) == 0.1
        @test length(restored.frequencies) == 101
        @test owner.numerical_input_sha256(restored) == record["input_sha256"]
        selected=record["formulation"]
        replay=compute(restored, Formulation(; selected.definitions..., options = selected.options))
        @test replay.Z.values == record["Z"]
        @test replay.Y.values == record["Y"]
        @test owner.resume_campaign(directory)
        @test read(path) == bytes
    end
end

@testitem "Gauntlet / explicit formulation axes checkpoint without replacement" tags=[:gauntlet] setup=[GauntletSupport] begin
    using SHA
    using JLD2
    using LineCableModels
    using .GauntletSupport
    owner=GauntletSupport
    mktempdir() do temporary
        root=joinpath(temporary, "paired")
        choices=(earth_impedance = Grid((:default, :Saad1996)),
            insulation_admittance = Grid((:default, :Ametani2004)),
            semicon_admittance = Grid((:default, :Ametani2004)))
        @test owner.run_campaign(root, [:two_bare_wires]; backends = (:coaxial,),
            catalogue = false, choices, combine = :zip)
        @test only(owner.campaign_status(root)).completed == 2
        _, plan=owner.campaign_plan(root)
        @test plan["combine"] == "zip"
        @test all(value -> haskey(value, "pipe_impedance"), plan["jobs"][1]["selections"])
        files=[joinpath(root, "two_bare_wires_coaxial", lpad(index, 4, '0')*".jld2")
               for index in 1:2]
        records=JLD2.load.(files)
        @test records[1]["selection"]["insulation_admittance"] == "default"
        @test records[2]["selection"]["insulation_admittance"] == "Ametani2004"
        @test records[1]["problem"] == records[2]["problem"]
        @test records[1]["formulation"].definitions.insulation_admittance === :default
        @test records[2]["formulation"].definitions.insulation_admittance === :Ametani2004
        @test records[1]["computation_signature"] != records[2]["computation_signature"]
        @test records[1]["Z"] != records[2]["Z"]
        model=owner.reference_case(:two_bare_wires)
        expected=compute(model.nominal_problem,
            Formulation(
                earth_impedance = :Saad1996,
                insulation_admittance = :Ametani2004, semicon_admittance = :Ametani2004,
                options = (reduce_bundle = false, kron_reduction = false,
                    ideal_transposition = false, temperature_correction = true)))
        @test records[2]["Z"] == expected.Z.values
        @test records[2]["Y"] == expected.Y.values
        before=[(sha256(read(path)), stat(path).mtime) for path in files]
        @test owner.resume_campaign(root)
        @test [(sha256(read(path)), stat(path).mtime) for path in files] == before
    end
end

@testitem "Gauntlet / CI cannot start or resume a manual campaign" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using LineCableModels
    using .GauntletSupport
    mktempdir() do root
        destination=joinpath(root, "must_not_be_created")
        withenv("CI"=>"true") do
            @test_throws ArgumentError GauntletSupport.run_campaign(destination,
                [:two_bare_wires]; backends = (:coaxial,), catalogue = false)
            @test_throws ArgumentError GauntletSupport.resume_campaign(destination)
        end
        @test !ispath(destination)
    end
end

@testitem "Gauntlet / FEM checkpoints fingerprint executable bytes" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using LineCableModels
    using .GauntletSupport
    mktempdir() do root
        path=joinpath(root, "getdp")
        write(path, "first executable fixture")
        formulation=Formulation(:LineCableModelsFEM; fem_options = (getdp_executable = path,))
        first_record=GauntletSupport.campaign_implementation(formulation)
        write(path, "changed executable fixture")
        second_record=GauntletSupport.campaign_implementation(formulation)
        @test first_record.selection.executable.path ==
              second_record.selection.executable.path
        @test first_record.selection.executable.sha256 !=
              second_record.selection.executable.sha256
        @test first_record.selection_sha256 != second_record.selection_sha256
        @test !isempty(first_record.selection.gmsh_version)
    end
end
