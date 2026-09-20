@testitem "Gmsh FEM / resume requires effective inputs and preserves completed runs" tags=[:extension] setup=[FormulaContractModels] begin
    using Gmsh
    using LineCableModels
    extension = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    artifact = withenv("LINECABLEMODELS_GETDP"=>nothing) do
        extension._getdp_selection(computation_options(LineCableModelsFEM, ComputationOptions((;))))
    end
    @test artifact.source === :artifact
    @test artifact.artifact_hash == string(extension.artifact_hash(
        "getdp", extension.GETDP_ARTIFACTS_TOML,
    ))
    @test isfile(artifact.path)
    copper = Material(kind=:conductor, rho=1.72e-8)
    dielectric = Material(kind=:insulator, rho=1.0e8, eps_r=2.3, tan_delta=0.025)
    design = build(CableDesign, "fem-resume-inputs", Stack(
        Group(:core, Region(:metal, Disk(0.005), copper)),
        Region(:insulation, Shell(0.005), dielectric)))
    system = build(LineCableSystem, design, (0.0, -0.1);
        connections = Dict(:core=>1), system_id = "fem-resume-inputs")
    problem = LineParametersProblem(system; frequencies = [50.0, 1000.0],
        earth_props = LineCableModels.Earth.EarthModel(100.0, 10.0, 1.0))
    formulation = Formulation(:LineCableModelsFEM;
        options=(ideal_transposition=false,))
    formulation_controls = (getdp_executable=artifact.path, gmsh_verbosity=0,)
    model = extension._resolved_fem_model(problem, formulation)
    inputs = extension._fem_input_record(model, formulation, computation_options(LineCableModelsFEM, ComputationOptions(formulation_controls)))
    @test inputs.schema_version == 7
    @test inputs.getdp_selection.source === :explicit
    @test inputs.getdp_selection.artifact_hash === nothing
    @test isfile(inputs.getdp_selection.path)
    @test !hasproperty(inputs.getdp_identity, :path)
    @test occursin(r"Version\s*:\s*3\.5\.0", inputs.getdp_identity.info)
    @test occursin("PETSc", inputs.getdp_identity.info)
    @test occursin("complex arithmetic", inputs.getdp_identity.info)
    @test_throws LineCableModelsFEMError extension._getdp_selection(
        computation_options(LineCableModelsFEM, ComputationOptions((getdp_executable=joinpath(tempdir(), "missing-getdp"),))),
    )
    withenv("LINECABLEMODELS_GETDP"=>joinpath(tempdir(), "missing-getdp")) do
        @test_throws LineCableModelsFEMError extension._getdp_selection(
            computation_options(LineCableModelsFEM, ComputationOptions((;))),
        )
    end
    @test inputs.mesh_fingerprint == extension._mesh_fingerprint(model, Gmsh.gmsh.GMSH_API_VERSION)
    smaller_controls = computation_options(LineCableModelsFEM, ComputationOptions((;formulation_controls..., domain_skin_depths=1.5)))
    smaller_model = extension._resolved_fem_model(problem, formulation, smaller_controls)
    smaller_inputs = extension._fem_input_record(smaller_model, formulation, smaller_controls)
    @test smaller_inputs.mesh_fingerprint != inputs.mesh_fingerprint
    # Dictionary iteration order can change between Julia versions. Resume
    # records the Julia version; the key must be repeatable within that runtime.
    recorded_key=extension._mesh_fingerprint(model, "recorded-gmsh-version")
    @test extension._mesh_fingerprint(deepcopy(model), "recorded-gmsh-version") == recorded_key
    @test extension._mesh_fingerprint(model, "different-gmsh-version") != recorded_key
    @test inputs.adapter_sources isa NamedTuple
    @test haskey(inputs.adapter_sources, Symbol("geometry.jl"))
    @test !haskey(inputs.adapter_sources, Symbol("formulations.jl"))
    other = Formulation(:LineCableModelsFEM; earth_properties = nothing,
        options = formulation.options)
    other_controls = formulation_controls
    other_model = extension._resolved_fem_model(problem, other)
    other_inputs = extension._fem_input_record(other_model, other, computation_options(LineCableModelsFEM, ComputationOptions(other_controls)))
    lossy = Formulation(:LineCableModelsFEM; insulation_admittance = :lossy,
        options = formulation.options)
    lossy_controls = formulation_controls
    lossy_model = extension._resolved_fem_model(problem, lossy)
    lossy_inputs = extension._fem_input_record(lossy_model, lossy, computation_options(LineCableModelsFEM, ComputationOptions(lossy_controls)))
    @test extension.JSON3.write(inputs) == extension.JSON3.write(other_inputs)
    @test extension.JSON3.write(inputs) != extension.JSON3.write(lossy_inputs)
    temperature_law=FormulaContractModels.ScaledResistivity(1.06)
    soil_law=FormulaContractModels.ScaledSoil(epsilon=2)
    thermal = LineCableModelsFEM(
        temperature_dependence=temperature_law,
        options=formulation.options)
    thermal_controls = formulation_controls
    hot_problem = LineParametersProblem(system;temperature=80.0,
        frequencies=problem.frequencies,earth_props=problem.earth_props)
    hot_model = extension._resolved_fem_model(hot_problem,thermal)
    hot_inputs = extension._fem_input_record(hot_model,thermal, computation_options(LineCableModelsFEM, ComputationOptions(thermal_controls)))
    @test first(hot_inputs.materials).sigma ≈ first(inputs.materials).sigma ./ 1.06
    dispersive = LineCableModelsFEM(
        earth_properties=soil_law,
        options=formulation.options)
    dispersive_controls = formulation_controls
    dispersive_model = extension._resolved_fem_model(problem,dispersive)
    dispersive_inputs = extension._fem_input_record(dispersive_model,dispersive, computation_options(LineCableModelsFEM, ComputationOptions(dispersive_controls)))
    @test dispersive_inputs.mesh_fingerprint == inputs.mesh_fingerprint
    @test getproperty.(dispersive_inputs.earth_materials,:eps_r) ==
        2 .* getproperty.(inputs.earth_materials,:eps_r)
    # Simulate a changed resolver with the declaration and mesh-size policy
    # held fixed. The domains handed to Gmsh are authoritative for reuse.
    changed_model = deepcopy(model)
    region = first(changed_model.region_plans)
    changed_model.region_plans[1] = extension.FEMRegionPlan(region.object_id,
        region.cable_index, region.region_index, region.terminal_index,
        region.material_index, Disk(region.shape.r * 0.9, region.shape.at), region.mesh_size)
    changed_inputs = extension._fem_input_record(changed_model, formulation, computation_options(LineCableModelsFEM, ComputationOptions(formulation_controls)))
    @test changed_inputs.mesh_fingerprint != inputs.mesh_fingerprint
    @test changed_inputs.materials == inputs.materials
    @test changed_inputs.region_mesh_sizes == inputs.region_mesh_sizes
    @test LineCableModels.ImportExport.serialize_value(changed_model.problem) ==
          LineCableModels.ImportExport.serialize_value(model.problem)
    ownership = deepcopy(model)
    ownership.region_plans[1] = extension.FEMRegionPlan(region.object_id,
        region.cable_index, region.region_index, 0, region.material_index,
        region.shape, region.mesh_size)
    @test extension._fem_input_record(ownership, formulation, computation_options(LineCableModelsFEM, ComputationOptions(formulation_controls))).mesh_fingerprint !=
          inputs.mesh_fingerprint
    @test_throws ArgumentError compute(problem, LineCableModelsFEM[])
    law=FormulaContractModels.ScaledSoil(rho=Inf)
    unsupported = Formulation(:LineCableModelsFEM;
        earth_properties = law,
        options = formulation.options)
    unsupported_controls = formulation_controls
    before = Bool(Gmsh.gmsh.is_initialized())
    @test_throws LineCableModelsFEMError compute(problem, [formulation, unsupported])
    @test Bool(Gmsh.gmsh.is_initialized()) == before
    mktempdir() do root
        run = extension._create_run(root)
        extension._prepare_run_inputs!(run, model)
        assets = extension._getdp_assets(joinpath(run.path, "input", "getdp"))
        @test keys(assets) ==
              (:model, :jacobian, :integration, :materials, :quasi_tem, :quasi_full)
        @test map(basename, values(assets)) ==
              ("model.pro", "jacobian.pro", "integration.pro", "materials.pro",
                  "quasi-tem.pro", "quasi-full.pro")
        captured = map(path -> (read(path), stat(path).mtime), assets)
        extension._prepare_run_inputs!(run, model)
        @test map(path -> (read(path), stat(path).mtime), assets) == captured
        # A changed equation file must not be silently repaired or reused.
        write(assets.quasi_tem, "// changed equation snapshot\n")
        @test_throws LineCableModelsFEMError extension._prepare_run_inputs!(run, model)
        @test read(assets.quasi_tem, String) == "// changed equation snapshot\n"
        write(assets.quasi_tem, captured.quasi_tem[1])
        @test !extension._resume_inputs_match(run.path, model, inputs)
        extension._write_json_atomic(joinpath(run.path, "input", "computation.json"), inputs)
        @test extension._resume_inputs_match(run.path, model, other_inputs)
        @test !extension._resume_inputs_match(run.path,hot_model,hot_inputs)
        @test !extension._resume_inputs_match(run.path,dispersive_model,dispersive_inputs)
        # Current missing metadata is corruption; obsolete development schema
        # examples do not establish a compatibility obligation.
        incomplete=Dict(String(k)=>v for (k,v) in pairs(extension.JSON3.read(extension.JSON3.write(inputs))))
        delete!(incomplete,"solver_protocol")
        extension._write_json_atomic(joinpath(run.path,"input","computation.json"),incomplete)
        @test !extension._resume_inputs_match(run.path,model,inputs)
        extension._write_json_atomic(joinpath(run.path,"input","computation.json"),inputs)
        @test !extension._resume_inputs_match(run.path, changed_model, changed_inputs)
        @test !extension._resume_inputs_match(run.path, lossy_model, lossy_inputs)
        @test_throws ArgumentError extension._resume_run(root, run.path, lossy_model, lossy_inputs)
        resumed = extension._resume_run(root, run.path, model, inputs)
        @test resumed.path == run.path
        extension._transition!(run, extension.completed, "test completion")
        preserved = read(joinpath(run.path, "run.json"))
        @test !extension._resume_inputs_match(run.path, model, inputs)
        @test_throws ArgumentError extension._resume_run(root, run.path, model, inputs)
        @test read(joinpath(run.path, "run.json")) == preserved
        next_run = extension._resume_run(root, :latest, model, inputs)
        @test next_run.path != run.path
        @test read(joinpath(run.path, "run.json")) == preserved
        # Completed results can be read only when their numerical files are
        # protected. This does not transition or rewrite their run state.
        for name in ("Z.tsv", "P.tsv", "scan_complete.tsv")
            write(joinpath(run.path, "raw", name), "checksum fixture $name")
        end
        scan = extension.FEMScan(zeros(ComplexF64, 1, 1, 2), zeros(ComplexF64, 1, 1, 2), String[])
        extension._write_scan_checksums(run, scan)
        retained_inputs = merge(inputs, (; getdp_identity=(sha256="fixture", info="fixture")))
        extension._write_json_atomic(joinpath(run.path, "input", "computation.json"), retained_inputs)
        @test extension._resume_inputs_match(run.path, model, retained_inputs)
        @test !extension._resume_inputs_match(
            run.path, model, merge(retained_inputs, (; getdp_identity = nothing)))
        retained = extension._resume_run(root, run.path, model, retained_inputs)
        @test retained.state === extension.completed
        @test retained.path == run.path
        @test read(joinpath(run.path, "run.json")) == preserved
        @test extension._check_scan_checksums(retained, scan) === nothing
        write(joinpath(run.path, "raw", "Z.tsv"), "corrupted numerical payload")
        @test_throws LineCableModelsFEMError extension._check_scan_checksums(retained, scan)
        write(joinpath(run.path, "raw", "checksums.json"), "broken JSON")
        @test_throws LineCableModelsFEMError extension._check_scan_checksums(retained, scan)
        changed = merge(retained_inputs, (;
            adapter_sources = Dict("geometry.jl"=>"different implementation")))
        @test !extension._resume_inputs_match(run.path, model, changed)
        @test !extension._resume_inputs_match(run.path, model, merge(retained_inputs, (;
            owned_gmsh = false)))
        for execution in
            (merge(inputs.execution, (ui=true,)),
            merge(inputs.execution, (mesh_policy=:remesh,)))
            excluded = merge(retained_inputs, (; execution))
            extension._write_json_atomic(joinpath(run.path, "input", "computation.json"), excluded)
            @test !extension._resume_inputs_match(run.path, model, excluded)
        end
        if Sys.isunix()
            executable = joinpath(root, "getdp-identity")
            write(executable, "#!/bin/sh\necho 'GetDP Version 3.6.0 fixture A'\n")
            chmod(executable, 0o700)
            configured = Formulation(:LineCableModelsFEM; options=formulation.options)
            configured_controls = (getdp_executable=executable, gmsh_verbosity=0,)
            first_record = extension._fem_input_record(model, configured, computation_options(LineCableModelsFEM, ComputationOptions(configured_controls)))
            write(executable, "#!/bin/sh\necho 'GetDP Version 3.6.0 fixture B'\n")
            second_record = extension._fem_input_record(model, configured, computation_options(LineCableModelsFEM, ComputationOptions(configured_controls)))
            @test first_record.getdp_selection.path == second_record.getdp_selection.path
            @test first_record.getdp_selection.source === :explicit
            @test first_record.getdp_identity.sha256 != second_record.getdp_identity.sha256
            @test first_record.getdp_identity.info != second_record.getdp_identity.info

            environment = joinpath(root, "getdp-environment")
            cp(executable, environment)
            chmod(environment, 0o700)
            relocated = Formulation(:LineCableModelsFEM; options=formulation.options)
            relocated_controls = (getdp_executable=environment, gmsh_verbosity=0,)
            relocated_inputs = extension._fem_input_record(model, relocated, computation_options(LineCableModelsFEM, ComputationOptions(relocated_controls)))
            @test relocated_inputs.getdp_identity == second_record.getdp_identity
            @test relocated_inputs.getdp_selection.path != second_record.getdp_selection.path
            relocation_run = extension._create_run(root)
            extension._prepare_run_inputs!(relocation_run, model)
            extension._write_json_atomic(
                joinpath(relocation_run.path, "input", "computation.json"), second_record)
            @test extension._resume_inputs_match(relocation_run.path, model, relocated_inputs)
            @test extension._resolve_getdp(computation_options(LineCableModelsFEM, ComputationOptions(relocated_controls)), relocation_run) == realpath(environment)
            withenv("LINECABLEMODELS_GETDP"=>environment) do
                selected = extension._getdp_selection(computation_options(LineCableModelsFEM, ComputationOptions((;))))
                @test selected.source === :environment
                @test selected.path == realpath(environment)
                explicit = extension._getdp_selection(computation_options(LineCableModelsFEM, ComputationOptions(configured_controls)))
                @test explicit.source === :explicit
                @test explicit.path == realpath(executable)
            end
            write(environment, "#!/bin/sh\necho 'GetDP Version 3.6.0 changed binary'\n")
            @test_throws LineCableModelsFEMError extension._resolve_getdp(computation_options(LineCableModelsFEM, ComputationOptions(relocated_controls)), relocation_run)
            @test !extension._resume_inputs_match(relocation_run.path, model,
                extension._fem_input_record(model, relocated, computation_options(LineCableModelsFEM, ComputationOptions(relocated_controls))))
        end
    end
end
