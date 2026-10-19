@testitem "PSCAD / transient settings, fitting and canvas XML" tags=[:integration, :pscad] setup=[TestFixtures] begin
    using LineCableModels
    using EzXML
    const P = LineCableModels.PSCAD

    parameters(node) = Dict(p["name"] => p["value"] for p in findall("./paramlist/param", node))
    design = TestFixtures.coaxial_design()
    system = build(LineCableSystem, [design, design], [(-0.08, -1.0), (0.08, -1.0)];
        connections = [Dict(:core => 1, :sheath => 2), Dict(:core => 3, :sheath => 4)],
        system_id = "transient-settings", line_length = 600.0)
    earth = homogeneous(rho = 100.0)
    formulation = Formulation()
    blueprints = P._pscad_blueprints(system)
    settings = (duration = 0.02, step = 2e-6, sample = 10e-6, output_stem = "channels")
    fit = (FS = 1.0, FE = 2e6, Numf = 200, YMaxP = 24, AMaxP = 26,
        YMaxE = 0.03, AMaxE = 0.04, DCenab = 1, DCCOR = 2, CPASS = 2, MaxRPtol = 3e6,
        ET_PE = 1e-6, MER_PE = 10.0, MIT_PE = 7, FSP = 0.01, FEP = 2e6,
        NFP = 500, FDIS = 1, Interp1 = 0)
    directory = mkpath(joinpath(@__DIR__, "..", "..", "..",
        ".linecablemodels", "test", "transient-settings"))
    mktempdir(directory) do destination
        path = joinpath(destination, "transient.pscx")
        document = P._pscad_project(system, earth, 50.0, blueprints, formulation, nothing;
            settings, fit)
        main = only(findall("/project/definitions/Definition[@name='Main']/schematic", document))
        canvas = P.PSCADCanvas(main)
        cells = [P.allocate_cell!(canvas) for _ in 1:7]
        @test cells == [(72, 72), (504, 72), (936, 72), (1368, 72),
            (1800, 72), (2232, 72), (72, 288)]
        interface = cells[1] .+ (144, 54)
        component = P.component!(canvas, "cable_interface", interface, ["NCAB" => "2"])
        @test (component["x"], component["y"]) == ("216", "126")
        # PowerImpedance.jl `src/PSCAD/emitters/cables.jl` line 137: interface ports.
        for side in (-1, 1), port in 1:2
            point = interface .+ (side * P._PSCAD_GRID, P._PSCAD_INTERFACE_PORT_SPACING * (port - 1))
            endpoint = point .+ (side * 36, 0)
            P.node_label!(canvas, "port-$(side)-$port", endpoint)
            P.wire!(canvas, point, endpoint)
        end
        for position in ((-18, 18), (18, -18), (14400, 18), (18, 14400))
            @test_throws r"coordinates must lie in \[0, 14382\]" P.component!(canvas, "ground", position)
            @test_throws r"coordinates must lie in \[0, 14382\]" P.wire!(canvas, (18, 18), position)
            @test_throws r"coordinates must lie in \[0, 14382\]" P.wire!(canvas, position, (18, 18))
        end
        @test P.validate((0, 14382), P.PSCADCanvas) == (0, 14382)
        @test_throws ArgumentError P.component!(canvas, "ground", (1, 18))
        @test_throws ArgumentError P.node_label!(canvas, "invalid", (18, 1))
        @test_throws ArgumentError P.wire!(canvas, (0, 0), (18, 18))
        @test_throws ArgumentError P.wire!(canvas, (0, 0), (0, 0))
        @test_throws ArgumentError P.wire!(canvas, (0, 0), (0, 1))
        write(path, document)
        parsed = readxml(path)
        @test all(parameters(c)["size"] == "5" for c in findall("//schematic", parsed))
        values = parameters(root(parsed))
        @test values["time_duration"] == "0.02"
        @test values["time_step"] == "2.0"
        @test values["sample_step"] == "10.0"
        @test values["PlotType"] == "1"
        @test values["output_filename"] == "channels.out"
        interface_component = only(findall("//User[@defn='master:cable_interface']", parsed))
        @test (interface_component["x"], interface_component["y"]) == ("216", "126")
        @test parameters(interface_component)["NCAB"] == "2"
        phase = only(findall("//User[@defn='master:Line_FrePhase_Options']", parsed))
        @test all(parse(Float64, parameters(phase)[string(key)]) == value for (key, value) in pairs(fit))
        @test isempty(findall("//User[@defn='master:Line_FreDep_Options']", parsed))
        for cable in findall("//User[@defn='master:Cable_Coax']", parsed)
            @test parameters(cable)["elim1"] == "0"
            @test parameters(cable)["CROSSBOND"] == "0"
            @test parameters(cable)["LC"] == "3"
        end
        labels = findall("//User[@defn='master:nodelabel']", parsed)
        @test [(p["x"], p["y"], parameters(p)["Name"]) for p in labels] ==
            [("162", "126", "port--1-1"), ("162", "180", "port--1-2"),
                ("270", "126", "port-1-1"), ("270", "180", "port-1-2")]
        wires = findall("//Wire[@classid='WireOrthogonal']", parsed)
        @test [(w["x"], w["y"], w["w"], w["h"]) for w in wires] ==
            [("198", "126", "36", "0"), ("198", "180", "36", "0"),
                ("234", "126", "36", "0"), ("234", "180", "36", "0")]
        @test [[(v["x"], v["y"]) for v in findall("./vertex", w)] for w in wires] ==
            [[("0", "0"), ("-36", "0")], [("0", "0"), ("-36", "0")],
                [("0", "0"), ("36", "0")], [("0", "0"), ("36", "0")]]
        identifiers = [node["id"] for node in findall("//*[@id]", parsed)]
        @test allunique(identifiers)

        defaults = P._pscad_project(system, earth, 50.0, blueprints, formulation, nothing; settings)
        write(path, defaults)
        phase = only(findall("//User[@defn='master:Line_FrePhase_Options']", readxml(path)))
        @test parameters(phase) == Dict("FS" => "0.5", "FE" => "1.0E6", "Numf" => "100",
            "DCCOR" => "1", "CPASS" => "0", "enablf" => "1", "shntcab" => "0.0",
            "YMaxP" => "20", "AMaxP" => "20", "YMaxE" => "0.2", "AMaxE" => "0.2",
            "DCenab" => "0", "MaxRPtol" => "2000000.0",
            "ET_PE" => "1.0E-10", "MER_PE" => "2", "MIT_PE" => "5",
            "FSP" => "0.001", "FEP" => "1.0E6", "NFP" => "1000",
            "FDIS" => "3", "Interp1" => "1")

        mode = P._pscad_project(system, earth, 50.0, blueprints, formulation, nothing;
            settings, model = :fd_mode)
        write(path, mode)
        parsed = readxml(path)
        @test isempty(findall("//User[@defn='master:Line_FrePhase_Options']", parsed))
        mode_component = only(findall("//User[@defn='master:Line_FreDep_Options']", parsed))
        @test parameters(mode_component) == Dict("FS" => "0.5", "FE" => "1.0E6", "FT" => "2000",
            "ZMaxP" => "20", "AMaxP" => "20", "ZMaxE" => "1.0", "AMaxE" => "1.0",
            "Interp1" => "1", "Output" => "0", "W1" => "1.0", "W2" => "1000.0",
            "W3" => "1.0", "Inflen" => "0")
        mode = P._pscad_project(system, earth, 50.0, blueprints, formulation, nothing;
            model = :fd_mode, fit = (FT = 500.0, ZMaxP = 30, ZMaxE = 0.5))
        write(path, mode)
        mode_component = only(findall("//User[@defn='master:Line_FreDep_Options']", readxml(path)))
        @test parameters(mode_component)["FT"] == "500.0"
        @test parse(Float64, parameters(mode_component)["ZMaxP"]) == 30
        @test parameters(mode_component)["ZMaxE"] == "0.5"

        line = LineParametersProblem(system; earth_props = earth, frequencies = [1.0, 1e6])
        problem = TransientAnalysisProblem(line;
            sending = (Step("cable:1:core", 1000.0; t0 = 2e-6),), receiving = (),
            window = (tmax = 20e-6, step = 1e-7))
        for (model, native_fit) in ((:fd_phase, fit), (:fd_mode, (FT = 500.0, ZMaxP = 30)))
            native = P.PSCADTransientFormulation(; model, fit = native_fit)
            inputs = P._pscad_inputs(problem, native, destination)
            controls = inputs.native_settings["frequency"]
            @test Set(keys(controls)) == Set(string.(keys(getproperty(P._PSCAD_FIT_FIELDS, model))))
            @test controls["Interp1"]["value"] isa Int
            @test controls["Interp1"]["readback"] == (model == :fd_phase ? "No" : "Yes")
            if model == :fd_phase
                @test controls["DCenab"] == Dict("value" => 1, "readback" => "Enabled")
                @test controls["DCCOR"] == Dict("value" => 2, "readback" => "Add Pole/Residue")
                @test controls["CPASS"] == Dict("value" => 2, "readback" => "Enforcement")
                @test controls["Numf"] == Dict("value" => 200, "readback" => "200")
            else
                @test controls["Output"] == Dict("value" => 0, "readback" => "No")
            end
            if model == :fd_phase
                @test controls["MER_PE"] == Dict("value" => 10.0, "readback" => 10.0)
                @test controls["MER_PE"]["readback"] isa Float64
            end
        end
        @test_throws ArgumentError P._pscad_project(system, earth, 50.0, blueprints, formulation,
            nothing; model = :unknown)
        @test_throws ArgumentError P._pscad_project(system, earth, 50.0, blueprints, formulation,
            nothing; fit = (unknown = 1,))
        @test_throws ArgumentError P._pscad_project(system, earth, 50.0, blueprints, formulation,
            nothing; model = :fd_mode, fit = (YMaxP = 1,))
    end
end

@testitem "PSCAD / fitting fields follow master-library definitions" tags=[:integration, :pscad] begin
    using LineCableModels.PSCAD: PSCADTransientFormulation
    const P = LineCableModels.PSCAD

    # Independent defaults from the two component definitions in m7-master.pslx.
    phase = (FS = 0.01, FE = 1e6, Numf = 100, DCCOR = 1, CPASS = 0,
        enablf = 1, shntcab = 1e-9, YMaxP = 20, AMaxP = 20, YMaxE = 0.2,
        AMaxE = 0.2, DCenab = 0, MaxRPtol = 2e6, ET_PE = 1e-10,
        MER_PE = 2, MIT_PE = 5, FSP = 0.001, FEP = 1e6,
        NFP = 1000, FDIS = 3, Interp1 = 1)
    mode = (FS = 0.5, FE = 1e6, FT = 2000, ZMaxP = 20, AMaxP = 20,
        ZMaxE = 1, AMaxE = 1, Interp1 = 1, Output = 0,
        W1 = 1, W2 = 1000, W3 = 1, Inflen = 0)
    for (model, defaults) in ((:fd_phase, phase), (:fd_mode, mode))
        @test keys(getproperty(P._PSCAD_FIT_FIELDS, model)) == keys(defaults)
        @test PSCADTransientFormulation(; model, fit = defaults).fit == defaults
        for (field, value) in pairs(defaults)
            fit = NamedTuple{(field,)}((value,))
            @test PSCADTransientFormulation(; model, fit).fit == fit
        end
    end
    @test_throws ArgumentError("PSCAD fd_phase fitting field Numf=300 must be one of (100, 200, 500, 1000)") PSCADTransientFormulation(fit = (Numf = 300,))
    @test_throws ArgumentError("PSCAD fd_phase fitting field YMaxP=60 must be Integer in [5, 50]") PSCADTransientFormulation(fit = (YMaxP = 60,))
    for fit in ((YMaxP = 20.5,), (MIT_PE = 5.5,), (NFP = 300,),
            (FDIS = 0,), (Interp1 = 0.5,), (DCCOR = 0,), (FS = 0,),
            (MER_PE = 101,), (shntcab = 0,), (FE = Inf,))
        @test_throws ArgumentError PSCADTransientFormulation(; fit)
    end
    @test_throws MethodError PSCADTransientFormulation(fit = (FS = "0.5",))
    @test_throws MethodError PSCADTransientFormulation(fit = (FS = 1 + 0im,))
    @test_throws r"ZMaxP=20.*allowed fields:.*FS" PSCADTransientFormulation(fit = (ZMaxP = 20,))
    @test_throws r"YMaxP=20.*allowed fields:.*FS" PSCADTransientFormulation(model = :fd_mode, fit = (YMaxP = 20,))
    @test_throws r"AMaxP=60.*\[5, 50\]" PSCADTransientFormulation(fit = (AMaxP = 60,))
    @test_throws r"ZMaxP=101.*\[5, 100\]" PSCADTransientFormulation(model = :fd_mode, fit = (ZMaxP = 101,))
    @test PSCADTransientFormulation(model = :fd_mode, fit = (AMaxP = 60, FS = 0)).fit == (AMaxP = 60, FS = 0)
    @test PSCADTransientFormulation(fit = (ET_PE = -1e-6, YMaxP = 20.0)).fit == (ET_PE = -1e-6, YMaxP = 20.0)
end

@testitem "PSCAD / constants project byte preservation" tags=[:integration, :pscad] setup=[TestFixtures] begin
    using LineCableModels
    using EzXML
    const P = LineCableModels.PSCAD

    # Load the unchanged project writer at the reference revision as an independent reference.
    repository = normpath(joinpath(@__DIR__, "..", "..", ".."))
    source = read(`git -C $repository show 2657ccf2:src/pscad/importexport/project.jl`, String)
    reference = Module(:PSCADConstantsReference)
    Base.include_string(reference, """
        using LineCableModels
        using LineCableModels.DataModel: LineCableSystem
        using LineCableModels.Earth: EarthModel
        import LineCableModels.DataModel
        import LineCableModels: constitutive
        using LineCableModels.DataModel.BaseParams: equivalent_eps
        using EzXML: XMLDocument, ElementNode, setroot!, addelement!
        using LineCableModels.PSCAD: _pscad_paramlist!, _pscad_attributes!,
            _pscad_canvas_params!, _pscad_wire!, _pscad_user!, _pscad_value,
            _PSCAD_FREQUENCY_BINDING, _PSCAD_CABLE_BINDING, _PSCAD_GROUND_BINDING
        """ * source)
    system = TestFixtures.three_phase_system()
    earth = homogeneous(rho = 100.0)
    formulation = Formulation()
    blueprints = P._pscad_blueprints(system)
    before = Base.invokelatest(reference._pscad_project,
        system, earth, 50.0, blueprints, formulation, nothing)
    after = P._pscad_project(system, earth, 50.0, blueprints, formulation, nothing)
    directory = mkpath(joinpath(repository, ".linecablemodels", "test", "transient-settings"))
    mktempdir(directory) do destination
        old_path, new_path = joinpath(destination, "before.pscx"), joinpath(destination, "after.pscx")
        for parameter in findall("//schematic/paramlist/param[@name='size']", before)
            @test parameter["value"] in ("0", "5")
            parameter["value"] = "5"
        end
        write(old_path, before)
        write(new_path, after)
        @test read(new_path) == read(old_path)
        output = export_data(:pscad, system, earth; file_name = joinpath(destination, "export.pscx"))
        @test read(output) == read(old_path)
    end
end
