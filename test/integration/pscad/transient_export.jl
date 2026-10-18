@testitem "PSCAD / transient export and native terminal wiring" tags=[:integration, :pscad] setup=[TestFixtures] begin
    using LineCableModels
    using EzXML
    using LineCableModels.PSCAD: PSCADTransientFormulation
    using LineCableModels.Commons: vacuum_permeability, vacuum_permittivity
    const P = LineCableModels.PSCAD

    parameters(node) = Dict(p["name"] => p["value"] for p in findall("./paramlist/param", node))
    components(document, name) = findall("//User[@defn='master:$name']", document)
    design = TestFixtures.coaxial_design()
    system = build(LineCableSystem, [design, design], [(-0.08, -1.0), (0.08, -1.0)];
        connections = [Dict(:core => 1, :sheath => 2), Dict(:core => 3, :sheath => 4)],
        system_id = "transient-export", line_length = 600.0)
    earth = homogeneous(rho = 100.0)
    line = LineParametersProblem(system; earth_props = earth, frequencies = [1.0, 1e6])
    window = (tmax = 0.02, step = 2e-6)
    sending = (Step("cable:1:core", 1000.0; t0 = 0.001),
        TerminalImpedance("cable:1:sheath"; R = 10),
        TerminalImpedance("cable:2:sheath"; R = 10))
    problem = TransientAnalysisProblem(line; sending, receiving = (), window)
    fit = (FS = 1.0, FE = 2e6, Numf = 200, YMaxP = 24, AMaxP = 26,
        YMaxE = 0.03, AMaxE = 0.04, DCenab = 1, DCCOR = 2, CPASS = 1, MaxRPtol = 3e6,
        enablf = 0, shntcab = 1e-8)
    formulation = PSCADTransientFormulation(; fit)
    directory = mkpath(joinpath(@__DIR__, "..", "..", "..",
        ".linecablemodels", "test", "transient-export"))

    # Reconstruct the electrical and signal stubs independently from written XML.
    function connections(document)
        schematic = only(findall("/project/definitions/Definition[@name='Main']/schematic", document))
        labels = Dict((parse(Int, c["x"]), parse(Int, c["y"])) =>
            (c["defn"] == "master:ground" ? "ground" : parameters(c)["Name"])
            for c in findall("./User", schematic)
            if c["defn"] in ("master:nodelabel", "master:datalabel", "master:ground"))
        stubs = Dict{Tuple{Int, Int}, String}()
        for wire in findall("./Wire[@classid='WireOrthogonal']", schematic)
            origin = (parse(Int, wire["x"]), parse(Int, wire["y"]))
            a, b = [origin .+ (parse(Int, v["x"]), parse(Int, v["y"]))
                    for v in findall("./vertex", wire)]
            if haskey(labels, a) && haskey(labels, b)
                @test (labels[a] == "ground") != (labels[b] == "ground")
                label = labels[a] == "ground" ? b : a
                @test length(findall("./User[@defn='master:nodelabel'][@x='$(label[1])'][@y='$(label[2])']", schematic)) == 1
            else
                @test xor(haskey(labels, a), haskey(labels, b))
            end
            port, label = haskey(labels, a) ? (b, a) : (a, b)
            @test !haskey(stubs, port)
            stubs[port] = labels[label]
        end
        return stubs
    end

    mktempdir(directory) do destination
        file_name = joinpath(destination, "transient.pscx")
        document = readxml(export_data(:pscad, problem, formulation; file_name))
        @test all(parameters(c)["size"] == "5" for c in findall("//schematic", document))
        interfaces = components(document, "cable_interface")
        @test length(interfaces) == 2
        @test all(c -> all(parameters(c)[field] == "0" for field in ("GroundS", "GroundA", "GroundOC")), interfaces)
        @test all(c -> parameters(c)["NCAB"] == "2" &&
            parameters(c)["C1T"] == parameters(c)["C2T"] == "2", interfaces)
        @test [parameters(c)["send_recv"] for c in interfaces] == ["1", "2"]
        @test Set(parameters(c)["Name"] for c in components(document, "pgb")) ==
            union(Set("$end_name:$terminal:v" for end_name in ("sending", "receiving")
                    for terminal in terminal_labels(system)),
                Set("sending:$(a.terminal):i" for a in sending))
        @test all(c -> parameters(c)["Scale"] == "1000", components(document, "pgb"))
        @test length(components(document, "voltmetergnd")) == 8
        @test length(components(document, "ammeter")) == length(sending)
        @test [parse(Float64, parameters(c)["R"]) for c in components(document, "resistor")] == [10, 10]
        settings = parameters(root(document))
        @test parse(Float64, settings["time_duration"]) == window.tmax
        @test parse(Float64, settings["time_step"]) == window.step * 1e6
        @test parse(Float64, settings["sample_step"]) == window.step * 1e6
        @test settings["PlotType"] == "1"
        @test settings["output_filename"] == "transient.out"
        phase = parameters(only(components(document, "Line_FrePhase_Options")))
        @test all(parse(Float64, phase[string(key)]) == value for (key, value) in pairs(fit))
        source = parameters(only(components(document, "source_1")))
        @test source["AC"] == "0"
        @test parse(Float64, source["Vm"]) == 1
        breaker = parameters(only(components(document, "breaker1")))
        @test breaker["BOpen"] == "1"
        comparator = parameters(only(components(document, "compare")))
        @test parse(Float64, comparator["X"]) == 0.001
        @test (comparator["OL"], comparator["OH"]) == ("1", "0")
        stubs = connections(document)
        for (interface, end_name) in zip(interfaces, ("sending", "receiving"))
            x, y = parse(Int, interface["x"]), parse(Int, interface["y"])
            @test stubs[(x - 18, y)] == "$(end_name)_c1_core"
            @test stubs[(x + 18, y + 18)] == "$(end_name)_c1_sheath"
            @test stubs[(x - 18, y + 54)] == "$(end_name)_c2_core"
            @test stubs[(x + 18, y + 72)] == "$(end_name)_c2_sheath"
        end
        for meter in components(document, "ammeter")
            x, y = parse(Int, meter["x"]), parse(Int, meter["y"])
            node = replace(parameters(meter)["Name"], r"_i$" => "")
            @test stubs[(x, y)] == node
            @test stubs[(x + 36, y)] == node * "_branch"
        end
        for meter in components(document, "voltmetergnd")
            x, y = parse(Int, meter["x"]), parse(Int, meter["y"])
            @test stubs[(x, y)] == replace(parameters(meter)["Name"], r"_v$" => "")
        end
        # Voltage probes connect to existing terminal nodes.
        electrical_labels = Set(parameters(c)["Name"] for c in components(document, "nodelabel"))
        @test electrical_labels == union(
            Set("$(end_name)_c$(t.cable)_$(t.terminal)" for end_name in ("sending", "receiving")
                for t in system.terminal_order),
            Set("sending_c$(t.cable)_$(t.terminal)_branch" for t in system.terminal_order[[1, 2, 4]]),
            Set(["sending_c1_core_branch_dc"]))
        @test allunique(n["id"] for n in findall("//*[@id]", document))
        main = only(findall("/project/definitions/Definition[@name='Main']/schematic", document))
        @test allunique((c["x"], c["y"]) for c in findall("./User", main))

        mode = readxml(export_data(:pscad, problem,
            PSCADTransientFormulation(model = :fd_mode, fit = (FT = 500.0, ZMaxP = 30),
                sample = 10e-6, output_stem = "mode_channels"); file_name))
        @test isempty(components(mode, "Line_FrePhase_Options"))
        @test parse(Float64, parameters(only(components(mode, "Line_FreDep_Options")))["FT"]) == 500
        @test parse(Float64, parameters(root(mode))["sample_step"]) == 10
        @test parameters(root(mode))["output_filename"] == "mode_channels.out"

        tmax = 0.01922048123456789
        precise_window = (tmax = tmax, step = tmax / 4000)
        precise_problem = TransientAnalysisProblem(line; sending, receiving = (),
            window = precise_window)
        for sample in (nothing, 3 * precise_window.step)
            document = readxml(export_data(:pscad, precise_problem,
                PSCADTransientFormulation(; sample); file_name))
            settings = parameters(root(document))
            @test parse(Float64, settings["time_duration"]) == precise_problem.window.tmax
            @test parse(Float64, settings["time_step"]) == precise_problem.window.step * 1e6
            @test parse(Float64, settings["sample_step"]) ==
                something(sample, precise_problem.window.step) * 1e6
        end

        waveforms = (Cosine("cable:1:core", -2000, 50; phase = pi / 6),
            DoubleExponential("cable:2:core", 3000, 1e3, 1e5),
            TerminalImpedance("cable:1:sheath"; L = 0.01),
            TerminalImpedance("cable:2:sheath"; R = 3, L = 0.02))
        waves = TransientAnalysisProblem(line; sending = waveforms,
            receiving = (TerminalImpedance("cable:1:core"),), window)
        document = readxml(export_data(:pscad, waves, PSCADTransientFormulation(); file_name))
        sources = parameters.(components(document, "source_1"))
        ac = only(filter(p -> p["AC"] == "1", sources))
        @test parse(Float64, ac["Vm"]) ≈ 2 / sqrt(2)
        @test parse(Float64, ac["f"]) == 50
        @test parse(Float64, ac["Ph"]) ≈ 120
        @test length(components(document, "surge1")) == 1
        surge = parameters(only(components(document, "surge1")))
        @test surge["SW"] == "-1"
        @test parse(Float64, surge["ALPHA"]) == 0.001
        @test parse(Float64, surge["BETA"]) ≈ 0.1
        peak_time = log(100) / 99000
        @test parse(Float64, surge["Ipk"]) ≈ 3(exp(-1000peak_time) - exp(-100000peak_time))
        @test length(components(document, "inductor")) == 2
        @test length(components(document, "ammeter")) == 5
        stubs = connections(document)
        resistor = only(components(document, "resistor"))
        x, y = parse(Int, resistor["x"]), parse(Int, resistor["y"])
        @test stubs[(x + 36, y)] == "sending_c2_sheath_branch_rl"

        matched = TransientAnalysisProblem(line; sending, receiving = Matched(), window)
        document = readxml(export_data(:pscad, matched,
            PSCADTransientFormulation(; fit, frequency = 60.0); file_name))
        interfaces = components(document, "cable_interface")
        @test length(interfaces) == 4
        @test length(components(document, "Cable_Coax")) == 2
        @test length(components(document, "Line_FrePhase_Options")) == 1
        wrappers = findall("//Wire[@classid='Cable']", document)
        @test length(wrappers) == 2
        @test wrappers[1]["defn"] == wrappers[2]["defn"]
        @test all(w -> parse(Float64, parameters(only(findall("./User", w)))["Freq"]) == 60.0, wrappers)
        extension = parameters(only(findall("./User", wrappers[2])))
        @test extension["Name"] == "receiving_extension"
        c = 1 / sqrt(vacuum_permeability(Float64) * vacuum_permittivity(Float64))
        @test parse(Float64, extension["Length"]) * 1000 ≈ 1.1 * c * window.tmax / 2
        @test all(c -> parameters(c)["Name"] == "receiving_extension", interfaces[3:4])
        stubs = connections(document)
        for (original, near) in zip(interfaces[2:2], interfaces[3:3])
            a = (parse(Int, original["x"]), parse(Int, original["y"]))
            b = (parse(Int, near["x"]), parse(Int, near["y"]))
            for offset in ((-18, 0), (18, 18), (-18, 54), (18, 72))
                @test stubs[a .+ offset] * "_branch" == stubs[b .+ offset]
            end
        end
        @test length(components(document, "ammeter")) == 7
        @test length(components(document, "voltmetergnd")) == 8
        @test length(components(document, "resistor")) == 2
        @test length(components(document, "source_1")) == 1
        @test Set(parameters(c)["Name"] for c in components(document, "pgb")) ==
            union(Set("$end_name:$terminal:v" for end_name in ("sending", "receiving")
                    for terminal in terminal_labels(system)),
                Set("sending:$(a.terminal):i" for a in sending),
                Set("receiving:$terminal:i" for terminal in terminal_labels(system)))
        for meter in components(document, "ammeter")
            x, y = parse(Int, meter["x"]), parse(Int, meter["y"])
            node = replace(parameters(meter)["Name"], r"_i$" => "")
            @test stubs[(x, y)] == node
            @test stubs[(x + 36, y)] == node * "_branch"
        end
        mode_matched = readxml(export_data(:pscad, matched,
            PSCADTransientFormulation(model = :fd_mode, fit = (FT = 700.0,)); file_name))
        @test length(components(mode_matched, "Line_FreDep_Options")) == 1
        @test length(findall("//Wire[@classid='Cable']", mode_matched)) == 2

        shifted = LineParametersProblem(system; earth_props = earth, frequencies = [1.0, 1e6], sigma = 1.0)
        shifted_problem = TransientAnalysisProblem(shifted; sending, receiving = (), window)
        @test_throws ArgumentError export_data(:pscad, shifted_problem, formulation; file_name)
        bundled = build(LineCableSystem, [design, design], [(-0.08, -1.0), (0.08, -1.0)];
            connections = [Dict(:core => 1, :sheath => 2), Dict(:core => 1, :sheath => 3)],
            system_id = "bundled", line_length = 600.0)
        bundled_problem = TransientAnalysisProblem(LineParametersProblem(bundled; earth_props = earth,
            frequencies = [1.0]); sending = (), receiving = (), window)
        @test_throws ArgumentError export_data(:pscad, bundled_problem, formulation; file_name)
        @test_throws ArgumentError PSCADTransientFormulation(model = :unknown)
        @test PSCADTransientFormulation().frequency == 50.0
        for frequency in (0.0, -1.0, Inf, -Inf, NaN)
            @test_throws ArgumentError PSCADTransientFormulation(; frequency)
        end
        @test_throws ArgumentError PSCADTransientFormulation(sample = 0)
        @test_throws ArgumentError PSCADTransientFormulation(sample = Inf)
        @test_throws ArgumentError PSCADTransientFormulation(output_stem = "a/b")
        @test_throws ArgumentError PSCADTransientFormulation(fit = (FS = NaN,))
        @test_throws ArgumentError export_data(:pscad, problem,
            PSCADTransientFormulation(fit = (unknown = 1,)); file_name)
    end
end

@testitem "PSCAD / transient cable dielectric selection" tags=[:integration, :pscad] begin
    using EzXML
    using LineCableModels.PSCAD: PSCADTransientFormulation
    const P = LineCableModels.PSCAD
    defaults = @inferred PSCADTransientFormulation()
    @test defaults.line isa LineCableModels.Engine.LineParametersFormulation
    @test formula_id(defaults.line.formulas.insulation_admittance) === :lossless
    @test formula_id(defaults.line.formulas.semicon_admittance) === :lossless
    copper = Material(:conductor, 1.72e-8)
    dielectric = Material(:insulator, 1e10, 2.3; tan_delta = 0.03)
    semicon = Material(:semicon, 1e4, 40; tan_delta = 0.02)
    design = build(CableDesign, "transient-dielectric", terminal(:core,
        solid(copper, Disk(0.004)), screen(semicon; t = 0.0005),
        insulation(dielectric; t = 0.002)))
    system = build(LineCableSystem, [design], [Pose2(0, -1)];
        connections = [Dict(:core => 1)], line_length = 600.0)
    earth = homogeneous(rho = 100.0)
    line = LineParametersProblem(system; earth_props = earth, frequencies = [1.0, 1e6],
        temperature = 60.0)
    problem = TransientAnalysisProblem(line; sending = (), receiving = (),
        window = (tmax = 20e-6, step = 1e-7))
    lossy = Formulation(insulation_admittance = :lossy, semicon_admittance = :lossy)
    selected = @inferred PSCADTransientFormulation(line = lossy, frequency = 60.0)
    @test selected.line === lossy
    @test P.validate(selected) === selected
    @test P.computation_details(selected).data.line == NamedTuple(lossy)
    unsupported = merge(lossy.formulas, (insulation_admittance = lossy.formulas.earth_admittance,))
    invalid = LineCableModels.Engine.LineParametersFormulation(unsupported, lossy.options, lossy.definitions)
    @test_throws r"PSCAD does not implement insulation_admittance" PSCADTransientFormulation(line = invalid)
    mktempdir() do directory
        for selection in (defaults.line, lossy,
                CableConstantsFormulation(insulation_admittance = :lossy, semicon_admittance = :lossy),
                Formulation(:pscad; insulation_admittance = :lossy, semicon_admittance = :lossy))
            transient = readxml(export_data(:pscad, problem,
                PSCADTransientFormulation(line = selection, frequency = 60.0);
                file_name = joinpath(directory, "transient.pscx")))
            constants = readxml(export_data(:pscad, system, earth; formulation = selection,
                base_freq = 60.0, temperature = 60.0, file_name = joinpath(directory, "constants.pscx")))
            losses = [p["value"] for p in findall("//User[@defn='master:Cable_Coax']/paramlist/param[starts-with(@name,'LT')]", transient)]
            @test losses == [p["value"] for p in findall("//User[@defn='master:Cable_Coax']/paramlist/param[starts-with(@name,'LT')]", constants)]
            @test selection === defaults.line ? all(iszero, parse.(Float64, losses)) :
                parse(Float64, first(losses)) > 0
        end
    end
end
