using LineCableModels.Engine: admittivity, reductions

@testset "Gmsh field model" begin
    extension = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    FieldModel, Strand = extension.FieldModel, extension.Strand
    formulation = FEM()
    facts = Dict(
        "two_bare_wires" => (2, 1, 0, Symbol[], [false, false]),
        "two_insulated_wires" => (4, 2, 2,
            [:core_insulation, :core_insulation], [false, false]),
        "coaxial_circular_strands" => (85, 5, 5,
            [:core_semicon_inner, :core_insulation, :core_semicon_outer,
                :sheath_bedding, :armour_jacket], [false]),
        "coaxial_rectangular_strands" => (85, 5, 5,
            [:core_semicon_inner, :core_insulation, :core_semicon_outer,
                :sheath_bedding, :armour_jacket], [false]),
        "overhead_bare_wires" => (2, 1, 0, Symbol[], [true, true]),
        "mixed_bare_wires" => (2, 1, 0, Symbol[], [true, false]),
        "one_insulated_sector" => (3, 2, 0, Symbol[], [false]),
        "milliken_solid_segments" => (81, 5, 7,
            [:binder, :semicon_inner, :xlpe, :semicon_outer, :bedding_inner,
                :bedding_outer, :jacket], [false]),
        "three_insulated_sectors" => (39, 3, 1, [:jacket], [false]),
        "three_cores_steel_pipe" => (9, 4, 4,
            [:a_insulation, :b_insulation, :c_insulation, :jacket], [false]))
    # Per case, the copies of a repeat by tag and count: the regions whose construction is a
    # `Strand`.
    coaxial = [:core_strands_1 => 6, :core_strands_2 => 12, :core_strands_3 => 18,
        :armour_wires => 40]
    strands = Dict("coaxial_circular_strands" => coaxial,
        "coaxial_rectangular_strands" => coaxial,
        "milliken_solid_segments" => [:segment => 6, :screen => 65],
        "three_insulated_sectors" => [:neutral_wires => 30])
    @testset "$case" for (case, reference) in REFERENCES
        regions, groups, layers, expected_layers, air = facts[case]
        problem = LineParametersProblem(reference.system, reference.earth;
            temperature = reference.temperature,
            frequencies = reference.frequencies[1:10:end])
        model = FieldModel(problem, formulation)
        soil = problem.earth.layers[2]
        @test length(model.earth) == length(problem.frequencies)
        @test all(state == EquivalentEarthLayer(soil.rho, soil.eps_r, soil.mu_r)
        for state in model.earth)
        @test length(model.material) == length(problem.system.regions) == regions
        @test model.label == reductions(problem.system, FEM().options).label
        @test length(model.label) == length(problem.system.terminals)
        @test model.air == air
        @test length(model.materials) == groups
        @test count(model.layer) == layers
        @test [region.source.tag for region in problem.system.regions[findall(model.layer)]] ==
              expected_layers
        for (r, region) in enumerate(problem.system.regions)
            material = region.source.material
            group = model.materials[model.material[r]]
            @test group.admittivity ==
                  [admittivity(formulation, material, f, problem.temperature)
                   for f in problem.frequencies]
        end
        copies = findall(==(Strand()), model.construction)
        tags = [region.source.tag for region in problem.system.regions[copies]]
        @test [tag => count(==(tag), tags) for tag in unique(tags)] ==
              get(strands, case, Pair{Symbol, Int}[])
    end

    copper = Material(:conductor, Inf)
    design = build(CableDesign, "wire",
        Terminal(:core, Region(:core_metal, Circle(0.0425), copper)))
    system = build(CableSystem, [design, design],
        [Placement(0.0, -1.0), Placement(1.0, -1.0)];
        connections = [(core = 1,), (core = 2,)], id = "two_bare_wires", line_length = 1.0)
    problem = LineParametersProblem(system, homogeneous(rho = 0.1);
        frequencies = last(first(REFERENCES)).frequencies[1:10:end])
    @test_throws DomainError FieldModel(problem, formulation)

    reference = last(first(REFERENCES))
    system = build(CableSystem, reference.system.designs,
        [Placement(0.0, -support(reference.system.designs[1].outer)),
            reference.system.placements[2]];
        connections = [(core = 1,), (core = 2,)], id = "two_bare_wires", line_length = 1.0)
    problem = LineParametersProblem(system, reference.earth;
        frequencies = reference.frequencies[1:10:end])
    @test_throws DomainError FieldModel(problem, formulation)
    @testset "a layered earth reaches the model as its equivalent earth" begin
        air, bottom = reference.earth.layers
        layered = LineParametersProblem(reference.system,
            EarthModel((air, EarthLayer(100.0, 1.0, 1.0, 5.0), bottom));
            frequencies = reference.frequencies[1:10:end])
        @test all(
            layer -> layer === EquivalentEarthLayer(bottom.rho, bottom.eps_r, bottom.mu_r),
            FieldModel(layered, formulation).earth)
    end
    @testset "a negative equivalent resistivity warns and builds the model" begin
        earth = EarthModel((first(reference.earth.layers), EarthLayer(50.0, 4.0, 1.0, 1.0),
            EarthLayer(2000.0, 30.0, 1.0, 20.0), EarthLayer(10.0, 60.0, 1.0, 3.0),
            EarthLayer(300.0, 8.0, 1.0)))
        problem = LineParametersProblem(reference.system, earth; frequencies = [1.0e6])
        model = @test_logs (:warn, r"the equivalent earth of :xue2021") FieldModel(
            problem, FEM(equivalent_earth = :xue2021))
        @test only(model.earth).rho < 0
    end
    @testset "a lossy insulation reaches the model with the value of its law" begin
        polymer = Material(:insulator, 1.0e14, 2.3; tan_delta = 1.0e-3)
        metal = Region(:core_metal, Circle(0.0425), Material(:conductor, 1.7241e-8))
        design = build(CableDesign, "insulated_wire", Terminal(:core, metal),
            Region(:core_insulation, Offset(1.0e-3), polymer))
        system = build(CableSystem, [design, design],
            [Placement(0.0, -1.0), Placement(1.0, -1.0)];
            connections = [(core = 1,), (core = 2,)], id = "two_insulated_wires",
            line_length = 1.0)
        problem = LineParametersProblem(system, homogeneous(rho = 0.1);
            frequencies = reference.frequencies[1:10:end])
        lossy = FEM(insulation_admittance = :lossy)
        model = FieldModel(problem, lossy)
        r = findfirst(region -> region.source.tag == :core_insulation, problem.system.regions)
        @test model.materials[model.material[r]].admittivity ==
              [lossy.insulation_admittance(polymer, f) for f in problem.frequencies]
    end
    @testset "a magnetic earth and a scan outside a law warn" begin
        polymer = Material(:insulator, 1.0e14, 2.3; tan_delta = 1.0e-3)
        metal = Region(:core_metal, Circle(0.0425), Material(:conductor, 1.7241e-8))
        design = build(CableDesign, "insulated_wire", Terminal(:core, metal),
            Region(:core_insulation, Offset(1.0e-3), polymer))
        system = build(CableSystem, [design, design],
            [Placement(0.0, -1.0), Placement(1.0, -1.0)];
            connections = [(core = 1,), (core = 2,)], id = "two_insulated_wires",
            line_length = 1.0)
        earth = EarthModel((EarthLayer(Inf, 1.0, 1.0), EarthLayer(100.0, 10.0, 2.0, 5.0),
            EarthLayer(50.0, 5.0, 1.0)))
        problem = LineParametersProblem(system, earth; frequencies = [1.0, 1.0e3])
        formulation = FEM(earth_properties = :alipio2014, equivalent_earth = :xue2021)
        @test_logs (:warn,
            "the formula :alipio2014 holds for frequencies in 100.0 .. 4.0e6; " *
            "the computation has frequencies in 1.0 .. 1000.0") (:warn,
            "the formula :xue2021 holds for mu_r_earth in 1.0 .. 1.0; " *
            "the computation has mu_r_earth in 1.0 .. 2.0") LineCableModels.validate(
            LineCableModels.domain(problem, formulation.options), formulation)
    end
    @testset "the laws of the earth reach the model" begin
        polymer = Material(:insulator, 1.0e14, 2.3; tan_delta = 1.0e-3)
        metal = Region(:core_metal, Circle(0.0425), Material(:conductor, 1.7241e-8))
        design = build(CableDesign, "insulated_wire", Terminal(:core, metal),
            Region(:core_insulation, Offset(1.0e-3), polymer))
        system = build(CableSystem, [design, design],
            [Placement(0.0, -1.0), Placement(1.0, -1.0)];
            connections = [(core = 1,), (core = 2,)], id = "two_insulated_wires",
            line_length = 1.0)
        layers = (EarthLayer(100.0, 10.0, 1.0, 5.0), EarthLayer(500.0, 20.0, 1.0, 10.0),
            EarthLayer(50.0, 5.0, 1.0))
        problem = LineParametersProblem(
            system, EarthModel((EarthLayer(Inf, 1.0, 1.0), layers...));
            frequencies = [1.0e3, 1.0e6])
        formulation = FEM(earth_properties = :alipio2014, equivalent_earth = :xue2021)
        model = FieldModel(problem, formulation)
        @test model.earth == [formulation.equivalent_earth(
                   map(layer -> formulation.earth_properties(layer, f), layers), f)
               for f in problem.frequencies]
        @test first(model.earth).eps_r < 0
    end
end
