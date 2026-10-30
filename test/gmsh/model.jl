using LineCableModels.Engine: admittivity, reductions

@testset "Gmsh field model" begin
    FieldModel = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).FieldModel
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
        "mixed_bare_wires" => (2, 1, 0, Symbol[], [true, false]))
    @testset "$case" for (case, reference) in REFERENCES
        regions, groups, layers, expected_layers, air = facts[case]
        problem = LineParametersProblem(reference.system, reference.earth;
            temperature = reference.temperature,
            frequencies = reference.frequencies[1:10:end])
        model = FieldModel(problem, formulation)
        earth = problem.earth.layers[2]
        @test length(model.earth) == length(problem.frequencies)
        @test all(state == earth for state in model.earth)
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
            if material isa Material{<:Real, :conductor}
                @test real.(group.admittivity) ==
                      fill(
                    inv(formulation.resistivity(material, problem.temperature)),
                    length(problem.frequencies))
            else
                expected = [admittivity(CoaxialFormulation(), material, f)
                            for f in problem.frequencies]
                @test all(isapprox.(group.admittivity, expected; rtol = 1e-14))
            end
        end
        if startswith(case, "coaxial")
            core = findall(
                region -> startswith(string(region.source.tag), "core_strands_"),
                problem.system.regions)
            armour = findall(region -> region.source.tag == :armour_wires,
                problem.system.regions)
            @test findall(model.strand) == sort([core; armour])
            @test length(intersect(findall(model.strand), core)) == 36
            @test length(intersect(findall(model.strand), armour)) == 40
        else
            @test isempty(findall(model.strand))
        end
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
end
