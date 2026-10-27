using LineCableModels.Engine: admittivity, resistivity, reductions

@testset "Gmsh field model" begin
    FieldModel = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).FieldModel
    formulation = FEM()
    @testset "$case" for ((case, reference), regions, groups, layers, expected_layers) in
                         zip(REFERENCES, (2, 4, 85, 85), (1, 2, 5, 5), (0, 2, 5, 5),
        (Symbol[], [:core_insulation, :core_insulation],
            [:core_semicon_inner, :core_insulation, :core_semicon_outer,
                :sheath_bedding, :armour_jacket],
            [:core_semicon_inner, :core_insulation, :core_semicon_outer,
                :sheath_bedding, :armour_jacket]))
        problem = LineParametersProblem(reference.system, reference.earth;
            temperature = reference.temperature,
            frequencies = reference.frequencies[1:10:end])
        model = FieldModel(problem, formulation)
        earth = problem.earth.layers[2]
        @test length(model.earth) == length(problem.frequencies)
        @test all(state == (rho = earth.rho, eps_r = earth.eps_r, mu_r = earth.mu_r)
        for state in model.earth)
        @test length(model.material) == length(problem.system.regions) == regions
        @test model.label == reductions(problem.system, FEM().options).label
        @test length(model.label) == length(problem.system.terminals)
        @test all(!, model.air)
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
                    inv(resistivity(formulation.resistivity, material,
                        problem.temperature)),
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
