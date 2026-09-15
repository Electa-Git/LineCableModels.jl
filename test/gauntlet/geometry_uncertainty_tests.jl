@testitem "Gauntlet / joint parameter ingestion samples each block once" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport.Gauntlet
    using Random, Measurements, Serialization
    joint = Gridspace{NamedTuple{(:core_radius,:insulation_thickness)}}(
        s -> (core_radius=0.01s,insulation_thickness=0.002s), (Grid(1.,10.),))
    record = (law=:fixture, independent_inputs=(scale=(nominal=1.,std=0.1,
        support=(1-sqrt(3)*0.1,1+sqrt(3)*0.1)),),
        dependencies=(core_radius=(:scale,),insulation_thickness=(:scale,)))
    variation = JointParameterGrids(joint;record)
    model = load_case(:two_insulated_wires;variation=compose_variations(
        variation,ParameterGrids(frequencies=Grid(([50.],[60.])))))
    @test length(model.problem) == 2
    @test [only(p.frequencies) for p in model.problem] == [50.,60.]
    calls = Ref(0)
    sampler = (rng,mu,sigma) -> (calls[]+=1;mu+sigma*(2rand(rng)-1))
    sampled = rand(Xoshiro(5),model.problem;distribution=sampler)
    @test calls[] == 1
    @test outer_radius(first(sampled.system.designs)) > 0
    measured = only(joint)
    @test Measurements.cov(measured.core_radius,measured.insulation_thickness) ≈ 2e-7
    inputs = correlation_record(model)
    @test inputs.rule === :joint_builders
    @test isempty(inputs.uncertain_primitives)
    @test inputs.joint_inputs == [record]
    @test variation_record(variation).assumptions == record
    exact = load_case(:two_bare_wires;variation=ExactOverrides(frequencies=[50.]))
    @test isempty(correlation_record(exact).uncertain_primitives)
    for later in (ExactOverrides(core_radius=0.02),
            ParameterGrids(core_radius=Grid((0.01,0.02))),variation,
            RelativeStandardUncertainty(10.;tags=(:geometry,:cable_layer)))
        @test_throws ArgumentError load_case(:two_insulated_wires;
            variation=compose_variations(variation,later))
    end
    unknown = JointParameterGrids(Gridspace{NamedTuple{(:unknown,)}}(
        x -> (unknown=x,), (Grid(1.,10.),));record=(law=:test,))
    @test_throws ArgumentError load_case(:two_insulated_wires;variation=unknown)
    @test_throws ArgumentError JointParameterGrids(Gridspace{NamedTuple}(x->(x=x,),
        (Grid(1.,10.),));record=(law=:test,))
    @test_throws ArgumentError JointParameterGrids(Gridspace{NamedTuple{()}}(
        ()->(;),());record=(law=:test,))
    io = IOBuffer(); serialize(io,model.problem); seekstart(io)
    restored = deserialize(io)
    a = rand(Xoshiro(99),model.problem;distribution=:uniform)
    b = rand(Xoshiro(99),restored;distribution=:uniform)
    @test LineCableModels.ImportExport.serialize_value(a.system) ==
        LineCableModels.ImportExport.serialize_value(b.system)
    @test a.frequencies == b.frequencies
    zipped = Gridspace{NamedTuple{(:core_radius,:insulation_thickness)}}(
        (r,t)->(core_radius=r,insulation_thickness=t),
        (Grid((0.01,0.02)),Grid((0.002,0.003)));combine=:zip)
    zip_model = load_case(:two_insulated_wires;variation=compose_variations(
        JointParameterGrids(zipped;record=(law=:deterministic_fixture,)),
        ParameterGrids(frequencies=Grid(([50.],[60.])))))
    @test length(zip_model.problem) == 4
    @test [only(p.frequencies) for p in zip_model.problem] == [50.,60.,50.,60.]
    @test [outer_radius(first(p.system.designs)) for p in zip_model.problem] ≈
        [0.012,0.012,0.023,0.023]
    earth = JointParameterGrids(Gridspace{NamedTuple{(:earth_rho,)}}(
        rho->(earth_rho=rho,), (Grid(100.,5.),));record=(law=:independent_earth,))
    two_blocks = load_case(:two_insulated_wires;variation=compose_variations(
        variation,earth,ExactOverrides(frequencies=[50.])))
    calls[] = 0
    rand(Xoshiro(5),two_blocks.problem;distribution=sampler)
    @test calls[] == 2
    @test length(correlation_record(two_blocks).joint_inputs) == 2
end

@testitem "Gauntlet / failed subsea cases retain 512 feasible draws without bedding repair" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport.Gauntlet
    using Random
    const G = GauntletSupport.Gauntlet
    for id in (:cable_525kv_subsea_armoured_ac_flat,:cable_525kv_subsea_armoured_dc_bipole)
        def = Base.include(G,G.case_index()[id])
        variation = G._catalogue_geometry_variation(def,ExactOverrides(frequencies=[50.]))
        model = load_case(id;variation)
        joint = only(filter(v->v isa JointParameterGrids,G._variation_leaves(variation)))
        for seed in 1:512
            # Replaying the same RNG from the joint source gives the physical
            # dimensions used by this one-point model's only uncertain block.
            inputs = rand(Xoshiro(seed),joint.source;distribution=:uniform)
            problem = rand(Xoshiro(seed),model.problem;distribution=:uniform)
            design = first(problem.system.designs)
            bedding = only(filter(r->r.source.tag===:sheath_bedding,design.geometry.regions))
            wires = filter(r->r.source.tag===:armor_wires,design.geometry.regions)
            @test thickness(bedding.primitive) ≈ inputs.bedding_thickness
            @test length(wires) == 68
            @test outer_radius(design)/inputs.core_diameter ≈
                0.0805/def.parameters.core_diameter.nominal
            first_wire, second_wire = wires[1:2]
            chord = hypot((centroid(first_wire).-centroid(second_wire))...)
            @test chord > inputs.armor_wire_diameter
        end
    end
end

@testitem "Gauntlet / bounded catalogue law and deterministic base variations" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport.Gauntlet
    using Measurements, Random
    const G = GauntletSupport.Gauntlet
    definition = Base.include(G,G.case_index()[:two_insulated_wires])
    variation = G._catalogue_geometry_variation(definition,compose_variations(
        ParameterGrids(core_radius=Grid((0.01,0.02))),
        ExactOverrides(insulation_thickness=0.003,frequencies=[50.])))
    joint = only(filter(v->v isa JointParameterGrids,G._variation_leaves(variation)))
    @test length(joint.source) == 2
    @test joint.record.law === :bounded_correlated_geometry_v1
    @test joint.record.roles == (core_radius=:length,insulation_thickness=:length)
    @test collect(keys(joint.record.independent_inputs)) == [:length_scale]
    for (p,radius) in zip(joint.source,(0.01,0.02))
        @test nominal(p.core_radius) == radius
        @test nominal(p.insulation_thickness) == 0.003
        @test uncertainty(p.core_radius) ≈ 0.1radius
        @test uncertainty(p.insulation_thickness) ≈ 0.0003
        @test Measurements.cov(p.core_radius,p.insulation_thickness) ≈ 0.01radius*0.003
        @test uncertainty(p.core_radius/p.insulation_thickness) < 1e-14
    end
    model = load_case(definition.id;variation)
    @test length(model.problem) == 2
    @test [nominal(outer_radius(first(p.system.designs))) for p in model.problem] ≈
        [0.013,0.023]
    @test_throws ArgumentError G._catalogue_geometry_variation(definition,
        ParameterGrids(core_radius=Grid(0.01,10.)))
    @test_throws ArgumentError G._catalogue_geometry_variation(definition,
        ExactOverrides(core_radius=measurement(0.01,0.001)))
    parameters = merge(definition.parameters,(core_radius=case_parameter(:core_radius,
        0.01;tags=(:geometry,:cable_layer)),))
    invalid = case_definition(definition.build,definition.id,parameters,definition.port_order)
    @test_throws ArgumentError G._catalogue_geometry_variation(invalid,NoVariation())
    certificate = G._catalogue_sector_support((core_sectors=6,core_fillet_factor=0.04,
        core_outer_radius=0.03115,core_wire_radius=0.001475))
    @test certificate.strands_per_sector == 61
    @test 60 < first(certificate.capacity_bounds) < last(certificate.capacity_bounds) < 90
    @test last(certificate.fillet_support) < 1/3
    @test_throws ArgumentError G._catalogue_sector_support((core_sectors=6,
        core_fillet_factor=0.3,core_outer_radius=0.03115,core_wire_radius=0.001475))
end

@testitem "Gauntlet / all catalogue joint supports and seeded construction" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport.Gauntlet
    using Measurements, Random
    const G = GauntletSupport.Gauntlet
    const DM = LineCableModels.DataModel
    for id in G.CATALOGUE_CASE_IDS
        @testset "$id" begin
            def = Base.include(G,G.case_index()[id])
            variation = G._catalogue_geometry_variation(def,ExactOverrides(frequencies=[50.]))
            joint = only(filter(v->v isa JointParameterGrids,G._variation_leaves(variation)))
            model = load_case(id;variation)
            reference_inventory = length(first(model.nominal_problem.system.designs).geometry.regions)
            @test joint.record.law === :bounded_correlated_geometry_v1
            @test length(model.problem) == 1
            for (name,role) in pairs(joint.record.roles)
                @test role in (:length,:area,:dimensionless)
            end
            # Uniform endpoint construction and seeded interior draws do not
            # substitute for the retained analytic sector support certificate.
            for q in (-sqrt(3),0.,sqrt(3))
                p = rand(Xoshiro(1),model.problem;distribution=(rng,mu,sigma)->mu+q*sigma)
                @test length(first(p.system.designs).geometry.regions) == reference_inventory
                @test all(r->nominal(area(r)) >= 0,first(p.system.designs).geometry.regions)
            end
            measured = only(model.problem)
            @test uncertainty(outer_radius(first(measured.system.designs))) > 0
            for seed in 1:4
                p = rand(Xoshiro(seed),model.problem;distribution=:uniform)
                @test length(first(p.system.designs).geometry.regions) == reference_inventory
            end
        end
    end
end
