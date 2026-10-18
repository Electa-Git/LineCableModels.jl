@testitem "ParametricBuilder / builders / one grammar method lifts every builder" tags=[:unit, :parametric] begin
    const DM = LineCableModels.DataModel
    cu = Material(:conductor, 1.7241e-8, 1.0, 1.0)
    al = Material(:conductor, 2.82e-8, 1.0, 1.0)
    xlpe = Material(:insulator, 1.97e14, 2.3, 1.0)
    pe = Material(:insulator, 1e15, 2.25, 1.0)

    # Each builder is the singleton of a builder type, and the grammar's one call method
    # takes all its calls.
    for builder in (layers, assembly, terminal, solid, shell, core, insulation, screen,
            sheath, bedding, jacket, filler, wires, stranded, milliken, rope, armor, tape,
            cores, pipe, duct)
        @test builder isa LineCableModels.Builder
        @test length(methods(builder)) == 1
    end

    # A zip pairs the builder's own parameters, its keywords included.
    zipped = insulation(Grid((pe, xlpe)); t = Grid((5.5e-3, 6e-3)), combine = :zip)
    @test [(region.material.rho, region.primitive.t) for region in zipped] ==
          [(1e15, 5.5e-3), (1.97e14, 6e-3)]
    paired = core(Grid((cu, al)); r = Grid((1e-3, 2e-3)), combine = :zip)
    @test [(region.material.rho, region.primitive.r) for region in paired] ==
          [(cu.rho, 1e-3), (al.rho, 2e-3)]
    member = terminal(:core, core(cu; r = 1e-3))
    ring = cores(member; n = 3, r = Grid((1e-2, 2e-2)), φ0 = Grid((0.0, 0.3)),
        names = (:a, :b, :c), combine = :zip)
    @test [(assembly.pattern.r, assembly.pattern.φ0) for assembly in ring] ==
          [(1e-2, 0.0), (2e-2, 0.3)]
    placed = at(member, Grid((0.1, 0.2)), Grid((0.3, 0.4)); combine = :zip)
    @test [(placement.at.x, placement.at.y) for placement in placed] == [(0.1, 0.3), (0.2, 0.4)]

    # A product varies the first source fastest, in the order the call writes them.
    lays = Grid((DM.LayRatio(10), DM.LayRatio(12)))
    courses = wires(cu; shape = Disk(1e-3), n = 6, lay = lays, r = Grid((3e-3, 4e-3)))
    @test [(group.path.lay.q, group.pattern.r) for group in courses] ==
          [(10, 3e-3), (12, 3e-3), (10, 4e-3), (12, 4e-3)]

    # The scalar call is the scalar body, and the calls infer.
    @test @inferred(insulation(xlpe; t = 1e-3)) == build(insulation, xlpe; t = 1e-3)
    @test @inferred(cores(member; n = 3, r = 1e-2, names = (:a, :b, :c))) isa Assembly
    @test @inferred(first(insulation(Grid((pe, xlpe)); t = 1e-3))) isa Region
    @test @inferred(at(member, Grid((0.1, 0.2)), 0.0)) isa Gridspace{DM.AssemblyMember}
end

@testitem "ParametricBuilder / builders / a build method extends a builder for an own type" tags=[:unit, :parametric] begin
    struct OwnInsulator
        rho::Float64
        eps_r::Float64
    end
    LineCableModels.build(::typeof(insulation), m::OwnInsulator; t, tag = :insulation) =
        Region(tag, Shell(t), Material(:insulator, m.rho, m.eps_r))
    mine, yours = OwnInsulator(1e15, 2.25), OwnInsulator(1.97e14, 2.3)

    # The scalar body receives scalars, and the grammar lifts every call over it.
    @test insulation(mine; t = 1e-3) isa Region
    paired = insulation(Grid((mine, yours)); t = Grid((5.5e-3, 6e-3)), combine = :zip)
    @test [(region.material.rho, region.primitive.t) for region in paired] ==
          [(1e15, 5.5e-3), (1.97e14, 6e-3)]
end
