@testitem "Engine / internal shunt / a coat of its layer's permittivity is the bare wire" tags=[:unit, :engine] begin
    using LinearAlgebra
    SM = LineCableModels.Engine.ShuntModel
    include(joinpath(pkgdir(LineCableModels), "test", "support", "internal_shunt.jl"))
    # 24 wires of metal radius 2 mm lie in a layer of permittivity 1. Coats of radius 2.5 mm
    # would cover 0.95 of their pitch circle.
    pitch = 24*5e-3/(2pi*0.95)
    ring(coat, epsilon) = [(x = pitch*cos(2pi*j/24), y = pitch*sin(2pi*j/24), r = 2e-3,
        coat, epsilon, terminal = 2) for j in 0:23]
    layer(wires) = annular_group(pitch - 3.5e-3, pitch + 4.5e-3, 1.0; wires)
    equal, bare = layer(ring(2.5e-3, 1.0)), layer(ring(2e-3, 1.0))
    boundary(wire) = SM.Formula(:boundary; options = (resolution = (wire = wire,),))
    capacitance(g, wire) = SM._shunt_capacitance(boundary(wire), g).C
    scale = maximum(abs, capacitance(bare, 64))
    # The discretizations differ, and they agree to the error of each.
    @test maximum(abs, capacitance(equal, 32)-capacitance(bare, 32)) < 1e-8scale
    @test maximum(abs, capacitance(equal, 64)-capacitance(bare, 64)) < 1e-12scale
    # A bare wire's coat operator vanishes, so its rows are the potential.
    plan = SM._shunt_plan(boundary(32), bare, SM._shunt_layout(bare, boundary(32).options.data.resolution))
    @test all(coat -> iszero(coat.operator) && iszero(coat.own), plan.coats)
    @test plan.layout.sectors == 24
end

@testitem "Engine / internal shunt / a coated wire on the axis has the capacitance of concentric layers" tags=[:unit, :engine] begin
    using LinearAlgebra
    SM = LineCableModels.Engine.ShuntModel
    include(joinpath(pkgdir(LineCableModels), "test", "support", "internal_shunt.jl"))
    epsilon0 = LineCableModels.Commons.vacuum_permittivity(Float64)
    a, c, ro, p = 2e-3, 2.5e-3, 0.01, 32
    angles = [2pi*j/(2p) for j in 0:(2p - 1)]
    placed = [(2pi*j/p, 0.75) for j in 0:(p - 1)]
    targets, sources = [c*cis(theta) for theta in angles], [f*c*cis(psi) for (psi, f) in placed]
    for (coat, host) in ((2.3, 1.0), (2.3, 4.0), (10.0, 1.0))
        wire = (x = 0.0, y = 0.0, r = a, coat = c, epsilon = coat, terminal = 2, layer = 1)
        operator, own = SM._shunt_coat(wire, host, angles, placed)
        # The Green function of a grounded disk of radius ro, by one image, replaces the
        # layered one. The coat's operator and table are those of the solver.
        G = [log(abs(ro^2-z*conj(s))/(ro*abs(z-s)))/host for z in targets, s in sources]
        C = 2pi*epsilon0*sum(qr(G+operator*G+own) \ ones(length(targets)))
        layers = layered_group([(ri = a, ro = c, epsilon = coat), (ri = c, ro, epsilon = host)])
        @test C ≈ 2pi*epsilon0/SM._shunt_kernel_coefficients(layers, 4).Rtotal rtol=1e-9
    end
end

@testitem "Engine / internal shunt / coated rings converge geometrically and stay rings" tags=[:unit, :engine] begin
    using LinearAlgebra
    SM = LineCableModels.Engine.ShuntModel
    include(joinpath(pkgdir(LineCableModels), "test", "support", "internal_shunt.jl"))
    boundary(wire) = SM.Formula(:boundary; options = (resolution = (wire = wire,),))
    # 24 wires of metal radius 2 mm with PE coats of radius 2.5 mm in a layer of permittivity 1.
    ring(pitch; layer = 1) = [(x = pitch*cos(2pi*j/24), y = pitch*sin(2pi*j/24), r = 2e-3,
        coat = 2.5e-3, epsilon = 2.3, terminal = 2, layer) for j in 0:23]
    spaced, touching = 24*5e-3/(2pi*0.6), 24*5e-3/(2pi)
    groups = (
        spaced = annular_group(spaced - 3.5e-3, spaced + 4.5e-3, 1.0; wires = ring(spaced)),
        # Each coat touches its two neighbours.
        touching = annular_group(touching - 3.5e-3, touching + 4.5e-3, 1.0;
            wires = ring(touching)),
        # Each coat touches the interface with a more permittive inner layer.
        interface = layered_group([(ri = spaced - 3.5e-3, ro = spaced - 2.5e-3, epsilon = 2.5),
            (ri = spaced - 2.5e-3, ro = spaced + 4.5e-3, epsilon = 1.0)];
            wires = ring(spaced; layer = 2)))
    for g in groups
        reference = SM._shunt_capacitance(boundary(128), g).C
        scale = maximum(abs, diag(reference))
        results = [SM._shunt_capacitance(boundary(wire), g) for wire in (16, 32, 64)]
        errors = [maximum(abs, result.C-reference)/scale for result in results]
        @test errors[3] < errors[2] < errors[1]
        @test errors[2] < 1e-6
        @test all(result -> result.diagnostic.sectors == 24, results)
        C = results[2].C
        @test abs(C[1, 2]-C[2, 1]) < 1e-8abs(C[2, 2])
    end
    # The ring's Fourier blocks give the dense solve of the same points.
    g = groups.touching
    selected = boundary(32)
    reduced = SM._shunt_capacitance(selected, g)
    dense = SM._shunt_plan(selected, g,
        merge(SM._shunt_layout(g, selected.options.data.resolution), (; sectors = 0, sector = 0)))
    buffers = SM.initialize_buffers(selected, Float64, g, dense, (;))
    unreduced = SM._shunt_capacitance!(selected, (; input = g, plan = dense, buffers))
    @test unreduced.diagnostic.sectors == 0
    @test unreduced.C ≈ reduced.C rtol=1e-12
end

@testitem "Engine / internal shunt / a coated design reaches the boundary shunt" tags=[:unit, :engine] setup=[CoatedArmour] begin
    SM = LineCableModels.Engine.ShuntModel
    DM = LineCableModels.DataModel
    design = CoatedArmour.armoured()
    groups = SM.dielectric_groups(SM.Formula(:boundary), design, DM.flatten(design),
        Formulation().formulas, 1).groups
    group = only(groups)
    @test all(wire -> wire.r ≈ 1.5e-3 && wire.coat ≈ 1.7e-3 && wire.epsilon ≈ 2.3, group.wires)
    # The coat radius and permittivity are descriptors, and they rebuild the group.
    values = Vector(group)
    @test length(values) == 3length(group.layers) + 5length(group.wires)
    @test SM.DielectricGroup(group, values).wires == group.wires
    coated = compute(CoatedArmour.problem(design), Formulation(shunt_model = :boundary))
    coaxial = compute(CoatedArmour.problem(design))
    @test all(isfinite, coated.Y)
    @test coated.Z == coaxial.Z
    # Coats touch but do not overlap, and a coat lies in its layer.
    layers = [(ri = 1.0, ro = 3.0, epsilon = 1.0)]
    wire(x, coat) = (x = x, y = 0.0, r = 0.2, coat, epsilon = 2.3, terminal = 2, layer = 1)
    touching = [wire(1.75, 0.25), wire(2.25, 0.25)]
    @test SM.validate(touching, SM.DielectricGroup, (cable_id = "coats",), layers, 3.0) ===
          touching
    @test_throws "overlapping coats" SM.validate([wire(1.75, 0.3), wire(2.25, 0.3)],
        SM.DielectricGroup, (cable_id = "coats",), layers, 3.0)
    @test_throws "a coat crosses a dielectric interface" SM.validate([wire(1.25, 0.3)],
        SM.DielectricGroup, (cable_id = "coats",), layers, 3.0)
    # Float32 descriptors solve in Float64.
    single = SM.DielectricGroup(group, Float32.(values))
    @test all(isfinite, SM.internal_capacitance(SM.Formula(:boundary), single).C)
    # A contextual Float32 ring has no fill. The boundary shunt refuses it with its error.
    contextual = CoatedArmour.problem(CoatedArmour.armoured(T = Float32, coat = 0.2f-3,
            coat_material = Material(:insulator, 1.0f14, 2.3f0), enclosed = false);
        frequencies = Float32[50])
    @test all(isfinite, compute(contextual).Y)
    @test_throws "neither a closed concentric tube" compute(contextual,
        Formulation(shunt_model = :boundary))
end

@testitem "ParametricBuilder / a Grid coat thickness computes through the boundary shunt" tags=[:unit, :parametric] setup=[CoatedArmour] begin
    designs = Gridspace{CableDesign}(coat -> CoatedArmour.armoured(; coat),
        (Grid((0.1e-3, 0.2e-3)),))
    results = compute(CoatedArmour.problem(designs), Formulation(shunt_model = :boundary))
    @test length(results) == 2
    @test all(result -> all(isfinite, result.Y), results)
    @test first(results).Y != last(results).Y
end

@testitem "Engine / Measurements / an uncertain coat permittivity reaches the boundary shunt" tags=[:unit, :measurements] setup=[CoatedArmour] begin
    using Measurements
    SM = LineCableModels.Engine.ShuntModel
    DM = LineCableModels.DataModel
    permittivity = measurement(2.3, 0.1)
    design = CoatedArmour.armoured(coat_material = Material(:insulator, 1.0e14, permittivity))
    result = compute(CoatedArmour.problem(design), Formulation(shunt_model = :boundary))
    @test uncertainty(imag(result.Y[3, 3, 1])) > 0
    @test all(isfinite, Measurements.value.(result.Y))
    # The implicit derivative through the coat's rows equals a central difference.
    group(design) = only(SM.dielectric_groups(SM.Formula(:boundary), design,
        DM.flatten(design), Formulation().formulas, 1).groups)
    capacitance(epsilon) = SM.internal_capacitance(SM.Formula(:boundary),
        group(CoatedArmour.armoured(coat_material = Material(:insulator, 1.0e14, epsilon)))).C
    propagated = SM.internal_capacitance(SM.Formula(:boundary), group(design)).C
    central = (capacitance(2.3001)-capacitance(2.2999))/0.0002
    @test Measurements.derivative.(propagated, Ref(permittivity)) ≈ central rtol=1e-3
end
