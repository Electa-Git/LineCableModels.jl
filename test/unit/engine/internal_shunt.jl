@testitem "Engine / internal shunt / geometry and local operator" tags=[:unit, :parametric] begin
    using LinearAlgebra
    E = LineCableModels.Engine
    include(joinpath(pkgdir(LineCableModels), "test", "support", "internal_shunt.jl"))
    design = internal_shunt_test_design()
    group, bp = internal_shunt_test_group(design)
    @test group.terminals == 1:3
    @test length(group.wires) == 6
    @test length(group.tapes) == 1
    @test group.epsilon ≈ 1
    @test only(group.left).epsilon ≈ 3.0
    @test only(group.right).epsilon ≈ 3.0
    renamed = internal_shunt_test_design(suffix = "_renamed")
    other, _ = internal_shunt_test_group(renamed)
    values = Vector(group)
    @test all(pair -> first(pair) === last(pair), zip(values, Vector(other)))
    @test values[1:3] == [group.a, group.b, group.epsilon]
    # The group rebuilt from its descriptors has the same descriptors.
    rebuilt = @inferred E.ShuntModel.DielectricGroup(group, values)
    @test Vector(rebuilt) == values && rebuilt.reference == group.reference
    host = only(filter(
        placed->placed.primitive isa LineCableModels.DataModel.DifferenceShape &&
                placed.primitive.outer isa Annulus,
        design.geometry.regions))
    incomplete = merge((primitive = host.primitive, source = host.source),
        (primitive = LineCableModels.DataModel.DifferenceShape(host.primitive.outer,
            Base.front(host.primitive.holes)),))
    extraction = E.ShuntModel._shunt_extract(design, bp, Float64, 1, Formulation().formulas)
    @test only(extraction.groups).terminals == group.terminals && isempty(extraction.rejected)
    # A host whose holes differ from its members reports why it is not a group.
    for primitive in (incomplete.primitive, LineCableModels.DataModel.DifferenceShape(
            host.primitive.outer, (
                first(host.primitive.holes), Base.front(host.primitive.holes)...)))
        groups = E.ShuntModel.DielectricGroup{Float64}[]
        rejected = E.ShuntModel.ShuntRejection[]
        E.ShuntModel._shunt_host_group!(groups, rejected, design, bp,
            merge(incomplete, (; primitive)), 1, Float64, Formulation().formulas)
        @test isempty(groups)
        @test occursin("coatings or dielectric islands", only(rejected).reason.message)
    end
    g = group
    level = (wire = 32, order = 16, quadrature = 128, modes = 256)
    selected = E.ShuntModel.Formula(:boundary; options = (resolution = level,))
    result = @inferred E.ShuntModel._shunt_capacitance(selected, g)
    @test result.C ≈ transpose(result.C) rtol=1e-5
    @test minimum(eigvals(Symmetric(result.C))) > 0
    @test sum(result.C[1, :]) > 0
    @test boundary_residual(selected, g, result).boundary < 0.06
    @test result.workspace.plan.rows == result.diagnostic.equations
    # The group's wires form one ring, which the solve separates by Fourier index. The same
    # solver without the ring, a dense solve of the same points, gives the same capacitance.
    @test result.diagnostic.ring == length(g.wires)
    @test result.diagnostic.condition >= 1
    dense = E.ShuntModel._shunt_plan(selected, g,
        merge(E.ShuntModel._shunt_layout(g, level), (; ring = 0)))
    buffers = E.ShuntModel.initialize_buffers(selected, Float64, g, dense, (;))
    unreduced = E.ShuntModel._shunt_capacitance!(selected, (; input = g, plan = dense, buffers))
    @test unreduced.diagnostic.ring == 0
    @test unreduced.C ≈ result.C rtol=1e-12
    doubled = E.ShuntModel.DielectricGroup(g.design, g.assembly, g.terminals, g.reference,
        g.a, g.b, 2g.epsilon, [merge(l, (epsilon = 2l.epsilon,)) for l in g.left],
        [merge(l, (epsilon = 2l.epsilon,)) for l in g.right], g.wires, g.tapes)
    @test E.ShuntModel._shunt_capacitance(selected, doubled).C ≈ 2result.C rtol=1e-8
    rotation = 0.317
    rotated = E.ShuntModel.DielectricGroup(g.design, g.assembly, g.terminals, g.reference,
        g.a, g.b, g.epsilon, g.left, g.right,
        [merge(w, (x = w.x*cos(rotation)-w.y*sin(rotation),
             y = w.x*sin(rotation)+w.y*cos(rotation))) for w in g.wires],
        [merge(t, (phi = t.phi+rotation,)) for t in g.tapes])
    @test E.ShuntModel._shunt_capacitance(selected, rotated).C ≈ result.C rtol=2e-4
    formulation=Formulation(shunt_model = formula(:boundary; options = (resolution = level,)))
    blueprints = [LineCableModels.DataModel.flatten(d, Float64) for d in (design, renamed)]
    pair = only(E.internal_potential_coefficients([design, renamed], blueprints, [formulation]))
    @test pair.details.solves == 1
    @test pair.blocks[1].C === pair.blocks[2].C
    # The per-group loop is inferred from the groups that it solves.
    groups = E.ShuntModel.dielectric_groups(formulation.formulas.shunt_model,
        design, first(blueprints), formulation.formulas, 1).groups
    solutions = @inferred E._internal_capacitances(formulation.formulas.shunt_model, groups)
    @test first(only(solutions)) == pair.blocks[1].C
    @test isempty(only(E.internal_potential_coefficients([design], blueprints[1:1],
        [Formulation(insulation_admittance = :lossy)])).blocks)
    cable = E.LocalCableData(first(blueprints))
    shunt = only(E.internal_potential_coefficients([design], blueprints[1:1], [formulation]))
    p, y = zeros(ComplexF64, 3, 3), zeros(ComplexF64, 3, 3)
    layers = zeros(ComplexF64, length(bp.dielectrics))
    coefficients, tails = zeros(ComplexF64, 3), zeros(ComplexF64, 3)
    local_shunt = shunt.blocks
    admittivity = similar(layers)
    formulas = Formulation().formulas
    rho = [constitutive(formulas.temperature_dependence, material, 20.0)
           for material in cable.dielectric_materials]
    E.dielectric!(admittivity, cable, rho, formulas, 50.0, 20.0)
    E.cable_potential!(p, cable, shunt, admittivity, 100pi*im,
        layers, coefficients, tails)
    E.cable_admittance!(y, cable, shunt, admittivity, 100pi*im, layers)
    @test y ≈ 100pi*im*inv(p) rtol=1e-10
    # Warm the measured entry points directly: the enclosing cable calculation
    # may inline them without compiling a standalone call.
    E._shunt_potential!(p, local_shunt)
    E._shunt_admittance!(y, local_shunt, 100pi*im)
    @test @allocated(E._shunt_potential!(p, local_shunt)) == 0
    @test @allocated(E._shunt_admittance!(y, local_shunt, 100pi*im)) == 0
end

@testitem "Engine / internal shunt / jacket reference and strict coverage" tags=[:unit, :engine] begin
    using LinearAlgebra
    E = LineCableModels.Engine
    DM = LineCableModels.DataModel
    metal = Material(:conductor, 2e-8)
    insulation = Material(:insulator, Inf, 3.0)
    semiconductor = Material(:semicon, 1000.0, 1000.0)
    host = Material(:insulator, Inf, 1.0)
    core = terminal(:core, Region(:core, Disk(0.002), metal))
    wires = Group(:screen, Region(:wire, Disk(0.0003), metal);
        pattern = Ring(4; r = 0.0043))
    enclosure = Enclosure(:host, wires; primitive = Annulus(0.004, 0.0046), fill = host)
    layers = Stack(core,
        Region(:inner_semicon, Shell(0.0005), semiconductor),
        Region(:inner_insulation, Shell(0.0015), insulation), enclosure,
        Region(:outer_semicon, Shell(0.0002), semiconductor),
        Region(:jacket, Shell(0.0008), insulation))
    exposed = build(CableDesign, "jacket-reference", layers)
    shielded = build(CableDesign, "metal-reference", Stack(layers,
        terminal(:shield, Region(:shield, Shell(0.0002), metal)),
        Region(:cover, Shell(0.0007), insulation)))
    selected = CableConstantsFormulation(shunt_model = formula(:boundary;
        options = (resolution = (wire = 16, order = 8, quadrature = 64, modes = 128),)))
    responses = [only(E.internal_potential_coefficients([design], [DM.flatten(design)], [selected]))
                 for design in (exposed, shielded)]
    open_block, closed_block = only(responses[1].blocks), only(responses[2].blocks)
    # A grounded jacket surface and a real closed shield at that same radius
    # impose the same local electrostatic boundary condition.
    @test open_block.C ≈ closed_block.C rtol=1e-12
    @test open_block.P ≈ closed_block.P rtol=1e-12
    @test open_block.terminals == 1:2 && open_block.assembly == 1:2
    @test open_block.reference == E.ShuntJacket(0.0056)
    @test closed_block.reference == E.ShuntShield(0.0056)
    @test convert(E.InternalCapacitance{Float32}, open_block).reference ==
          E.ShuntJacket(Float32(open_block.reference.radius))
    @test responses[1].covered == [true, true]
    @test occursin("ShuntJacket", only(responses[1].details.domains).message)
    @test only(responses[1].details.domains).terminals == 1:2
    @test responses[1].details.requested === :boundary
    @test all(responses[1].covered)
    group = only(E.ShuntModel.dielectric_groups(selected.formulas.shunt_model,
        exposed, DM.flatten(exposed), selected.formulas, 1).groups)
    @test length(group.left) == length(group.right) == 2
    @test [layer.epsilon for layer in group.left] ≈ [1000.0, 3.0]
    @test [layer.epsilon for layer in group.right] ≈ [1000.0, 3.0]
    p, y = zeros(ComplexF64, 2, 2), zeros(ComplexF64, 2, 2)
    E._shunt_potential!(p, responses[1].blocks)
    E._shunt_admittance!(y, responses[1].blocks, 100pi*im)
    @test y ≈ 100pi*im*inv(p) rtol=1e-12
    # Separate cable assemblies keep independent jacket references.
    pair = build(CableDesign, "two-jackets",
        assembly(at(exposed.tree, -0.01, 0.0),
            at(build(CableDesign, "other-jacket", Stack(
                terminal(:other_core, Region(:other_core, Disk(0.002), metal)),
                Region(:other_inner, Shell(0.002), insulation),
                Enclosure(:other_host,
                    Group(:other_screen, Region(:other_wire, Disk(0.0003), metal);
                        pattern = Ring(4; r = 0.0043));
                    primitive = Annulus(0.004, 0.0046), fill = host),
                Region(:other_cover, Shell(0.001), insulation))).tree, 0.01, 0.0)))
    pair_response = only(E.internal_potential_coefficients([pair], [DM.flatten(pair)], [selected]))
    @test all(block -> issubset(block.terminals, block.assembly), pair_response.blocks)
    combined = zeros(ComplexF64, 4, 4)
    E._shunt_admittance!(combined, pair_response.blocks, 100pi*im)
    @test combined[1:2, 1:2] ≈ y rtol=1e-12
    @test iszero(combined[1:2, 3:4]) && iszero(combined[3:4, 1:2])
    for T in (Float32, Float64)
        response = only(E.internal_potential_coefficients([exposed], [DM.flatten(exposed, T)], [selected]))
        @test response.blocks isa Vector{E.InternalCapacitance{T}}
        @test all(isfinite, only(response.blocks).C)
    end
    scalar = @inferred compute(CableConstantsProblem(exposed), selected)
    gridded = compute(CableConstantsProblem(Grid((exposed,)); frequency = Grid((50.0, 60.0))),
        CableConstantsFormulation(shunt_model = Grid((selected.definitions.shunt_model,))))
    @test length(gridded) == 2
    @test all(result -> isapprox(result.C, scalar.C; rtol = 1e-12), gridded)
    @test all(result -> result.C == compute(
        CableConstantsProblem(exposed; frequency = result.frequency), selected).C, gridded)
    # A closed annulus uses the exact radial solution for the selected model.
    radial = build(CableDesign, "radial", Stack(core,
        Region(:radial_insulation, Shell(0.002), insulation)))
    exact = compute(CableConstantsProblem(radial), selected)
    @test exact.C == compute(CableConstantsProblem(radial)).C
    @test details(exact).data.shunt_model.requested === :boundary
    @test details(exact).data.shunt_model.solves == 0
    @test only(details(exact).data.shunt_model.domains).reason === :radial
    # The selected formulation rejects a terminal in a second host.
    unsupported = build(CableDesign, "two-hosts", Stack(layers,
        Enclosure(:second_host,
            Group(:outer_screen, Region(:outer_wire, Disk(0.0003), metal);
                pattern = Ring(4; r = 0.0059));
            primitive = Annulus(0.0056, 0.0062), fill = host),
        Region(:second_jacket, Shell(0.001), insulation)))
    @test_throws BoundarySolveError compute(CableConstantsProblem(unsupported), selected)
    @test_throws r"next conductor outward is not a closed concentric shield" compute(
        CableConstantsProblem(unsupported), selected)
    @test_throws r"lossless dielectric laws" compute(CableConstantsProblem(exposed),
        CableConstantsFormulation(shunt_model = selected.definitions.shunt_model,
            semicon_admittance = :lossy))
end

@testitem "Engine / internal shunt / touching and overlapping wires expose only outer surfaces" tags=[:unit, :engine] begin
    using LinearAlgebra
    E = LineCableModels.Engine
    include(joinpath(pkgdir(LineCableModels), "test", "support", "internal_shunt.jl"))
    SM = E.ShuntModel
    metal = Material(:conductor, 2e-8)
    dielectric = Material(:insulator, Inf, 3.0)
    host = Material(:insulator, Inf, 1.0)
    pitch, radius = sqrt(2)*0.001, 0.001
    design = build(CableDesign, "capacity-armour", Stack(
        terminal(:core, Region(:core, Disk(0.0002), metal)),
        Region(:inner, Shell(pitch-radius-0.0002), dielectric),
        Enclosure(:wire_host,
            Group(:armour, Region(:wire, Disk(radius), metal);
                pattern = Ring(capacity(); r = pitch));
            primitive = Annulus(pitch-radius, pitch+radius), fill = host),
        Region(:jacket, Shell(0.001), dielectric)))
    blueprint = LineCableModels.DataModel.flatten(design)
    selected = CableConstantsFormulation(shunt_model = formula(:boundary;
        options = (resolution = (wire = 32, order = 8, quadrature = 64, modes = 256),)))
    response = only(E.internal_potential_coefficients([design], [blueprint], [selected]))
    @test response.details.requested === :boundary && response.details.solves == 1
    @test only(response.blocks).terminals == 1:2
    @test all(isfinite, only(response.blocks).C)
    @test minimum(eigvals(Symmetric(only(response.blocks).C))) > 0
    # A common scale perturbation preserves both tangency and 2D capacitance.
    group = only(SM.dielectric_groups(selected.formulas.shunt_model, design,
        blueprint, selected.formulas, 1).groups)
    values = Vector(group)
    direction = copy(values)
    direction[[3, 6, 9]] .= 0
    sensitivity = SM.internal_capacitance(SM.Formula(:boundary;
            options = (resolution = (wire = 16, order = 8, quadrature = 64, modes = 256),)),
        group, reshape(direction, :, 1))
    @test norm(sensitivity.tangents) < 1e-6norm(sensitivity.C)
    # Tangency removes isolated contact points. Overlap removes buried arcs.
    for separation in (0.4, 0.3)
        wires = [(x = 2-separation/2, y = 0.0, r = 0.2, terminal = 2),
            (x = 2+separation/2, y = 0.0, r = 0.2, terminal = 2)]
        g = annular_group(1.0, 3.0, 1.0; wires)
        boundary = selected.formulas.shunt_model
        plan = SM._shunt_plan(boundary, g,
            SM._shunt_layout(g, boundary.options.data.resolution))
        @test plan.ring == 0
        @test length(plan.points) < 4plan.level.wire == 2length(plan.sources)
        @test length(plan.points) == length(plan.terminals)
        @test all(point -> minimum(abs(point-complex(wire.x, wire.y))-wire.r
            for wire in wires) >= -SM._shunt_tol(0.2), plan.points)
        result = SM._shunt_capacitance(selected.formulas.shunt_model, g)
        @test all(isfinite, result.C)
        @test minimum(eigvals(Symmetric(result.C))) > 0
        @test boundary_residual(selected.formulas.shunt_model, g, result).wire < 0.06
        independent = [wires[1], merge(wires[2], (terminal = 3,))]
        @test_throws r"touching wire terminals" SM.validate(
            independent, SM.DielectricGroup, (cable_id = "different-terminals",), 3.0)
    end
end

@testitem "Engine / internal shunt / nested shields and reduction" tags=[:unit, :engine] begin
    using LinearAlgebra
    E = LineCableModels.Engine
    # 2 dielectric groups share a closed shield (terminal 3). The outer group's
    # inner-anchor charge includes all terminals enclosed by that shield.
    a = (Diagonal([2.0, 3.0])+[1.0, -1.0]*[1.0 -1.0]) .* 1e-9
    b = (Diagonal([5.0, 7.0])+2*[1.0, -1.0]*[1.0 -1.0]) .* 1e-9
    blocks = [E.InternalCapacitance(1:5, 1:3, E.ShuntShield(0.004), a, inv(a)),
        E.InternalCapacitance(1:5, 3:5, E.ShuntShield(0.008), b, inv(b))]
    # Outermost jacket to the external reference remains a radial contribution.
    p = fill(inv(7e-9), 5, 5)
    y = zeros(ComplexF64, 5, 5)
    s = 100pi*im
    y[5, 5] = s*7e-9
    E._shunt_potential!(p, blocks)
    E._shunt_admittance!(y, blocks, s)
    @test y ≈ s*inv(p) rtol=1e-12
    # Bundling terminals sums their charges and imposes equal potentials.
    incidence = [1 0 0; 0 1 0; 0 1 0; 0 0 1; 0 0 1]
    bundled = transpose(incidence)*y*incidence
    @test bundled ≈ transpose(bundled)
    @test minimum(eigvals(Symmetric(imag.(bundled)))) > 0
end

@testitem "Engine / internal shunt / multiple independent open terminals" tags=[:unit, :parametric] begin
    using LinearAlgebra
    E = LineCableModels.Engine
    include(joinpath(pkgdir(LineCableModels), "test", "support", "internal_shunt.jl"))
    original, _ = internal_shunt_test_group(internal_shunt_test_design(tapes = false))
    # Independent test of the numerical terminal map, not permission to bypass
    # the existing public coaxial lowering rules for radial conductor ordering.
    wires = [merge(wire, (terminal = iseven(index) ? 2 : 3,))
             for (index, wire) in pairs(original.wires)]
    group = E.ShuntModel.DielectricGroup(
        original.design, 1:4, 1:4, original.reference, original.a, original.b,
        original.epsilon, original.left, original.right, wires, original.tapes)
    @test group.terminals == 1:4
    @test sort(unique(getproperty.(group.wires, :terminal))) == [2, 3]
    result = E.ShuntModel._shunt_capacitance(E.ShuntModel.Formula(:boundary;
        options = (resolution = (wire = 32, order = 16, quadrature = 128, modes = 256),)), group)
    @test size(result.C) == (3, 3)
    @test result.C ≈ transpose(result.C) rtol=1e-6
    @test result.C[2, 3] < 0
    @test minimum(eigvals(Symmetric(result.C))) > 0
end
