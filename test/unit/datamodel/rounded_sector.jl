@testitem "DataModel / RoundedSector / exact geometry and conformal shells" tags=[:unit] begin
    const DM=LineCableModels.DataModel

    primitive=RoundedSector(
        span = deg2rad(119.0),
        r_base = 1.10e-3,
        r_back = 10.24e-3,
        fillet = 1.02e-3
    )
    @test fieldnames(typeof(primitive)) == (:span, :r_base, :r_back, :fillet)
    for forbidden in (
        :material, :terminal, :n_sectors, :d_insulation, :at,
        :vertices, :centroid, :area, :resistance, :gmr
    )
        @test forbidden ∉ fieldnames(typeof(primitive))
    end

    shape=DM.resolve(DM.EmptyBoundary(), primitive)
    @test shape isa DM.RoundedSectorShape
    @test DM.area(shape) ≈ 9.207305021593468e-5 rtol=4eps(Float64)
    @test DM.perimeter(shape) ≈ 3.780313381817682e-2 rtol=4eps(Float64)
    @test DM.centroid(shape)[1] ≈ 6.131368623545798e-3 rtol=4eps(Float64)
    @test abs(DM.centroid(shape)[2]) <= 16eps(Float64)

    exact=(DM.area(shape), DM.perimeter(shape), DM.centroid(shape))
    coarse=DM.tessellate(shape; points_per_arc = 3)
    dense=DM.tessellate(shape; points_per_arc = 257)
    @test exact == (DM.area(shape), DM.perimeter(shape), DM.centroid(shape))
    @test length(dense) > length(coarse)
    @test all(point -> point isa NTuple{2, Float64}, coarse)
    @test all(point -> point isa NTuple{2, Float64}, dense)
    @test !occursin("GeometryBasics", string(typeof(coarse)))
    @test !occursin("Makie", string(typeof(coarse)))

    for angle in range(-pi, pi; length = 37)
        sampled=maximum(point -> point[1] * cos(angle) + point[2] * sin(angle), dense)
        @test DM.support(shape, angle) + 2e-7 >= sampled
    end
    coarse_angle=0.31
    coarse_support=maximum(
        point -> point[1] * cos(coarse_angle) + point[2] * sin(coarse_angle),
        coarse
    )
    @test DM.support(shape, coarse_angle) - coarse_support > 4e-4

    pose=Pose2(0.2, -0.1, pi / 3)
    placed=DM.resolve(pose, shape)
    local_centre=DM.centroid(shape)
    expected=(
        pose.x + cos(pose.φ) * local_centre[1] - sin(pose.φ) * local_centre[2],
        pose.y + sin(pose.φ) * local_centre[1] + cos(pose.φ) * local_centre[2]
    )
    @test collect(DM.centroid(placed)) ≈ collect(expected)
    @test primitive == shape.primitive == placed.primitive

    thickness=0.5e-3
    shell=DM.resolve(shape, Shell(thickness))
    @test shell isa DM.ShellShape
    @test shell.outer.primitive.span == primitive.span
    @test shell.outer.primitive.r_base == primitive.r_base - thickness
    @test shell.outer.primitive.r_back == primitive.r_back + thickness
    @test shell.outer.primitive.fillet == primitive.fillet + thickness
    @test DM.area(shell) == DM.area(shell.outer) - DM.area(shell.inner)
    @test DM.perimeter(shell) ==
          DM.perimeter(shell.inner) + DM.perimeter(shell.outer)
    @test DM.area(shell.outer) ≈
          DM.area(shape) + DM.perimeter(shape) * thickness + pi * thickness^2
    @test DM.perimeter(shell.outer) ≈ DM.perimeter(shape) + 2pi * thickness

    inner_area=DM.area(shell.inner)
    outer_area=DM.area(shell.outer)
    inner_centre=DM.centroid(shell.inner)
    outer_centre=DM.centroid(shell.outer)
    expected_shell_centre=(
        (outer_area * outer_centre[1] - inner_area * inner_centre[1]) /
        (outer_area - inner_area),
        (outer_area * outer_centre[2] - inner_area * inner_centre[2]) /
        (outer_area - inner_area)
    )
    @test collect(DM.centroid(shell)) ≈ collect(expected_shell_centre)

    second=DM.resolve(DM.boundary(shell), Shell(0.25e-3))
    @test second.inner == shell.outer
    @test second.outer.primitive.r_back == primitive.r_back + 0.75e-3
    @test second.outer.primitive.fillet == primitive.fillet + 0.75e-3
    @test_throws DomainError DM.resolve(shape, Shell(primitive.r_base + eps()))
    @test_throws DomainError RoundedSector(0.2, 0.0, 1.0, 0.95)
end

@testitem "DataModel / RoundedSector / terminal assembly and Gridspace" tags=[:unit] begin
    const DM=LineCableModels.DataModel

    copper=Material(kind = :conductor, rho = 1.7241e-8, mu_r = 0.999994,
        alpha = 0.00393)
    insulation_material=Material(kind = :insulator, rho = 1.97e14, eps_r = 2.5)
    primitive=RoundedSector(
        span = deg2rad(119.0),
        r_base = 1.10e-3,
        r_back = 10.24e-3,
        fillet = 1.02e-3
    )
    phase=terminal(
        :phase,
        Region(:core, primitive, copper),
        Region(:insulation, Shell(0.5e-3), insulation_material)
    )
    root=assembly(
        phase;
        pattern = Ring(3; r = 0),
        names = (:a, :b, :c)
    )
    design=build(CableDesign, "sectorized", root)
    @test design.terminal_order == [:a, :b, :c]
    @test count(!isnothing, getproperty.(design.geometry.regions, :terminal)) == 3
    @test count(region -> region.primitive isa DM.RoundedSectorShape,
        design.geometry.regions) == 3
    @test count(region -> region.primitive isa DM.ShellShape,
        design.geometry.regions) == 3
    @test !(design.geometry.outer isa DM.Disk)
    @test_throws MethodError DM.resolve(design.geometry.outer, Shell(0.5e-3))

    first_region=first(design.geometry.regions)
    @test first_region.source.primitive == primitive
    preview_shape=only(DM.preview_shapes(
        first_region, (
            include_label = false,
            label = nothing,
            group = :rounded_sector
        )))
    @test length(preview_shape.geometry.exterior) > 6
    @test all(point -> all(isfinite, point), preview_shape.geometry.exterior)

    conductor_shapes=[region.primitive
                      for region in design.geometry.regions
                      if region.terminal !== nothing]
    @test getproperty.(getproperty.(conductor_shapes, :at), :φ) ≈
          [0.0, 2pi / 3, 4pi / 3]

    varied=RoundedSector(
        span = Grid((deg2rad(118.0), deg2rad(119.0))),
        r_base = Grid((1.00e-3, 1.10e-3)),
        r_back = Grid((10.24e-3, 10.34e-3)),
        fillet = Grid((0.9e-3, 1.02e-3))
    )
    varied_shell=Shell(Grid((0.4e-3, 0.5e-3)))
    varied_phase=terminal(
        :phase,
        Region(:core, varied, copper),
        Region(:insulation, varied_shell, insulation_material)
    )
    designs=build(CableDesign, "sector-grid", varied_phase)
    @test designs isa Gridspace{CableDesign}
    @test eltype(designs) === CableDesign
    @test length(designs) == 32
    @test all(design -> design isa CableDesign, designs)
    @test all(design -> design.root isa DM.Group, designs)
end

@testitem "DataModel / RoundedSector / member-local equivalent-area flattening" tags=[:unit] begin
    const DM=LineCableModels.DataModel

    copper=Material(kind = :conductor, rho = 1.7241e-8, mu_r = 0.999994,
        alpha = 0.00393)
    xlpe=Material(kind = :insulator, rho = 1.97e14, eps_r = 2.5)
    primitive=RoundedSector(
        span = deg2rad(119.0),
        r_base = 1.10e-3,
        r_back = 10.24e-3,
        fillet = 1.02e-3
    )
    phase=terminal(
        :phase,
        Region(:core, primitive, copper),
        Region(:insulation, Shell(0.5e-3), xlpe)
    )
    design=build(CableDesign, "sector-input", assembly(
        phase;
        pattern = Ring(3; r = 0),
        names = (:a, :b, :c)
    ))

    components=DM.flatten(design, 50.0)
    @test getproperty.(components, :name) == [:a, :b, :c]
    source_shape=first(design.geometry.regions).primitive
    shell_shape=design.geometry.regions[2].primitive
    source_area=DM.area(source_shape)
    equivalent_radius=sqrt(source_area / pi)
    outer_radius=sqrt(DM.area(shell_shape.outer) / pi)
    expected_position=DM.centroid(source_shape)
    component=first(components)
    @test component.conductor.cross_section ≈ source_area
    @test component.conductor.r_in == 0
    @test component.conductor.r_ex ≈ equivalent_radius
    @test component.conductor.resistance ≈ copper.rho / source_area
    @test component.conductor.gmr ≈
          equivalent_radius * exp(-copper.mu_r / 4)
    @test collect(component.conductor.position) ≈ collect(expected_position)
    @test component.dielectric.r_in ≈ equivalent_radius
    @test component.dielectric.r_ex ≈ outer_radius
    @test component.dielectric.cross_section ≈ DM.area(shell_shape)
    @test component.dielectric.shunt_capacitance ≈
          DM.shunt_capacitance(equivalent_radius, outer_radius, xlpe.eps_r)
    @test component.dielectric.shunt_conductance ≈
          DM.shunt_conductance(equivalent_radius, outer_radius, xlpe.rho)
    @test DM.tubular_resistance(
        component.conductor.r_in,
        component.conductor.r_ex,
        component.conductor.material.rho
    ) ≈ component.conductor.resistance
    @test DM.tubular_gmr(
        component.conductor.r_ex,
        component.conductor.r_in,
        component.conductor.material.mu_r
    ) ≈ component.conductor.gmr
    @test DM.shunt_capacitance(
        component.dielectric.r_in,
        component.dielectric.r_ex,
        component.dielectric.material.eps_r
    ) ≈ component.dielectric.shunt_capacitance
    @test DM.shunt_conductance(
        component.dielectric.r_in,
        component.dielectric.r_ex,
        component.dielectric.material.rho
    ) ≈ component.dielectric.shunt_conductance
    @test component.dielectric.material.mu_r ≈ xlpe.mu_r

    magnetic_inner=Material(kind = :insulator, rho = 1e12, eps_r = 2.0, mu_r = 1.5)
    magnetic_outer=Material(kind = :insulator, rho = 2e12, eps_r = 3.0, mu_r = 2.5)
    magnetic_layers=[
        (r_in = 0.01, r_ex = 0.015, material = magnetic_inner),
        (r_in = 0.015, r_ex = 0.02, material = magnetic_outer)
    ]
    radial_mu=(
        magnetic_inner.mu_r*log(0.015/0.01) +
        magnetic_outer.mu_r*log(0.02/0.015)
    )/log(0.02/0.01)
    @test DM.equivalent_dielectric_permeability(
        magnetic_layers,
        4.0,
        0.01,
        0.02
    ) ≈ radial_mu*DM.solenoid_factor(4.0, 0.01, 0.02)

    problem=LineParametersProblem(
        design,
        Pose2(0, -1);
        connections = (a = 1, b = 2, c = 3),
        earth_props = Earth(rho = 100),
        frequencies = [50.0]
    )
    flattened=homogenize(design)
    @test flattened.root isa DM.Assembly
    @test flattened.terminal_order == design.terminal_order
    flattened_conductors=filter(
        region->region.source.material.kind===:conductor,
        flattened.geometry.regions
    )
    @test all(zip(
        DM.centroid.(getproperty.(flattened_conductors, :primitive)),
        getproperty.(getproperty.(components, :conductor), :position)
    )) do (actual, expected)
        isapprox(actual[1], expected[1]) && isapprox(actual[2], expected[2])
    end
    @test_throws ArgumentError compute(problem)

    encoded=LineCableModels.ImportExport.serialize_value(design)
    @test encoded["root"]["item"]["item"]["items"][1]["primitive"]["kind"] ==
          "rounded_sector"
    serialized=sprint(show, encoded)
    for derived in (
        "contacts", "geometry", "equivalent_radius", "resistance", "gmr",
        "shunt_capacitance", "shunt_conductance"
    )
        @test !occursin(derived, serialized)
    end
    restored=LineCableModels.ImportExport.deserialize_value(encoded)
    @test restored == design
end

@testitem "DataModel / RoundedSector / legacy equivalence equations" tags=[:unit] begin
    const DM=LineCableModels.DataModel

    # The sector fork used these material values and dimensions. Its scientific
    # contract was the equivalent-area calculation; sampled polygon vertices
    # were only its approximation of the same physical boundary.
    aluminum=Material(
        kind = :conductor,
        rho = 2.8264e-8,
        eps_r = 1.0,
        mu_r = 1.000022,
        T0 = 20.0,
        alpha = 0.00429
    )
    pvc=Material(
        kind = :insulator,
        rho = Inf,
        eps_r = 8.0,
        mu_r = 1.0,
        T0 = 20.0,
        alpha = 0.1
    )
    phase=terminal(
        :phase,
        Region(:core,
            RoundedSector(
                span = deg2rad(119.0),
                r_base = 1.10e-3,
                r_back = 10.24e-3,
                fillet = 1.02e-3
            ),
            aluminum),
        Region(:insulation, Shell(1.1e-3), pvc)
    )
    design=build(CableDesign, "legacy-sector-equivalence", phase)
    component=only(DM.flatten(design, 20.0))

    # Fixed values prevent the compatibility test from merely recomputing its
    # expectations through the same implementation under test.
    @test component.conductor.cross_section ≈ 9.207305021593469e-5 rtol=8eps()
    @test component.conductor.resistance ≈ 3.069736468349179e-4 rtol=8eps()
    @test component.conductor.gmr ≈ 4.216142877891434e-3 rtol=8eps()
    @test component.dielectric.r_in ≈ 5.413664390671869e-3 rtol=8eps()
    @test component.dielectric.r_ex ≈ 6.614694587068139e-3 rtol=8eps()
    @test component.dielectric.cross_section ≈ 4.538477431083816e-5 rtol=8eps()
    @test component.dielectric.shunt_capacitance ≈ 2.221219434950488e-9 rtol=8eps()
    @test component.dielectric.shunt_conductance == 0
end
