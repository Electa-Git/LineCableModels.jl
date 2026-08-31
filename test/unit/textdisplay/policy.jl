@testitem "TextDisplay / engineering values and target summaries" tags=[:unit] begin
    const TD=LineCableModels.TextDisplay
    const EP=LineCableModels.EarthProps

    copper=Material(:conductor, 1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
    constants=CableConstants(1.33803e-5, 4.04835e-6, 1.83687e-10)
    earth=EP.EarthModel(100.0, 10.0, 1.0)
    frequencies=[1.0, 1.0e6]
    parameters=LineParameters(
        zeros(ComplexF64, 2, 2, 2),
        zeros(ComplexF64, 2, 2, 2),
        frequencies
    )

    @test TD.engineering(0.00183245, :meter) == "1.83245 mm"
    @test TD.engineering(1.7241e-8, :ohm_meter) == "17.241 nΩ·m"
    @test TD.value(Inf) == "∞"
    @test TD.value(-Inf) == "−∞"
    @test TD.angle(2π) == "2π"

    @test sprint(show, copper) ==
          "Material(:conductor; ρ=17.241 nΩ·m, εᵣ=1, μᵣ=0.999994)"
    @test sprint(show, constants) ==
          "CableConstants(assemblies=1, frequency=0.001)"
    @test sprint(show, MIME"text/plain"(), earth) == join((
        "EarthModel · homogeneous",
        "├─ air        ρ=∞  εᵣ=1  μᵣ=1",
        "└─ earth      ρ=100 Ω·m  εᵣ=10  μᵣ=1",
    ), '\n')
    @test sprint(show, MIME"text/plain"(), parameters) == join((
        "LineParameters · phase domain",
        "├─ f  2 points · 1 Hz … 1 MHz",
        "├─ Z  2×2×2 · Ω/m",
        "└─ Y  2×2×2 · S/m",
    ), '\n')

    inactive=Material(:insulator, Inf, 2.3, 1.0, 20.0, 0.0)
    inactive_text=sprint(show, MIME"text/plain"(), inactive)
    @test !occursin("α", inactive_text)
    @test !occursin("tanδ", inactive_text)
    @test !endswith(inactive_text, '\n')
end

@testitem "TextDisplay / bounded structural families" tags=[:unit] setup=[
    TestFixtures
] begin
    const EP=LineCableModels.EarthProps

    conductor=Material(:conductor, 1.7241e-8)
    insulator=Material(:insulator, Inf, 2.3)
    wire=Region(:wire, Disk(1.0e-3), conductor)
    group=Group(
        :core,
        wire;
        pattern = Ring(6; r = 2.0e-3),
        path = Helix(LayRatio(11.0))
    )
    insulating_layers=ntuple(index ->
        Region(Symbol(:layer_, index), Annulus(index * 1.0e-3, (index + 1) * 1.0e-3), insulator),
        10)
    stack=Stack(group, insulating_layers...)
    repeated=Assembly(
        group;
        pattern = Ring(2; r = 6.0e-3),
        names = (:left, :right)
    )
    enclosure=Enclosure(
        :duct,
        group;
        primitive = Disk(10.0e-3),
        fill = insulator
    )
    design=TestFixtures.mv_cable_design()
    system=TestFixtures.three_phase_system()
    earth=EP.EarthModel(100.0, 10.0, 1.0; thickness = 2.0)
    add!(earth, EP.EarthLayer(30.0, 15.0, 1.0, 5.0))
    add!(earth, EP.EarthLayer(500.0, 8.0, 1.0))
    grid=Grid((1.0, 2.0, 3.0, 4.0, 5.0))
    build_calls=Ref(0)
    space=Gridspace{CableConstants}(
        (r, l, c) -> begin
            build_calls[] += 1
            CableConstants(r, l, c)
        end,
        (Grid((1.0, 2.0, 3.0)), Grid((2.0,)), Grid((3.0, 4.0)))
    )
    materials=MaterialsLibrary()
    cables=CablesLibrary()
    add!(cables, design)
    parameters=TestFixtures.two_conductor_results()
    parametric=ParametricResult(Combinatorial(Formulation()), [parameters])
    monte_carlo=TestFixtures.cable_monte_carlo_result()

    families=(
        group,
        stack,
        repeated,
        enclosure,
        design.geometry,
        design,
        system,
        earth,
        grid,
        space,
        materials,
        cables,
        parametric,
        parameters,
        monte_carlo,
    )
    family_names=(
        :group,
        :stack,
        :assembly,
        :enclosure,
        :geometry,
        :design,
        :system,
        :earth,
        :grid,
        :gridspace,
        :materials_library,
        :cables_library,
        :parametric_result,
        :line_parameters,
        :monte_carlo_result,
    )

    function structural_snapshot(names, objects, display_size)
        sections = map(names, objects) do name, object
            rendered = sprint(
                show,
                MIME"text/plain"(),
                object;
                context = IOContext(
                    IOBuffer(), :compact => false, :limit => true,
                    :displaysize => display_size
                )
            )
            return string("## ", name, '\n', rendered)
        end
        return join(sections, "\n\n") * '\n'
    end

    snapshot_root=joinpath(pkgdir(LineCableModels), "test", "fixtures", "textdisplay")
    @test structural_snapshot(family_names, families, (40, 120)) ==
          read(joinpath(snapshot_root, "wide.txt"), String)
    @test structural_snapshot(family_names, families, (6, 48)) ==
          read(joinpath(snapshot_root, "narrow.txt"), String)

    for object in families
        compact=sprint(show, object)
        @test !occursin('\n', compact)
        @test !occursin(r"\{[^}]*\}", compact)
        @test !endswith(compact, '\n')

        wide_context=IOContext(
            IOBuffer(), :compact => false, :limit => true,
            :displaysize => (40, 120)
        )
        narrow_context=IOContext(
            IOBuffer(), :compact => false, :limit => true,
            :displaysize => (6, 48)
        )
        wide=sprint(show, MIME"text/plain"(), object; context = wide_context)
        repeated_wide=sprint(show, MIME"text/plain"(), object; context = wide_context)
        narrow=sprint(show, MIME"text/plain"(), object; context = narrow_context)

        @test wide == repeated_wide
        @test !endswith(wide, '\n')
        @test !endswith(narrow, '\n')
        @test length(split(narrow, '\n')) <= 6
        @test all(textwidth(line) <= 48 for line in split(narrow, '\n'))
    end

    @test build_calls[] == 0
    @test sprint(show, space) == "CableConstants parameter space · 6 points"
    @test occursin("⋮", sprint(
        show,
        MIME"text/plain"(),
        stack;
        context = IOContext(
            IOBuffer(), :limit => true, :displaysize => (5, 48)
        )
    ))
    @test !occursin("0.0 + 0.0im", sprint(show, MIME"text/plain"(), parameters))
end
