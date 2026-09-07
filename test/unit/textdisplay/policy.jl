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
          "CableConstants(assemblies=1, frequency=50.0)"
    @test sprint(show, MIME"text/plain"(), earth) == join(
        (
            "EarthModel · homogeneous",
            "├─ air        ρ=∞  εᵣ=1  μᵣ=1",
            "└─ earth      ρ=100 Ω·m  εᵣ=10  μᵣ=1"
        ),
        '\n')
    @test sprint(show, MIME"text/plain"(), parameters) == join(
        (
            "LineParameters · phase domain",
            "├─ f  2 points · 1 Hz … 1 MHz",
            "├─ Z  2×2×2 · Ω/m",
            "└─ Y  2×2×2 · S/m"
        ),
        '\n')

    inactive=Material(:insulator, Inf, 2.3, 1.0, 20.0, 0.0)
    inactive_text=sprint(show, MIME"text/plain"(), inactive)
    @test !occursin("α", inactive_text)
    @test !occursin("tanδ", inactive_text)
    @test !endswith(inactive_text, '\n')
end

@testitem "TextDisplay / scientific values preserve signs and physical units" tags=[:unit] begin
    const TD = LineCableModels.TextDisplay
    const U = LineCableModels.Units

    @test TD.value(0.0) == "0"
    @test TD.value(-0.0) == "0"
    @test TD.value(NaN) == "NaN"
    @test TD.value(1.23456789; sigdigits=4) == "1.235"
    @test TD.value(1.25 + 2.5im) == "1.25 + 2.5im"
    @test TD.value(1.25 - 2.5im) == "1.25 − 2.5im"
    @test TD.value(complex(1.0, -0.0)) == "1 − 0im"
    @test TD.value(complex(Inf, -Inf)) == "∞ − ∞im"
    @test TD.value(:unrecorded) == ":unrecorded"
    @test_throws ArgumentError TD.value(1.0; sigdigits=0)
    @test_throws ArgumentError TD.value(1.0 + im; sigdigits=0)

    for (number, unit, expected) in (
            (0.0, :meter, "0 m"),
            (-2e-3, :meter, "-2 mm"),
            (0.1, :hertz, "100 mHz"),
            (1e6, :hertz, "1 MHz"),
            (20.0, :celsius, "20 °C"),
            (0.00393, :kelvin_inverse, "0.00393 K⁻¹"),
            (30.0, :degree, "30°"),
            (0.001, :dimensionless, "0.001"),
            (Inf, :ohm_meter, "∞ Ω·m"),
            (-Inf, :celsius, "−∞ °C"),
            (NaN, :dimensionless, "NaN"),
            (Inf, :degree, "∞°"),
        )
        @test TD.engineering(number, unit) == expected
    end
    @test TD.angle(0.0) == "0°"
    @test TD.angle(pi) == "π"
    @test TD.angle(-pi) == "−π"
    @test TD.angle(-2pi) == "-2π"
    @test TD.angle(pi / 6) == "30°"

    # These independent unit expectations prevent a report from changing only
    # its caption while retaining values on the wrong per-length scale.
    for (observed, selector, basis, expected) in (
            (1e-4, R, :pul, "0.1 Ω/km"),
            (0.1, R, :total, "0.1 Ω"),
            (2e-7, L, :pul, "0.2 mH/km"),
            (2e-4, L, :total, "0.2 mH"),
            (3e-10, C, :pul, "0.3 μF/km"),
            (3e-7, C, :total, "0.3 μF"),
            (4e-9, G, :pul, "4.0e-6 S/km"),
            (0.001 + 0.002im, Z, :pul, "1 + 2im Ω/km"),
            (0.003 - 0.004im, Y, :total, "0.003 − 0.004im S"),
        )
        @test TD.quantity(observed, U.quantity(selector), basis) == expected
    end
    @test TD.quantity(pi / 2, U.quantity(Z, angle)) == "90°"
    @test TD.quantity(0.25, U.Quantity{:dimensionless}()) == "0.25"
    @test_throws ArgumentError TD.quantity(1.0, U.quantity(R), :unknown)

    for (object, compact, summary_text) in (
            (U.Unit(:ohm, :milli), "mΩ", "Physical unit mΩ"),
            (U.units(:micro, :farad; per=(:kilo, :meter)), "μF/km", "Unit expression μF/km"),
            (U.UnitExpr(), "", "Unit expression "),
            (U.quantity(R), "R · Series resistance", "Series resistance"),
            (U.Quantity{:dimensionless}(), "Dimensionless", "Dimensionless"),
        )
        @test sprint(show, object) == compact
        @test sprint(show, MIME"text/plain"(), object) == compact
        @test sprint(summary, object) == summary_text
    end
    @test TD.name(typeof(formula(:default))) == "FormulaDefinition"
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
    insulating_layers=ntuple(
        index->Region(Symbol(:layer_, index), Annulus(index*1.0e-3, (index+1)*1.0e-3), insulator),
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
    earth=build(EP.EarthModel, (
        EP.EarthLayer(100.0, 10.0, 1.0, 2.0),
        EP.EarthLayer(30.0, 15.0, 1.0, 5.0),
        EP.EarthLayer(500.0, 8.0, 1.0)
    ))
    grid=Grid((1.0, 2.0, 3.0, 4.0, 5.0))
    build_calls=Ref(0)
    space=Gridspace{CableConstants}(
        (r, l, c)->begin
            build_calls[]+=1
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
        monte_carlo
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
        :monte_carlo_result
    )

    function structural_snapshot(names, objects, display_size)
        sections=map(names, objects) do name, object
            rendered=sprint(
                show,
                MIME"text/plain"(),
                object;
                context = IOContext(
                    IOBuffer(), :compact=>false, :limit=>true,
                    :displaysize=>display_size
                )
            )
            return string("## ", name, '\n', rendered)
        end
        return join(sections, "\n\n")*'\n'
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
            IOBuffer(), :compact=>false, :limit=>true,
            :displaysize=>(40, 120)
        )
        narrow_context=IOContext(
            IOBuffer(), :compact=>false, :limit=>true,
            :displaysize=>(6, 48)
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
    @test occursin("⋮",
        sprint(
            show,
            MIME"text/plain"(),
            stack;
            context = IOContext(
                IOBuffer(), :limit => true, :displaysize => (5, 48)
            )
        ))
    @test !occursin("0.0 + 0.0im", sprint(show, MIME"text/plain"(), parameters))
end
