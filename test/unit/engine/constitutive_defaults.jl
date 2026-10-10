@testitem "Engine / explicit lossless defaults and lossy radial admittance" tags=[:unit, :parametric] begin
    using Test
    using LineCableModels
    import EzXML
    const E = LineCableModels.Engine

    @test formula(:default) isa LineCableModels.FormulaDefinition
    ε0 = 8.8541878128e-12
    for (owner, kind) in ((E.InsulationAdmittance, :insulator),
        (E.SemiconAdmittance, :semicon))
        material = Material(kind, 1.0, 100.0; tan_delta = 0.02)
        lossless = owner.Formula(:default)
        lossy = owner.Formula(:lossy)
        @test formula_id(lossless) === :lossless
        for frequency in (0.1, 50.0, 1.0e6)
            ω = 2π * frequency
            actual = @inferred constitutive(lossless, material, frequency, 20.0)
            @test real(actual) == 0
            @test imag(actual) ≈ ω * ε0 * material.eps_r
            actual = @inferred constitutive(lossy, material, frequency, 20.0)
            εstar = ε0 * material.eps_r * (1 - im * material.tan_delta) +
                    1 / (im * ω * material.rho)
            @test actual ≈ im * ω * εstar
        end
    end

    # Ametani (2004), Fig. 4 geometry and Eqs. (14)-(15), not a regenerated baseline.
    b, c, r0 = 30.45e-3, 35.45e-3, 71.15e-3
    screen = Material(:semicon, 1.0, 1000.0)
    insulation_material = Material(:insulator, Inf, 3.1)
    design = build(CableDesign,
        "ametani-radial",
        terminal(:core,
            solid(Material(:conductor, 1.82e-8), Disk(b)),
            LineCableModels.screen(screen; t = c - b),
            insulation(insulation_material; t = r0 - c)))
    selected = E.CableConstantsFormulation(
        insulation_admittance = :lossy, semicon_admittance = :lossy)
    for frequency in (50.0, 60.0)
        ω = 2π * frequency
        εs = ε0 * screen.eps_r + 1 / (im * ω * screen.rho)
        ys = im * ω * 2π * εs / log(c / b)
        yi = im * ω * 2π * ε0 * insulation_material.eps_r / log(r0 / c)
        expected = inv(inv(ys) + inv(yi))
        result = compute(E.CableConstantsProblem(design; frequency), selected)
        @test only(result.G) ≈ real(expected)
        @test only(result.C) ≈ imag(expected) / ω
        lossless = compute(E.CableConstantsProblem(design; frequency))
        @test only(lossless.G) == 0
    end

    # The exporters evaluate the group's admittivity at the export frequency, where it
    # matches the original layers, not at every frequency.
    polar = Material(:insulator, 2.0e9, 3.1; tan_delta = 0.04)
    screen = Material(:semicon, 1.0e5, 1000.0; tan_delta = 0.02)
    source = build(CableDesign,
        "polarization-equivalence",
        terminal(:core,
            solid(Material(:conductor, 1.82e-8), Disk(b)),
            LineCableModels.screen(screen; t = c - b),
            insulation(polar; t = r0 - c)))
    system = build(LineCableSystem, [source], [(0.0, -1.0)]; connections = [Dict(:core => 1)])
    for frequency in (50.0, 60.0)
        equivalent = homogenize(source)
        actual = compute(E.CableConstantsProblem(source; frequency), selected)
        reconstructed = compute(E.CableConstantsProblem(equivalent; frequency), selected)
        @test reconstructed.C ≈ actual.C
        @test reconstructed.G ≈ actual.G
        # The ATP export writes the group's capacitance through epsI and its conductance.
        epsI, Gext = mktempdir() do directory
            path = export_data(:atp, system, EarthModel(100.0, 10.0, 1.0);
                formulation = selected, base_freq = frequency,
                file_name = joinpath(directory, "polarization.xml"))
            core = first(EzXML.findall("//conductor", EzXML.readxml(path)))
            parse(Float64, core["epsI"]), parse(Float64, core["Gext"])
        end
        @test 2π * ε0 * epsI / log(r0 / b) ≈ only(actual.C)
        @test Gext ≈ only(actual.G)
    end
end
