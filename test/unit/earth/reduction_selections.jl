@testitem "Analytical selections / explicit earth reduction preserves formula and indexed geometry" tags=[:unit] begin
    const EP = LineCableModels.Earth
    copper = Material(:conductor, 1.72e-8, 1.0)
    design = build(CableDesign, "explicit-reduction", terminal(:core,
        solid(copper, Disk(0.004)), insulation(Material(:insulator, 1e14, 2.3); t = 0.002)))
    earth = build(EP.EarthModel, (EP.EarthLayer(100.0, 10.0, 1.0, 0.5),
        EP.EarthLayer(500.0, 20.0, 1.0)))
    for (height, identifiers) in ((2.0, (:default, :Gary1976, :Carson1926)),
            (-1.0, (:default, :Pollaczek1926, :Saad1996, :WedepohlWilcox1973)))
        system = build(LineCableSystem, [design, design], [Pose2(0, height), Pose2(0.75, height)];
            connections = [Dict(:core => 1), Dict(:core => 2)])
        original = LineParametersProblem(system; earth_props = earth, frequencies = [50.0])
        equivalent = LineParametersProblem(system;
            earth_props = homogeneous(rho = 500.0, eps_r = 20.0), frequencies = [50.0])
        admittance_identifiers = height > 0 ? (:default,) : (:default, :Pollaczek1926)
        for identifier in identifiers, admittance_id in admittance_identifiers
            requested = Formulation(earth_impedance = formula(identifier;
                equivalent_earth = formula(:default)),
                earth_admittance = formula(admittance_id; equivalent_earth = formula(:default)))
            reduced = compute(original, requested)
            direct = compute(equivalent, Formulation(earth_impedance = identifier,
                earth_admittance = admittance_id))
            @test Z(reduced) ≈ Z(direct) rtol = 1e-10
            @test Y(reduced) ≈ Y(direct) rtol = 1e-10
            @test details(reduced).formulations.effective.earth_impedance === identifier
            @test details(reduced).formulations.effective.earth_admittance === admittance_id
            @test length(original.earth_props.layers) == 3
            @test original.earth_props === earth
        end
    end
end
