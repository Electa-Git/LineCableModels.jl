@testitem "ImportExport / PSCAD / selected constitutive laws and operating temperature" tags=[:integration] begin
    using EzXML
    const IE = LineCableModels.ImportExport
    copper = Material(:conductor, 1.72e-8, 1, 1, 20, 0.004)
    semicon = Material(:semicon, 1e4, 40, 1, 20, -0.005; tan_delta=0.02)
    dielectric = Material(:insulator, 1e8, 2.3, 1, 20, -0.003; tan_delta=0.03)
    design = build(CableDesign, "selected-dielectrics", terminal(:core,
        solid(copper, Disk(0.004)), screen(semicon; t=0.0005),
        insulation(dielectric; t=0.002)))
    system = build(LineCableSystem, [design], [Pose2(0, -1)];
        connections=[Dict(:core=>1)])
    earth = homogeneous(rho=100.0)
    frequency = 60.0
    omega = 2pi * frequency
    epsilon0 = 8.8541878128e-12
    admittivity(material, id, temperature=60.0) = im * omega * epsilon0 * material.eps_r +
        (id === :default ? 0 : inv(material.rho * (1 + material.alpha * (temperature - material.T0))) +
         omega * epsilon0 * material.eps_r * material.tan_delta)

    mktempdir() do directory
        for constructor in (Formulation, CableConstantsFormulation),
            insulation_id in (:default, :Ametani2004), semicon_id in (:default, :Ametani2004)
            selected = constructor(insulation_admittance=insulation_id,
                semicon_admittance=semicon_id)
            component = only(IE._pscad_components(design, frequency, selected, 60.0))
            @test component.conductor.material.rho ≈ copper.rho * 1.16
            expected = inv(log(0.0045 / 0.004) / (2pi * admittivity(semicon, semicon_id)) +
                log(0.0065 / 0.0045) / (2pi * admittivity(dielectric, insulation_id)))
            @test component.dielectric.shunt_capacitance ≈ imag(expected) / omega
            @test component.dielectric.shunt_conductance ≈ real(expected)
            path = export_data(:pscad, system, earth; base_freq=frequency,
                formulation=selected, temperature=60.0, file_name=joinpath(directory, "selected.pscx"))
            document = readxml(path)
            cable = only(findall("//User[@defn='master:Cable_Coax']", document))
            values = Dict(node["name"]=>node["value"] for node in findall("./paramlist/param", cable))
            @test parse(Float64, values["FLT"]) == frequency
            @test parse(Float64, values["RHOC"]) ≈ copper.rho * 1.16 rtol=1e-5
            @test parse(Float64, values["LT1"]) ≈ real(expected) / imag(expected) atol=5e-5
            @test parse(Float64, values["EPS1"]) ≈
                imag(expected) / omega * log(0.0065 / 0.004) / (2pi * epsilon0) rtol=1e-5
        end
        selected = Formulation()
        @test only(IE._pscad_components(design, frequency, selected, nothing)).conductor.material.rho == copper.rho
        uncorrected = Formulation(options=(temperature_correction=false,))
        @test only(IE._pscad_components(design, frequency, uncorrected, 60.0)).conductor.material.rho == copper.rho
        reference_shunt = inv(log(0.0045 / 0.004) / (2pi * admittivity(semicon, :Ametani2004, 20.0)) +
            log(0.0065 / 0.0045) / (2pi * admittivity(dielectric, :Ametani2004, 20.0)))
        for (temperature, correction) in ((nothing, true), (60.0, false))
            selected = Formulation(insulation_admittance=:Ametani2004,
                semicon_admittance=:Ametani2004, options=(temperature_correction=correction,))
            component = only(IE._pscad_components(design, frequency, selected, temperature))
            @test component.dielectric.shunt_conductance ≈ real(reference_shunt)
            @test component.dielectric.shunt_capacitance ≈ imag(reference_shunt) / omega
        end
        @test_throws ArgumentError export_data(:pscad, system, earth;
            temperature=60.0, file_name=joinpath(directory, "unspecified.pscx"))
        @test !isfile(joinpath(directory, "$(system.system_id)_unspecified.pscx"))
        frequency_error = try
            export_data(:pscad, system, earth;
                base_freq=Inf, file_name=joinpath(directory, "infinite.pscx"))
        catch exception
            exception
        end
        @test frequency_error isa DomainError
        @test occursin("base frequency must be positive and finite", sprint(showerror, frequency_error))
        @test !isfile(joinpath(directory, "$(system.system_id)_infinite.pscx"))

        lossy = build(CableDesign, "loss-cap", terminal(:core, solid(copper, Disk(0.004)),
            insulation(Material(:insulator, 1.0, 2.3); t=0.002)))
        component = only(IE._pscad_components(lossy, 50.0,
            Formulation(insulation_admittance=:Ametani2004), 20.0))
        requested = component.dielectric.shunt_conductance /
            (2pi * 50 * component.dielectric.shunt_capacitance)
        @test requested > 10
        emitted = Dict(IE._pscad_part_parameters(component, 1, 2pi * 50))
        @test parse(Float64, emitted["LT1"]) == 10
        @test component.dielectric.shunt_conductance /
            (2pi * 50 * component.dielectric.shunt_capacitance) == requested

        bare = build(CableDesign, "bare-wire", terminal(:core, solid(copper, Disk(0.0425))))
        bare_system = build(LineCableSystem, [bare], [Pose2(0, -1)];
            connections=[Dict(:core=>1)])
        for formulation in (nothing, Formulation())
            path = export_data(:pscad, bare_system, earth; formulation,
                file_name=joinpath(directory, "bare.pscx"))
            document = EzXML.readxml(path)
            cable = only(EzXML.findall("//User[@defn='master:Cable_Coax']", document))
            values = Dict(node["name"]=>node["value"] for node in EzXML.findall("./paramlist/param", cable))
            @test values["LL"] == "0"
            @test parse(Float64, values["R2"]) == parse(Float64, values["R3"]) == 0.0425
            _, imported = import_data(:pscad, path)
            @test length(only(imported.designs).geometry.regions) == 1
        end
        insulated_values = Dict(IE._pscad_cable_parameters(design, Pose2(0, -1), [1], 1, frequency))
        @test insulated_values["LL"] == "1"
    end
end
