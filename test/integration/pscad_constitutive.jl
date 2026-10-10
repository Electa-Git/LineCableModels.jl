@testitem "ImportExport / PSCAD / selected constitutive laws and operating temperature" tags=[:integration, :pscad] begin
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
    blueprint = LineCableModels.DataModel.flatten(design)
    # The admittance of the core's dielectric group: its layers in series, each under the
    # law of its class.
    function group_admittance(selected, temperature)
        laws = (selected.methods.insulation_admittance, selected.methods.semicon_admittance)
        inv(sum(blueprint.dielectrics) do layer
            kappa = constitutive(laws, layer.material, frequency,
                temperature === nothing ? layer.material.T0 : temperature;
                temperature_dependence = selected.methods.temperature_dependence)
            log(layer.r_ex / layer.r_in) / (2pi * kappa)
        end)
    end
    # The parameters of the one Cable_Coax of an export.
    function exported_values(path)
        cable = only(findall("//User[@defn='master:Cable_Coax']", readxml(path)))
        Dict(node["name"]=>node["value"] for node in findall("./paramlist/param", cable))
    end

    mktempdir() do directory
        for constructor in (Formulation, CableConstantsFormulation),
            insulation_id in (:default, :lossy), semicon_id in (:default, :lossy)
            selected = constructor(insulation_admittance=insulation_id,
                semicon_admittance=semicon_id)
            @test constitutive(selected.methods.temperature_dependence,
                only(blueprint.conductors).material, 60.0) ≈ copper.rho * 1.16
            expected = inv(log(0.0045 / 0.004) / (2pi * admittivity(semicon, semicon_id)) +
                log(0.0065 / 0.0045) / (2pi * admittivity(dielectric, insulation_id)))
            @test group_admittance(selected, 60.0) ≈ expected
            path = export_data(:pscad, system, earth; base_freq=frequency,
                formulation=selected, temperature=60.0, file_name=joinpath(directory, "selected.pscx"))
            values = exported_values(path)
            @test parse(Float64, values["FLT"]) == frequency
            @test parse(Float64, values["RHOC"]) ≈ copper.rho * 1.16 rtol=1e-5
            @test parse(Float64, values["LT1"]) ≈ real(expected) / imag(expected) atol=5e-5
            @test parse(Float64, values["EPS1"]) ≈
                imag(expected) / omega * log(0.0065 / 0.004) / (2pi * epsilon0) rtol=1e-5
        end
        # The reference temperature, or no temperature law, keeps the reference resistivity.
        for (selected, temperature) in ((Formulation(), nothing),
            (Formulation(temperature_dependence=nothing), 60.0))
            path = export_data(:pscad, system, earth; base_freq=frequency,
                formulation=selected, temperature, file_name=joinpath(directory, "reference.pscx"))
            @test parse(Float64, exported_values(path)["RHOC"]) == round(copper.rho; sigdigits=6)
        end
        reference_shunt = inv(log(0.0045 / 0.004) / (2pi * admittivity(semicon, :lossy, 20.0)) +
            log(0.0065 / 0.0045) / (2pi * admittivity(dielectric, :lossy, 20.0)))
        for (temperature, correction) in ((nothing, true), (60.0, false))
            selected = Formulation(insulation_admittance=:lossy,
                semicon_admittance=:lossy, temperature_dependence=correction ? formula(:default) : nothing)
            @test group_admittance(selected, temperature) ≈ reference_shunt
        end
        default_path = export_data(:pscad, system, earth;
            temperature=60.0, file_name=joinpath(directory, "unspecified.pscx"))
        @test isfile(default_path)
        default_document = readxml(default_path)
        default_cable = only(findall("//User[@defn='master:Cable_Coax']", default_document))
        default_values = Dict(node["name"]=>node["value"] for node in findall("./paramlist/param", default_cable))
        @test parse(Float64,default_values["LT1"]) == 0
        @test parse(Float64,default_values["RHOC"]) ≈ copper.rho * 1.16 rtol=1e-5
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
        lossy_blueprint = LineCableModels.DataModel.flatten(lossy)
        selected = Formulation(insulation_admittance=:lossy)
        part = LineCableModels.PSCAD._pscad_part_parameters(lossy_blueprint, 1, 1, 50.0,
            selected, 20.0)
        @test part.loss > 10
        @test parse(Float64, Dict(part.parameters)["LT1"]) == 10

        bare = build(CableDesign, "bare-wire", terminal(:core, solid(copper, Disk(0.0425))))
        bare_system = build(LineCableSystem, [bare], [Pose2(0, -1)];
            connections=[Dict(:core=>1)])
        for keywords in ((;), (formulation=Formulation(),))
            path = export_data(:pscad, bare_system, earth; keywords...,
                file_name=joinpath(directory, "bare.pscx"))
            document = EzXML.readxml(path)
            cable = only(EzXML.findall("//User[@defn='master:Cable_Coax']", document))
            values = Dict(node["name"]=>node["value"] for node in EzXML.findall("./paramlist/param", cable))
            @test values["LL"] == "0"
            @test parse(Float64, values["R2"]) == parse(Float64, values["R3"]) == 0.0425
            _, imported = import_data(:pscad, path)
            @test length(only(imported.designs).geometry.regions) == 1
        end
        insulated_values = Dict(LineCableModels.PSCAD._pscad_cable_parameters(design,
            blueprint, 1, Pose2(0, -1), [1], 1, frequency, Formulation(), nothing))
        @test insulated_values["LL"] == "1"
    end
end
