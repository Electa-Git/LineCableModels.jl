@testitem "Engine / homogenized dielectric layers retain constitutive selection across frequency" tags=[:unit, :measurements] begin
    using Test
    using LineCableModels
    using Measurements: measurement, value, uncertainty
    const E = LineCableModels.Engine
    const DM = LineCableModels.DataModel
    const IE = LineCableModels.ImportExport
    using EzXML: findall, readxml
    using Logging: NullLogger, with_logger
    # The capacitance and conductance of the core's dielectric group, read back from the
    # ATP export of a design alone.
    function exported_group(design, formulation, frequency, temperature)
        system = build(LineCableSystem, [design], [(0.0, -1.0)];
            connections = [Dict(:core => 1)])
        return mktempdir() do directory
            path = with_logger(NullLogger()) do
                export_data(:atp, system, EarthModel(100.0, 10.0, 1.0); formulation,
                    base_freq = frequency, temperature,
                    file_name = joinpath(directory, "group.xml"))
            end
            document = readxml(path)
            core = first(findall("//conductor", document))
            outer = parse(Float64, first(findall("//cable", document))["Rout"])
            (C = 2π * LineCableModels.Commons.vacuum_permittivity(Float64) *
                 parse(Float64, core["epsI"]) / log(outer / parse(Float64, core["Rout"])),
                G = parse(Float64, core["Gext"]))
        end
    end
    metal = Material(:conductor, 1.72e-8)
    semicon = Material(:semicon, 1e5, 1000.0, 1.0, 20.0, 0.01; tan_delta = 0.02)
    dielectric = Material(:insulator, 2e9, 3.1, 1.0, 20.0, 0.02; tan_delta = 0.04)
    radii = (0.01, 0.012, 0.02, 0.021)
    source = @cable "radial-mixture" begin
        @terminal :core begin
            solid(metal, Disk(radii[1]))
            screen(semicon; t = radii[2]-radii[1])
            insulation(dielectric; t = radii[3]-radii[2])
            screen(semicon; t = radii[4]-radii[3])
        end
    end
    # homogenize reduces the conductor and keeps each dielectric layer at its radii.
    homogeneous = homogenize(source)
    @test length(homogeneous.geometry.regions) == 4
    layers = [region.source.material for region in homogeneous.geometry.regions[2:end]]
    @test all(material -> material isa Material, layers)
    # A solid core has no lay, and the permeability scale is one up to rounding.
    for (actual, expected) in zip(layers, (semicon, dielectric, semicon))
        @test actual.kind === expected.kind
        for field in (:rho, :eps_r, :T0, :alpha, :tan_delta)
            @test getproperty(actual, field) == getproperty(expected, field)
        end
        @test actual.mu_r ≈ expected.mu_r rtol=1e-15
    end
    blueprint = @inferred(LineCableModels.DataModel.flatten(homogeneous, Float64))
    @test blueprint isa LineCableModels.DataModel.CableBlueprint{Float64}
    original_blueprint = LineCableModels.DataModel.flatten(source)
    @test getproperty.(blueprint.dielectrics, :r_in) ==
          getproperty.(original_blueprint.dielectrics, :r_in)
    @test getproperty.(blueprint.dielectrics, :r_ex) ==
          getproperty.(original_blueprint.dielectrics, :r_ex)
    @test blueprint.mu_r_ins ≈ original_blueprint.mu_r_ins rtol=1e-15
    @test only(DM.preview_materials(last(homogeneous.geometry.regions))) == last(layers)
    @test length(DM.preview_shapes(last(homogeneous.geometry.regions))) == 1
    reconstructed = IE.deserialize_value(IE.serialize_value(homogeneous))
    @test [region.source.material for region in reconstructed.geometry.regions[2:end]] ==
          layers

    frequencies = 10.0 .^ range(-1, 6; length = 31)
    for ins in (:default, :lossy), semi in (:default, :lossy)

        formulation = CableConstantsFormulation(insulation_admittance = ins,
            semicon_admittance = semi)
        relations = (formulation.methods.insulation_admittance,
            formulation.methods.semicon_admittance)
        for f in frequencies
            original,
            reduced = map((source, homogeneous)) do design
                blueprint = LineCableModels.DataModel.flatten(design)
                input = E.LocalCableData(blueprint)
                shunt = only(E.internal_shunt_response([design], [blueprint], [formulation]))
                y = zeros(ComplexF64, 1, 1)
                admittivity = zeros(ComplexF64, length(input.dielectric_materials))
                rho = [constitutive(formulation.methods.temperature_dependence, material, 20.0)
                       for material in input.dielectric_materials]
                E.dielectric!(admittivity, input, rho, formulation.methods, f, 20.0)
                E.cable_admittance!(y, input, shunt, admittivity,
                    complex(0.0, 2π*f), zeros(ComplexF64, length(input.dielectric_materials)))
                only(y)
            end
            @test imag(reduced) ≈ imag(original) rtol=2e-13
            @test real(reduced) ≈ real(original) rtol=2e-13
            # The layers in series, each under the law of its class.
            y = inv(sum(zip(layers, radii[1:(end - 1)], radii[2:end])) do (m, a, b)
                law = m.kind === :semicon ? last(relations) : first(relations)
                log(b / a) / (2π * @inferred(constitutive(law, m, f, 20.0)))
            end)
            @test imag(reduced) ≈ imag(y) rtol=2e-13
            @test real(reduced) ≈ real(y) rtol=2e-13
            if ins === semi === :default
                @test real(reduced) == 0
            end
        end
        for f in (50.0, 60.0)
            original = compute(CableConstantsProblem(source; frequency = f), formulation)
            reduced = compute(CableConstantsProblem(homogeneous; frequency = f), formulation)
            @test reduced.C ≈ original.C rtol=2e-13
            @test reduced.G ≈ original.G rtol=2e-13
            # The complex admittance agrees in its last bits. G alone is its small real
            # part and carries the absolute error of the whole.
            exported = exported_group(homogeneous, formulation, f, 20.0)
            @test complex(exported.G, 2π * f * exported.C) ≈
                  complex(only(original.G), 2π * f * only(original.C)) rtol=2e-13
        end
    end

    # Independent series circuit, not a baseline regenerated by the implementation.
    ε0 = 8.8541878128e-12
    weights = log.(collect(radii[2:end]) ./ collect(radii[1:(end - 1)]))
    physical = (semicon, dielectric, semicon)
    for correction in (false, true), f in (50.0, 60.0)

        formulation = CableConstantsFormulation(insulation_admittance = :lossy,
            semicon_admittance = :lossy, temperature_dependence = correction ?
                                                                  formula(:default) :
                                                                  nothing)
        original = compute(CableConstantsProblem(source; frequency = f, temperature = 60.0), formulation)
        reduced = compute(CableConstantsProblem(homogeneous; frequency = f, temperature = 60.0), formulation)
        @test original.C ≈ reduced.C rtol=2e-13
        @test original.G ≈ reduced.G rtol=2e-13
        expected = inv(sum(q / (2π * (inv(m.rho * (correction ? 1+m.alpha*40 : 1)) +
                             2π*f*ε0*m.eps_r*(m.tan_delta+im)))
        for (q, m) in zip(weights, physical)))
        @test only(original.G) ≈ real(expected) rtol=2e-13
        @test only(original.C) ≈ imag(expected)/(2π*f) rtol=2e-13
        @test only(reduced.G) ≈ real(expected) rtol=2e-13
        @test only(reduced.C) ≈ imag(expected)/(2π*f) rtol=2e-13
        exported = exported_group(homogeneous, formulation, f, 60.0)
        @test complex(exported.G, 2π * f * exported.C) ≈ expected rtol=2e-13
    end
    @test semicon.T0 == dielectric.T0 == 20.0
    other = Material(:insulator, 4e9, 2.5, 1.0, 20.0, 0.02; tan_delta = 0.01)
    regions = Region(:equivalent, Annulus(0.01, 0.02), Grid((dielectric, other)))
    @test regions isa Gridspace{Region}
    @test length(regions) == 2
    @test isconcretetype(eltype(regions))
    @test @inferred(first(regions)) isa Region
    transported = Gridspace{CableDesign}(homogenize, (Grid((source, source)),))
    @test length(transported) == 2
    @test first(transported) isa CableDesign
    @test last(first(transported).geometry.regions).source.material == last(layers)
    uncertain = @cable "uncertain-radial-mixture" begin
        @terminal :core begin
            solid(metal, Disk(0.01))
            insulation(Material(:insulator, 1e10, measurement(3.1, 0.1)); t = 0.005)
            insulation(dielectric; t = 0.005)
        end
    end
    measured = homogenize(uncertain)
    before = compute(CableConstantsProblem(uncertain), CableConstantsFormulation())
    after = compute(CableConstantsProblem(measured), CableConstantsFormulation())
    @test value.(before.C) ≈ value.(after.C)
    @test uncertainty.(before.C) ≈ uncertainty.(after.C)
    @test only(uncertainty.(after.C)) > 0
    roundtrip = IE.deserialize_value(IE.serialize_value(measured))
    replay = compute(CableConstantsProblem(roundtrip), CableConstantsFormulation())
    @test value.(after.C) ≈ value.(replay.C)
    @test uncertainty.(after.C) ≈ uncertainty.(replay.C)
end

@testitem "ImportExport / a homogenized design library retains the dielectric layers" tags=[:unit, :importexport] begin
    using JSON3
    using JSONSchema
    const IE = LineCableModels.ImportExport
    metal = Material(:conductor, 1.72e-8)
    semicon = Material(:semicon, 1e5, 1000.0; tan_delta = 0.02)
    dielectric = Material(:insulator, 2e9, 3.1; tan_delta = 0.04)
    source = @cable "stored-radial-mixture" begin
        @terminal :core begin
            solid(metal, Disk(0.01))
            screen(semicon; t = 0.002)
            insulation(dielectric; t = 0.008)
        end
    end
    design = homogenize(source)
    library = CablesLibrary()
    add!(library, design)
    schema_path = joinpath(pkgdir(LineCableModels), "schemas", "cable-design.schema.json")
    schema = JSONSchema.Schema(JSON3.read(read(schema_path, String), Dict{String, Any}))
    mktempdir() do directory
        path = joinpath(directory, "homogeneous.json")
        save(library; file_name = path)
        document = JSON3.read(read(path, String), Dict{String, Any})
        @test JSONSchema.validate(schema, document) === nothing
        restored = CablesLibrary()
        load!(restored; file_name = path)
        actual = restored[design.cable_id]
        @test IE.serialize_value(actual) == IE.serialize_value(design)
        materials = [region.source.material for region in actual.geometry.regions[2:end]]
        @test getproperty.(materials, :kind) == [:semicon, :insulator]
        for selector in (:default, :lossy), frequency in (50.0, 60.0)

            formulation = CableConstantsFormulation(insulation_admittance = selector,
                semicon_admittance = selector)
            expected = compute(CableConstantsProblem(design; frequency), formulation)
            replay = compute(CableConstantsProblem(actual; frequency), formulation)
            @test replay.C == expected.C
            @test replay.G == expected.G
        end
    end
end
