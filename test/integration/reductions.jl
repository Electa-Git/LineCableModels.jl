# A cable design reduces to one blueprint of coaxial layers. The engine and the ATP, PSCAD
# and TRALIN exporters read it, and homogenize builds a new design from it. These items
# check that the consumers agree: the homogenized design reduces to the same equivalents,
# and each export reproduces the cable constants of the engine.

@testitem "DataModel / reductions / conductor equivalents agree across consumers" tags=[:integration, :parametric] setup=[TestFixtures] begin
    const DM = LineCableModels.DataModel
    designs = (TestFixtures.coaxial_design(),
        TestFixtures.internal_shunt_test_design(; tapes = false),
        TestFixtures.internal_shunt_test_design(), TestFixtures.three_core_design())
    for design in designs
        blueprint = LineCableModels.DataModel.flatten(design, Float64)
        homogenized = DM.flatten(homogenize(design), Float64)
        @test length(homogenized) == length(blueprint.conductors)
        @test homogenized.mu_r_ins ≈ blueprint.mu_r_ins rtol = 1e-14
        for (index, row) in pairs(blueprint.conductors)
            h = homogenized.conductors[index]
            @test row.terminal === h.terminal
            @test row.r_in == h.r_in
            @test row.r_ex == h.r_ex
            for field in (:rho, :mu_r, :alpha, :T0)
                @test getfield(h.material, field) ≈ getfield(row.material, field) rtol = 1e-14
            end
            # The reduction of a homogenized design starts from its equivalent material.
            @test DM.tubular_resistance(h.r_in, h.r_ex, h.material.rho) ≈
                  DM.tubular_resistance(row.r_in, row.r_ex, row.material.rho) rtol = 1e-14
            @test DM.tubular_gmr(h.r_ex, h.r_in, h.material.mu_r) ≈
                  DM.tubular_gmr(row.r_ex, row.r_in, row.material.mu_r) rtol = 1e-14
        end
    end
end

@testitem "DataModel / reductions / dielectric groups agree within each law" tags=[:integration, :parametric] setup=[TestFixtures] begin
    const EN = LineCableModels.Engine
    frequency = 50.0
    s = complex(0.0, 2π * frequency)
    # Read the admittance of each interval, in S/m, from the nodal matrix of the engine at
    # 20 °C. The permeability and radii of the interval come with it.
    function intervals(design, formulation)
        blueprint = LineCableModels.DataModel.flatten(design, Float64)
        input = EN.LocalCableData(blueprint)
        shunt = only(EN.internal_shunt_response([design], [blueprint], [formulation]))
        admittivity = Vector{ComplexF64}(undef, length(input.dielectric_materials))
        EN.dielectric!(admittivity, input, formulation.methods, frequency, 20.0)
        count = length(input.terminals)
        nodal = EN.cable_admittance!(Matrix{ComplexF64}(undef, count, count), input, shunt,
            admittivity, s, similar(admittivity))
        admittance = zeros(ComplexF64, count)
        for conductors in input.assemblies, (position, index) in pairs(conductors)
            admittance[index] = position < length(conductors) ?
                                -nodal[index, conductors[position + 1]] :
                                nodal[index, index] +
                                (position > 1 ? nodal[conductors[position - 1], index] : 0)
        end
        return (; admittance, mu_r = input.mu_r_ins, r_in = input.r_ins_in,
            r_ex = input.r_ins_ext)
    end
    lossless = Formulation()
    lossy = Formulation(insulation_admittance = :lossy, semicon_admittance = :lossy)
    designs = (TestFixtures.coaxial_design(),
        TestFixtures.internal_shunt_test_design(; tapes = false),
        TestFixtures.internal_shunt_test_design(), TestFixtures.three_core_design())
    for design in designs
        blueprint = LineCableModels.DataModel.flatten(design, Float64)
        homogenized = homogenize(design)
        for formulation in (lossless, lossy)
            engine = intervals(design, formulation)
            again = intervals(homogenized, formulation)
            laws = (formulation.methods.insulation_admittance,
                formulation.methods.semicon_admittance)
            # The exporters put the layers of each group in series, each under its law.
            for (index, row) in enumerate(blueprint)
                (; conductor, layers, mu_r_ins) = row
                isempty(layers) && continue
                admittance = inv(sum(layers) do layer
                    log(layer.r_ex / layer.r_in) /
                    (2π * constitutive(laws, layer.material, frequency, 20.0;
                        temperature_dependence = formulation.methods.temperature_dependence))
                end)
                @test engine.admittance[index] ≈ admittance rtol = 1e-12
                @test again.admittance[index] ≈ admittance rtol = 1e-12
                @test engine.mu_r[index] ≈ mu_r_ins rtol = 1e-14
                @test again.mu_r[index] ≈ mu_r_ins rtol = 1e-14
                @test first(layers).r_in == conductor.r_ex == engine.r_in[index]
                @test last(layers).r_ex == engine.r_ex[index]
            end
        end
    end
end

@testitem "DataModel / reductions / a terminal of mixed materials combines alpha to first order" tags=[:integration, :parametric] setup=[TestFixtures] begin
    const DM = LineCableModels.DataModel
    const EN = LineCableModels.Engine
    # The screen terminal of the 18 kV cable holds 49 copper wires and an aluminum tape.
    design = TestFixtures.screened_cable_design(tape_material = :aluminum)
    blueprint = LineCableModels.DataModel.flatten(design, Float64)
    row = only(filter(conductor -> conductor.terminal === :sheath, blueprint.conductors))
    members = filter(region -> region.terminal === :sheath, design.geometry.regions)
    alpha = [region.source.material.alpha for region in members]
    @test length(unique(alpha)) == 2
    resistance = map(members) do region
        factor = prod(entry -> overlength(entry.path, entry.radius), region.paths;
            init = 1.0)
        region.source.material.rho * factor / DM.area(region.source.primitive)
    end
    weights = inv.(resistance) ./ sum(inv, resistance)
    mean = sum(weights .* alpha)
    equivalent = row.material.rho / (π * (row.r_ex^2 - row.r_in^2))
    @test equivalent ≈ inv(sum(inv, resistance)) rtol = 1e-12
    @test row.material.alpha ≈ mean rtol = 1e-12
    # The linear law exceeds the parallel of the corrected members by a second-order term.
    # The weighted variance of the member coefficients times dT^2 bounds this term. A wrong
    # coefficient changes the result at first order in dT.
    dT = 90.0 - row.material.T0
    exact = inv(sum(inv, resistance .* (1 .+ alpha .* dT)))
    linear = equivalent * (1 + row.material.alpha * dT)
    variance = sum(weights .* (alpha .- mean) .^ 2)
    @test -1e-12 <= linear / exact - 1 <= variance * dT^2
end

@testitem "DataModel / reductions / each export reproduces the cable constants of its cores" tags=[:integration, :pscad] setup=[TestFixtures] begin
    using EzXML
    # The item reads each format back and forms the capacitance and conductance of each core
    # group as that format does. It compares the admittance of each group, G + jωC, with
    # CableConstants.
    # The formats write the resistivity of every conductor too.
    attribute(node, name) = parse(Float64, node[name])
    epsilon0 = LineCableModels.Commons.vacuum_permittivity(Float64)
    capacitance(a, b, eps_r) = 2π * epsilon0 * eps_r / log(b / a)
    earth = homogeneous(rho = 100.0)
    temperature = 90.0
    omega = 2π * 50.0
    lossy = CableConstantsFormulation(insulation_admittance = :lossy,
        semicon_admittance = :lossy)
    designs = (TestFixtures.coaxial_design(),
        TestFixtures.internal_shunt_test_design(; tapes = false),
        TestFixtures.internal_shunt_test_design(), TestFixtures.screened_cable_design(),
        TestFixtures.three_core_design())
    mktempdir() do directory
        for design in designs, formulation in (CableConstantsFormulation(), lossy)
            constants = compute(CableConstantsProblem(design; temperature), formulation)
            blueprint = LineCableModels.DataModel.flatten(design)
            rho = [constitutive(formulation.methods.temperature_dependence,
                       row.conductor.material, temperature) for row in blueprint]
            expected = complex.(constants.G, omega .* constants.C)
            system = TestFixtures.single_cable_system(design)
            keywords = (; formulation, temperature)
            @test length(constants.C) == length(blueprint.assembly_ranges)

            path = export_data(:atp, system, earth; keywords...,
                file_name = joinpath(directory, "export.xml"))
            cables = findall("//cable", readxml(path))
            @test length(cables) == length(blueprint.assembly_ranges)
            for (assembly, cable) in enumerate(cables)
                conductors = findall("conductor", cable)
                core = first(conductors)
                outer = length(conductors) > 1 ? attribute(conductors[2], "Rin") :
                        attribute(cable, "Rout")
                # ATPDraw adds Cext to the capacitance of epsI.
                C = capacitance(attribute(core, "Rout"), outer, attribute(core, "epsI")) +
                    attribute(core, "Cext")
                @test complex(attribute(core, "Gext"), omega * C) ≈ expected[assembly] rtol = 1e-12
                @test attribute(core, "Cext") == 0
                for (node, index) in zip(conductors, blueprint.assembly_ranges[assembly])
                    @test attribute(node, "rho") ≈ rho[index] rtol = 1e-12
                end
            end

            path = export_data(:pscad, system, earth; keywords...,
                file_name = joinpath(directory, "export.pscx"))
            cables = findall("//User[@defn='master:Cable_Coax']", readxml(path))
            @test length(cables) == length(blueprint.assembly_ranges)
            for (assembly, cable) in enumerate(cables)
                values = Dict(node["name"] => parse(Float64, node["value"])
                for node in findall("./paramlist/param", cable)
                if !isnothing(tryparse(Float64, node["value"])))
                # PSCAD receives G as ωC·LT, below the cap of ten.
                C = capacitance(values["R2"], values["R3"], values["EPS1"])
                @test C ≈ constants.C[assembly] rtol = 1e-4
                values["LT1"] < 10 &&
                    @test omega * C * complex(values["LT1"], 1) ≈ expected[assembly] rtol = 1e-4
                @test values["RHOC"] ≈ rho[first(blueprint.assembly_ranges[assembly])] rtol = 1e-5
            end

            path = export_data(:tralin, system, earth; keywords...,
                file_name = joinpath(directory, "export.f05"))
            lines = readlines(path)
            @test count(startswith("CABLE,"), lines) == length(blueprint.assembly_ranges)
            at = findall(line -> any(label -> startswith(line, label * ","),
                    ("CORE", "SHEATH", "ARMOUR")), lines)
            @test length(at) == length(blueprint.conductors)
            for (index, line) in enumerate(at)
                row = split(lines[line], ",")
                @test parse(Float64, row[6]) * 1.7241e-8 ≈ rho[index] rtol = 1e-5
                layers = blueprint[index].layers
                isempty(layers) && continue
                inner, outer = first(layers).r_in, last(layers).r_ex
                sigma = parse(Float64, row[8]) * 1e-14
                assembly = findfirst(conductors -> first(conductors) == index,
                    blueprint.assembly_ranges)
                # Fields 7 and 8 on the radii of the interval.
                assembly === nothing ||
                    @test complex(2π * sigma / log(outer / inner),
                        omega * capacitance(inner, outer, parse(Float64, row[9]))) ≈
                          expected[assembly] rtol = 1e-5
                # The INSULATION subcommand repeats the inline fields and adds mu_r_ins.
                insulation = split(lines[line + 1], ",")
                @test first(insulation) == "INSULATION"
                @test parse(Float64, insulation[2]) ≈ inv(max(sigma, 1e-20)) rtol = 1e-4
                @test parse(Float64, insulation[3]) ≈ outer - inner atol = 1e-6
                @test insulation[4] == row[9]
                @test parse(Float64, insulation[5]) ≈ blueprint[index].mu_r_ins rtol = 1e-5
            end
        end
    end
end
