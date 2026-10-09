# Each consumer reduces a cable design to coaxial layers. The consumers are the engine
# blueprint, which the ATP and TRALIN exporters also read, a homogenized design and the
# PSCAD components. They compute the same equivalents. These items check that the
# agreements remain while the reductions converge.

@testitem "DataModel / reductions / conductor equivalents agree across consumers" tags=[:integration, :pscad] setup=[TestFixtures] begin
    const DM = LineCableModels.DataModel
    const EN = LineCableModels.Engine
    const PS = LineCableModels.PSCAD
    designs = (TestFixtures.coaxial_design(),
        TestFixtures.internal_shunt_test_design(; tapes = false),
        TestFixtures.internal_shunt_test_design(), TestFixtures.three_core_design())
    for design in designs
        blueprint = LineCableModels.DataModel.flatten(design, Float64)
        homogenized = DM.flatten(homogenize(design), Float64)
        exported = PS._pscad_components(blueprint, 50.0, Formulation(), nothing)
        @test length(homogenized) == length(exported) == length(blueprint.conductors)
        @test homogenized.mu_r_ins ≈ blueprint.mu_r_ins rtol = 1e-14
        for (index, row) in pairs(blueprint.conductors)
            h, p = homogenized.conductors[index], exported[index].conductor
            @test row.terminal === h.terminal === exported[index].name
            @test row.r_in == h.r_in == p.r_in
            @test row.r_ex == h.r_ex == p.r_ex
            @test p.material.rho == row.material.rho
            @test p.material.mu_r == row.material.mu_r
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

@testitem "DataModel / reductions / dielectric groups agree within each law" tags=[:integration, :pscad] setup=[TestFixtures] begin
    const DM = LineCableModels.DataModel
    const EN = LineCableModels.Engine
    const PS = LineCableModels.PSCAD
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
            exported = PS._pscad_components(blueprint, frequency, formulation, nothing)
            for (index, component) in pairs(exported)
                dielectric = component.dielectric
                admittance = complex(dielectric.shunt_conductance,
                    imag(s) * dielectric.shunt_capacitance)
                @test engine.admittance[index] ≈ admittance rtol = 1e-12
                @test again.admittance[index] ≈ admittance rtol = 1e-12
                @test engine.mu_r[index] ≈ dielectric.material.mu_r rtol = 1e-14
                @test again.mu_r[index] ≈ dielectric.material.mu_r rtol = 1e-14
                @test engine.r_ex[index] == dielectric.r_ex
                # The ATP and TRALIN exporters combine the layers with the lossy circuit.
                formulation === lossy || continue
                conductor = blueprint.conductors[index]
                layers = blueprint.dielectrics[blueprint.dielectric_ranges[index]]
                circuit = DM.dielectric_circuit(layers, imag(s), Float64)
                interval = (r_in = isempty(layers) ? conductor.r_ex : first(layers).r_in,
                    r_ex = isempty(layers) ? conductor.r_ex : last(layers).r_ex, layers,
                    shunt_capacitance = circuit.capacitance,
                    shunt_conductance = circuit.conductance)
                material = DM.equivalent_dielectric_material(interval,
                    blueprint.mu_r_ins[index], conductor.material.T0, imag(s))
                @test complex(interval.shunt_conductance,
                    imag(s) * interval.shunt_capacitance) ≈ admittance rtol = 1e-12
                @test material.eps_r ≈ dielectric.material.eps_r rtol = 1e-12
                @test material.mu_r ≈ dielectric.material.mu_r rtol = 1e-14
                @test interval.r_in == engine.r_in[index]
                @test interval.r_ex == engine.r_ex[index]
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
