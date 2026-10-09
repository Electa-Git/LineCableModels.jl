# Each consumer reduces a cable design to coaxial layers. The consumers are the engine
# blueprint, the reduction of the exporters, a homogenized design and the PSCAD components.
# They compute the same equivalents. These items check that the agreements remain while the
# reductions converge.

@testitem "DataModel / reductions / conductor equivalents agree across consumers" tags=[:integration, :pscad] setup=[TestFixtures] begin
    const DM = LineCableModels.DataModel
    const EN = LineCableModels.Engine
    const PS = LineCableModels.PSCAD
    designs = (TestFixtures.coaxial_design(),
        TestFixtures.internal_shunt_test_design(; tapes = false),
        TestFixtures.internal_shunt_test_design(), TestFixtures.three_core_design())
    for design in designs
        blueprint = EN.flatten(LineCableModelsCoaxial(), design, Float64)
        reduced = DM.flatten(design, 50.0)
        homogenized = DM.radial_components(homogenize(design), Float64)
        exported = PS._pscad_components(blueprint, 50.0, Formulation(), nothing)
        @test length(reduced) == length(homogenized) == length(exported) ==
              length(blueprint.conductors)
        for (index, row) in pairs(blueprint.conductors)
            @test row.terminal === reduced[index].name === homogenized[index].name ===
                  exported[index].name
            d, h, p = reduced[index].conductor, homogenized[index].conductor,
            exported[index].conductor
            @test row.r_in == d.r_in == h.r_in == p.r_in
            @test row.r_ex == d.r_ex == h.r_ex == p.r_ex
            @test row.alpha == d.alpha == h.alpha
            @test p.material.rho == row.material.rho
            @test p.material.mu_r == row.material.mu_r
            for field in (:rho, :mu_r, :alpha, :T0)
                expected = getfield(row.material, field)
                @test getfield(d.material, field) == expected
                @test getfield(h.material, field) ≈ expected rtol = 1e-14
            end
            # The reduction of a homogenized design starts from its equivalent material.
            for field in (:resistance, :gmr)
                @test getfield(d, field) == getfield(row, field)
                @test getfield(h, field) ≈ getfield(row, field) rtol = 1e-14
            end
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
        blueprint = EN.flatten(LineCableModelsCoaxial(), design, Float64)
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
        blueprint = EN.flatten(LineCableModelsCoaxial(), design, Float64)
        reduced = DM.flatten(design, frequency)
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
                # The exporters' reduction combines the layers with the lossy circuit.
                formulation === lossy || continue
                interval = reduced[index].dielectric
                @test complex(interval.shunt_conductance,
                    imag(s) * interval.shunt_capacitance) ≈ admittance rtol = 1e-12
                @test interval.material.eps_r ≈ dielectric.material.eps_r rtol = 1e-12
                @test interval.material.mu_r ≈ dielectric.material.mu_r rtol = 1e-14
                @test interval.r_in == engine.r_in[index]
                @test interval.r_ex == engine.r_ex[index]
            end
        end
    end
end
