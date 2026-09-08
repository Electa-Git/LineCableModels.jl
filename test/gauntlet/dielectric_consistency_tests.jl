@testitem "Gauntlet / 18 kV radial equivalence does not select dielectric losses" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using LinearAlgebra
    using LineCableModels.Engine
    using .GauntletSupport
    const E = LineCableModels.Engine
    model = GauntletSupport.reference_case(:cable_18kv_1000mm2_trefoil)
    source = first(model.nominal_problem.system.designs)
    reduced = homogenize(source)
    inputs = map(design -> E.LocalCableData(E.flatten(LineCableModelsCoaxial(), design)),
        (source, reduced))
    for identifier in (:default, :Ametani2004)
        formulation = Formulation(insulation_admittance=identifier, semicon_admittance=identifier)
        for frequency in model.nominal_problem.frequencies
            matrices = map(inputs) do input
                y = zeros(ComplexF64, length(input.terminals), length(input.terminals))
                E.cable_admittance!(y, input, formulation.methods, frequency,
                    model.nominal_problem.temperature, complex(0.0, 2π*frequency),
                    zeros(ComplexF64, length(input.dielectric_materials)))
            end
            @test matrices[1] ≈ matrices[2] rtol=1e-12
        end
        constants = CableConstantsFormulation(insulation_admittance=identifier,
            semicon_admittance=identifier)
        a, b = map(design -> compute(CableConstantsProblem(design), constants), (source, reduced))
        @test a.C ≈ b.C rtol=1e-12
        @test a.G ≈ b.G rtol=1e-12
        @test a.R ≈ b.R rtol=1e-10
        @test a.L ≈ b.L rtol=1e-10
    end
end
