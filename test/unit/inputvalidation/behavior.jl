@testitem "InputValidation / materialized inputs / actionable failures" tags=[:unit] setup=[
    TestFixtures
] begin
    using RequiredInterfaces

    const IV=LineCableModels.InputValidation

    material=Material(
        kind = :conductor,
        rho = 1.7241e-8,
        eps_r = 1.0,
        mu_r = 1.0,
        T0 = 20.0,
        alpha = 0.00393
    )
    @test validate(material) === material
    @test parentmodule(validate) === IV
    @test_throws MethodError validate((answer = 42,))

    struct UnvalidatedEarthDescription<:AbstractEarthModel end
    @test_throws RequiredInterfaces.NotImplementedError validate(
        UnvalidatedEarthDescription()
    )

    material_error=try
        Material(
            kind = :conductor,
            rho = -1.0,
            eps_r = 1.0,
            mu_r = 1.0,
            T0 = 20.0,
            alpha = 0.0
        )
        nothing
    catch error
        error
    end
    @test material_error isa DomainError
    @test occursin("Material.rho", sprint(showerror, material_error))
    @test occursin("-1.0", sprint(showerror, material_error))

    layer_error=try
        layer(rho = 100.0, thickness = 0.0)
        nothing
    catch error
        error
    end
    @test layer_error isa DomainError
    @test occursin("EarthLayer.thickness", sprint(showerror, layer_error))
    @test occursin("0.0", sprint(showerror, layer_error))

    design=TestFixtures.mv_cable_design()
    @test validate(design) === design
    damaged_design=TestFixtures.mv_cable_design()
    pop!(damaged_design.terminal_map)
    design_error=try
        validate(damaged_design)
        nothing
    catch error
        error
    end
    @test design_error isa DimensionMismatch
    @test occursin("CableDesign.terminal_map", sprint(showerror, design_error))

    system=TestFixtures.three_phase_system()
    @test validate(system) === system
    system.connection_order[1]=99
    system_error=try
        validate(system)
        nothing
    catch error
        error
    end
    @test system_error isa DimensionMismatch
    @test occursin("LineCableSystem.connection_order[1]", sprint(showerror,
        system_error))

    struct CheckedProblem<:AbstractProblemDefinition
        value::Int
    end
    function LineCableModels.validate(problem::CheckedProblem)
        problem.value>0||throw(DomainError(
            problem.value,
            "CheckedProblem.value must be positive"
        ))
        return problem
    end
    space=Gridspace{CheckedProblem}(CheckedProblem, (Grid((1, -1)),))
    points=LineCableModels.points(space)
    first_point, state=iterate(points)
    @test LineCableModels.materialize(first_point) == CheckedProblem(1)
    second_point=first(iterate(points, state))
    point_error=try
        LineCableModels.materialize(second_point)
        nothing
    catch error
        error
    end
    @test point_error isa DomainError
    @test occursin("CheckedProblem.value", sprint(showerror, point_error))
end
