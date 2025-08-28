@testitem "Validation(Tubular): rule order unit test" setup = [defaults] begin
    # Use fully-qualified names; do not add extra `using` here.
    V = LineCableModels.Validation
    T = LineCableModels.DataModel.Tubular
    M = LineCableModels.Materials.Material

    r = V._rules(T)

    expected = (
        V.Normalized(:radius_in), V.Normalized(:radius_ext),
        V.Finite(:radius_in), V.Nonneg(:radius_in),
        V.Finite(:radius_ext), V.Nonneg(:radius_ext),
        V.Less(:radius_in, :radius_ext),
        V.Finite(:temperature),
        V.IsA{M}(:material_props),
    )

    if r != expected
        @error "[Validation] Rule set for Tubular is wrong. Someone ‘helpfully’ changed the bundle order or duplicated rules.\n" *
               "Expected exact structural equality with the generated bundle. Fix your traits/extra_rules and stop being clever."
        @show expected
        @show r
    end

    @test r == expected
end
