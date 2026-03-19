@testitem "Validation(Tubular): rule order unit test" setup = [defaults] begin
	# Use fully-qualified names; do not add extra `using` here.
	V = LineCableModels.Validation
	T = LineCableModels.DataModel.Tubular
	M = LineCableModels.Materials.Material

	r = V._rules(T)

	expected = (
		V.Normalized(:r_in), V.Normalized(:r_ex),
		V.Finite(:r_in), V.Nonneg(:r_in),
		V.Finite(:r_ex), V.Nonneg(:r_ex),
		V.Less(:r_in, :r_ex),
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
