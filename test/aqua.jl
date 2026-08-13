@testitem "Aqua tests" tags=[:skipci] begin
	using Aqua
	Aqua.test_all(LineCableModels)
end

