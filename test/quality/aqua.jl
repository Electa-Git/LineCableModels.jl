@testitem "Quality / Aqua / package hygiene" tags = [:aqua, :importexport] begin
    using Aqua
    Aqua.test_all(LineCableModels)
end
