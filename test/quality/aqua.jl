@testitem "Quality / Aqua / package hygiene" tags = [:quality] begin
    using Aqua
    Aqua.test_all(LineCableModels)
end
