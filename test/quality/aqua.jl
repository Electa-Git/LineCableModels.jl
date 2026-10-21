@testitem "Quality / Aqua / package hygiene" tags = [:aqua, :pscad] begin
    using Aqua
    Aqua.test_all(LineCableModels)
end
