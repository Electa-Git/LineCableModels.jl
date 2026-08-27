@testitem "Quality / explicit imports / package ownership" tags = [:quality] begin
    using ExplicitImports: test_explicit_imports

    test_explicit_imports(LineCableModels)
end
