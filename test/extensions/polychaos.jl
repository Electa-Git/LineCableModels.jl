@testitem "Extensions / PolyChaos boundary / activation" tags = [:extension] begin
    import LineCableModels
    @test Base.get_extension(
        LineCableModels, :LineCableModelsPolyChaosExt) === nothing

    using PolyChaos

    extension = Base.get_extension(LineCableModels, :LineCableModelsPolyChaosExt)
    @test extension !== nothing
    @test pkgversion(PolyChaos) == v"2.1.0"
    @test any(
        method -> method.module === extension,
        methods(LineCableModels.compute)
    )
end
