using Test
using Aqua
using LineCableModels

@testset "LineCableModels" begin
    @testset "Aqua" begin
        Aqua.test_all(LineCableModels)
    end
    include("architecture/architecture.jl")
end
