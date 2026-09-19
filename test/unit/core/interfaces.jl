@testitem "Core / interfaces / domain and deterministic uncertainty fallbacks" tags=[:unit] begin
    @test domain(Int) === nothing
    @test domain(1) === nothing
    @test PhaseDomain() isa LineCableModels.LineParamsDomain
    @test_throws MethodError ModalDomain()
    @test PhaseDomain !== ModalDomain

    @test nominal(3.0 + 4.0im) == 3.0 + 4.0im
    @test nominal([1.0, 2.0]) == [1.0, 2.0]
    @test uncertainty(3.0) == 0.0
    @test uncertainty("not numeric") == 0.0
    value=(answer = 42,)
    @test_throws MethodError validate(value)
end

@testitem "Core / docstrings / sanitized method-list provenance" tags=[:unit] begin
    using DocStringExtensions

    method=which(LineCableModels.domain, (Int,))
    @test !isabspath(LineCableModels._method_path(method))
    @test !occursin(pkgdir(LineCableModels), LineCableModels._method_path(method))

    binding=Docs.Binding(LineCableModels.Engine, :domain)
    doc=first(values(Docs.meta(LineCableModels.Engine)[binding].docs))
    buffer=IOBuffer()
    @test DocStringExtensions.format(LineCableModels.METHODLIST, buffer, doc) === nothing
    rendered=String(take!(buffer))
    @test occursin("domain", rendered)
    @test occursin(LineCableModels._method_path(method), rendered)
    @test !occursin(pkgdir(LineCableModels), rendered)
end

@testitem "Core / owner-local numerics / transforms and conductivity" tags=[:unit] begin
    using LinearAlgebra
    const Engine=LineCableModels.Engine
    const MatrixOps=Engine

    matrix=[Float64(2i+3j+i*j) for i in 1:3, j in 1:3]
    circulant=copy(matrix)
    @test MatrixOps.ideal_transposition!(circulant) === circulant
    @test all(circulant[i, j] == circulant[mod1(i + 1, 3), mod1(j + 1, 3)]
    for i in 1:3, j in 1:3)
    @test sum(circulant) ≈ sum(matrix)
    @test_throws DimensionMismatch MatrixOps.ideal_transposition!(ones(2, 3))
    @test LineCableModels.Transforms.offdiagonal_ratio(Diagonal([1.0, 2.0])) == 0.0
    @test_throws DimensionMismatch LineCableModels.Transforms.offdiagonal_ratio(zeros(2, 3))

    @test Engine.conductivity(Inf) == 0.0
    @test isinf(Engine.conductivity(0.0))
    @test Engine.conductivity(4.0) == 0.25
end
