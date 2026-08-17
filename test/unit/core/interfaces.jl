@testitem "Core / interfaces / domain and deterministic uncertainty fallbacks" tags=[:unit] begin
    @test domain(Int) === nothing
    @test domain(1) === nothing
    @test PhaseDomain() isa LineCableModels.LineParamsDomain
    @test ModalDomain() isa LineCableModels.LineParamsDomain
    @test PhaseDomain !== ModalDomain

    @test nominal(3.0 + 4.0im) == 3.0 + 4.0im
    @test nominal([1.0, 2.0]) == [1.0, 2.0]
    @test standard_uncertainty(3.0) == 0.0
    @test standard_uncertainty("not numeric") == 0.0
    value=(answer = 42,)
    @test validate(value) === value
end

@testitem "Core / retired features / actionable tombstones" tags=[:unit] begin
    const DM=LineCableModels.DataModel
    const Engine=LineCableModels.Engine

    retired_calls=(
        () -> Engine.Formulation(Val(:FEM)),
        () -> Engine.FEM.Darwin(),
        () -> Engine.FEM.Electrodynamics(),
        () -> Engine.FEM.MeshTransition(),
        () -> Engine.FEM.calc_domain_size(),
        () -> Engine.FEM.preview_results(),
        () -> DM.SectorParams(),
        () -> DM.Sector(),
        () -> DM.SectorInsulator()
    )
    for invoke in retired_calls
        exception=try
            invoke()
            nothing
        catch caught
            caught
        end
        @test exception isa ArgumentError
        message=sprint(showerror, exception)
        @test occursin("was removed from LineCableModels", message)
        @test occursin("legacy/fem-sector", message)
        @test occursin("Pkg.add", message)
    end

    exception=try
        LineCableModels.retired_legacy_json()
        nothing
    catch caught
        caught
    end
    @test exception isa ArgumentError
    @test occursin("a71bdfe1ac832f27a0c88b1d02596194aac46ec7", sprint(showerror, exception))
end

@testitem "Core / docstrings / sanitized method-list provenance" tags=[:unit] begin
    using DocStringExtensions

    method=which(LineCableModels.domain, (Int,))
    @test LineCableModels._method_path(method) == joinpath("src", "interfaces.jl")
    @test !occursin(pkgdir(LineCableModels), LineCableModels._method_path(method))

    binding=Docs.Binding(LineCableModels, :domain)
    doc=first(values(Docs.meta(LineCableModels)[binding].docs))
    buffer=IOBuffer()
    @test DocStringExtensions.format(LineCableModels.METHODLIST, buffer, doc) === nothing
    rendered=String(take!(buffer))
    @test occursin("domain", rendered)
    @test occursin("src/interfaces.jl:", rendered)
    @test !occursin(pkgdir(LineCableModels), rendered)
end

@testitem "Core / owner-local numerics / transforms and earth kernels" tags=[:unit] begin
    using LinearAlgebra
    const Engine=LineCableModels.Engine
    const Transforms=Engine.Transforms

    matrix=[1.0 2.0 9.0; 4.0 5.0 8.0; 3.0 7.0 6.0]
    symmetrized=Transforms.reciprocity(matrix)
    @test symmetrized == transpose(symmetrized)
    @test diag(symmetrized) == diag(matrix)
    destination=copy(matrix)
    @test Transforms.reciprocity!(destination) === destination
    @test destination == symmetrized

    circulant=copy(matrix)
    @test Transforms.ideal_transposition!(circulant) === circulant
    @test all(circulant[i, j] == circulant[mod1(i + 1, 3), mod1(j + 1, 3)]
    for i in 1:3, j in 1:3)
    @test sum(circulant) ≈ sum(matrix)
    @test_throws DimensionMismatch Transforms.ideal_transposition!(ones(2, 3))
    @test Transforms.offdiagonal_ratio(Diagonal([1.0, 2.0])) == 0.0
    @test_throws DimensionMismatch Transforms.offdiagonal_ratio(zeros(2, 3))

    @test Engine.conductivity(Inf) == 0.0
    @test isinf(Engine.conductivity(0.0))
    @test Engine.conductivity(4.0) == 0.25
    @test Engine.bessel_difference(0.0 + 0.0im, 2.0, 4.0) ≈ log(2.0)
    @test isfinite(Engine.bessel_difference(1.0 + 1.0im, 2.0, 4.0))
end
