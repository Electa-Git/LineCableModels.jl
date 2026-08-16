@testitem "Commons / domain tags / fallback and identity" tags = [:unit] begin
    const Commons = LineCableModels.Commons

    @test Commons.domain(Int) === nothing
    @test Commons.domain(1) === nothing
    @test Commons.PhaseDomain() isa Commons.LineParamsDomain
    @test Commons.ModalDomain() isa Commons.LineParamsDomain
    @test Commons.PhaseDomain !== Commons.ModalDomain
end

@testitem "Commons / retired APIs / subsystem tombstones remain actionable" tags = [:unit] begin
    const DM = LineCableModels.DataModel
    const Engine = LineCableModels.Engine

    retired_calls = (
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
        exception = try
            invoke()
            nothing
        catch caught
            caught
        end
        @test exception isa ArgumentError
        message = sprint(showerror, exception)
        @test occursin("was removed from LineCableModels", message)
        @test occursin("legacy/fem-sector", message)
    end
end

@testitem "Commons / retired FEM API / actionable tombstone" tags = [:unit] begin
    const Commons = LineCableModels.Commons

    exception = try
        Commons.retired_fem_sector("SectorGeometry")
        nothing
    catch caught
        caught
    end
    @test exception isa ArgumentError
    message = sprint(showerror, exception)
    @test occursin("SectorGeometry was removed", message)
    @test occursin(Commons.LEGACY_FEM_SECTOR_BRANCH, message)
    @test occursin(Commons.LEGACY_FEM_SECTOR_COMMIT, message)
    @test occursin("Pkg.add", message)
end

@testitem "Commons / docstrings / sanitized method-list provenance" tags = [:unit] begin
    using DocStringExtensions

    const Commons = LineCableModels.Commons
    binding = Docs.Binding(Commons, :domain)
    doc = first(values(Docs.meta(Commons)[binding].docs))
    buffer = IOBuffer()

    @test DocStringExtensions.format(Commons.METHODLIST, buffer, doc) === nothing
    rendered = String(take!(buffer))
    @test occursin("domain", rendered)
    @test occursin("Commons.jl:", rendered)
    @test !occursin(pkgdir(LineCableModels), rendered)
end
