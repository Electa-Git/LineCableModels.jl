# Ownership is checked on actual loaded bindings by quality/explicit_imports.jl.
# Lowering and reuse counts run through compute in integration/formulation_grid.jl
# and integration/line_parameters.jl. Source spelling, directory layout and
# absence of superseded helpers are not current architecture requirements.
@testitem "Core / architecture / earth models are immutable values" tags=[:unit] begin
    earth=@earth begin
        layer(rho = 100.0, thickness = 5.0)
        layer(rho = 500.0)
    end
    earth_layer=layer(rho = 50.0)

    @test !ismutabletype(typeof(earth))
    @test eltype(typeof(earth)) === Float64
    @test fieldtype(typeof(earth), :layers) ===
          NTuple{3, EarthLayer{Float64}}
    @test !hasmethod(add!, Tuple{typeof(earth), typeof(earth_layer)})
    @test_throws MethodError setindex!(earth.layers, earth_layer, 2)
end

@testitem "Core / selectors / owner-local symbol facades" tags=[:unit] begin
    const Engine=LineCableModels.Engine
    const ImportExport=LineCableModels.ImportExport
    const WirePatterns=LineCableModels.ParametricBuilder.WirePatterns

    @test which(Engine.Formulation, (Symbol,)).module === Engine
    @test which(ImportExport.export_data, (Symbol, Nothing)).module === ImportExport
    @test which(ImportExport.import_data, (Symbol, Nothing)).module === ImportExport
    @test which(Engine.LineParameters, (Symbol, String)).module === ImportExport
    @test which(
        getindex,
        (WirePatterns.WireEstimate{Float64, WirePatterns.HexaPattern{Float64}}, Symbol)
    ).module === WirePatterns

    @test Engine.Formulation() isa Engine.LineParametersFormulation
    @test_throws MethodError Engine.Formulation(:analytical)
    @test_throws MethodError Engine.Formulation(:line_cable_models)
    @test_throws MethodError ImportExport.export_data(:unregistered, nothing)
    @test_throws MethodError ImportExport.import_data(:unregistered, nothing)
end
