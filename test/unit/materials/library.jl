@testitem "Materials / Material / numeric normalization and conversion" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport
] begin
    ordinary=Material(
        :conductor, Float32(1.7e-8), big"2.3", 1, 20, Float32(0.004);
        rho_thermal = Float32(0.25), theta_max = 95,
        tan_delta = Float32(0.01), sigma_solar = 0.5
    )
    @test ordinary isa Material{BigFloat}
    @test eltype(ordinary) === BigFloat
    @test ordinary.kind === :conductor

    converted=convert(Material{Float32}, ordinary)
    @test converted isa Material{Float32}
    @test converted.rho ≈ Float32(ordinary.rho)
    @test converted.kind === ordinary.kind

    defaulted=Material(:conductor, Float32(1.7e-8))
    @test defaulted isa Material{Float32}
    @test defaulted.kind === :conductor

    uncertain=Material(
        :conductor,
        measurement(1.7e-8, 0.1e-8),
        2.3,
        1.0,
        20.0,
        0.004
    )
    @test eltype(uncertain) <: Measurement
    @test uncertainty(uncertain.rho) > 0
    @test_throws ArgumentError Material(Symbol(""), 1.0)
end

@testitem "Materials / MaterialsLibrary / dictionary and presentation contracts" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport
] begin
    empty_library=MaterialsLibrary(add_defaults = false)
    @test empty_library isa AbstractDict{String, Material}
    @test keytype(typeof(empty_library)) === String
    @test valtype(typeof(empty_library)) === Material
    @test eltype(typeof(empty_library)) === Pair{String, Material}
    @test pairs(empty_library) === empty_library
    @test isempty(empty_library)
    @test sprint(show, MIME("text/plain"), empty_library) ==
          "MaterialsLibrary · 0 materials"

    copper=Material(:conductor, 1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
    @test add!(empty_library, :copper, copper) === empty_library
    @test empty_library["copper"] === copper
    @test collect(empty_library) isa Vector{Pair{String, Material}}
    @test collect(keys(empty_library)) == ["copper"]
    @test_throws ArgumentError add!(empty_library, "copper", copper)

    fallback=Material(:insulator, Inf, 1.0, 1.0, 20.0, 0.0)
    @test get(empty_library, "missing", fallback) === fallback
    @test get(() -> fallback, empty_library, "missing") === fallback
    @test get(() -> fallback, empty_library, "copper") === copper
    @test_throws MethodError get(empty_library, "missing")
    @test !haskey(empty_library, :copper)
    @test get(empty_library, :copper, fallback) === fallback
    @test_throws KeyError empty_library[:copper]

    replacement=Material(:conductor, 1.8e-8, 1.0, 1.0, 20.0, 0.004)
    @test setindex!(empty_library, replacement, "copper") === empty_library
    @test empty_library["copper"] === replacement
    @test (empty_library["copper"] = copper) === copper
    @test empty_library["copper"] === copper
    empty_library["copper"]=replacement

    copied=copy(empty_library)
    @test copied isa MaterialsLibrary
    @test copied !== empty_library
    @test copied.data !== empty_library.data
    @test copied["copper"] === replacement
    @test delete!(copied, "copper") === copied
    @test haskey(empty_library, "copper")

    blank=empty(empty_library)
    @test blank isa MaterialsLibrary
    @test isempty(blank)
    @test empty!(copied) === copied
    @test isempty(copied)

    @test delete!(empty_library, "copper") === empty_library
    @test isempty(empty_library)
    @test delete!(empty_library, "missing") === empty_library
    @test delete!(empty_library, :missing) === empty_library

    defaults=MaterialsLibrary()
    @test all(name -> haskey(defaults, name), ("air", "copper", "xlpe", "steel"))
    shown=sprint(show, MIME("text/plain"), defaults)
    @test contains(shown, "copper")
    @test contains(shown, "$(length(defaults)) materials")
    @test occursin("ρ", sprint(show, MIME("text/plain"), defaults["copper"]))
    @test occursin("nΩ·m", sprint(show, MIME("text/plain"), defaults["copper"]))
end

@testitem "Materials / presentation / bounded collection summaries" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport
] begin
    library=MaterialsLibrary(add_defaults = false)
    for index in 1:7
        add!(library, "material-$index", Material(:conductor, index, 1.0, 1.0, 20.0, 0.0))
    end
    @test eltype(typeof(first(values(library)))) === Float64
    @test length(collect(values(library))) == 7

    summary=sprint(show, MIME"text/plain"(), library)
    @test occursin("MaterialsLibrary · 7 materials", summary)
    @test occursin("material-1", summary)

    reverse_library=MaterialsLibrary(add_defaults = false)
    for index in 7:-1:1
        add!(
            reverse_library,
            "material-$index",
            Material(:conductor, index, 1.0, 1.0, 20.0, 0.0)
        )
    end
    @test sprint(show, MIME"text/plain"(), reverse_library) == summary

    narrow=sprint(show, MIME"text/plain"(), library;
        context = IOContext(stdout, :limit=>true, :displaysize=>(5, 80)))
    @test occursin("⋮", narrow)
    @test occursin("more materials", narrow)

    # MaterialsLibrary owns deterministic semantic display. Its backing Dict
    # retains ordinary Base display; the package does not pirate Dict methods.
    dictionary_summary=sprint(show, MIME"text/plain"(), library.data)
    @test occursin("Dict{String", dictionary_summary)
    @test parentmodule(which(show, (
        IO,
        MIME"text/plain",
        typeof(library.data)
    ))) === Base
end
