@testitem "Materials / Material / numeric normalization and conversion" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport
] begin
    ordinary=Material(Float32(1.7e-8), big"2.3", 1, 20, Float32(0.004))
    @test ordinary isa Material{Float64}
    @test eltype(ordinary) === Float64

    converted=convert(Material{Float32}, ordinary)
    @test converted isa Material{Float32}
    @test converted.rho ≈ Float32(ordinary.rho)

    uncertain=Material(
        measurement(1.7e-8, 0.1e-8),
        2.3,
        1.0,
        20.0,
        0.004
    )
    @test eltype(uncertain) <: Measurement
    @test uncertainty(uncertain.rho) > 0
end

@testitem "Materials / MaterialsLibrary / dictionary and presentation contracts" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport
] begin
    empty_library=MaterialsLibrary(add_defaults = false)
    @test isempty(empty_library)
    @test sprint(show, MIME("text/plain"), empty_library) ==
          "MaterialsLibrary with 0 materials"

    copper=Material(1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
    @test add!(empty_library, :copper, copper) === empty_library
    @test empty_library["copper"] === copper
    @test collect(keys(empty_library)) == ["copper"]
    @test_throws ErrorException add!(empty_library, "copper", copper)

    fallback=Material(Inf, 1.0, 1.0, 20.0, 0.0)
    @test get(empty_library, "missing", fallback) === fallback
    @test_logs (:warn, r"not found") get(empty_library, "missing")

    @test_logs (:info, r"removed") delete!(empty_library, "copper")
    @test isempty(empty_library)
    @test_logs (:error, r"not found") @test_throws KeyError delete!(empty_library, "missing")

    defaults=@test_logs (:info, r"Initializing default materials") MaterialsLibrary()
    @test all(name -> haskey(defaults, name), ("air", "copper", "xlpe", "steel"))
    table=DataFrame(defaults)
    @test nrow(table) == length(defaults)
    @test Set(names(table)) == Set(["name", "rho", "eps_r", "mu_r", "T0", "alpha"])
    @test occursin("rho=", sprint(show, MIME("text/plain"), defaults["copper"]))
end

@testitem "Materials / presentation / bounded collection summaries" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport
] begin
    library=MaterialsLibrary(add_defaults = false)
    for index in 1:7
        add!(library, "material-$index", Material(index, 1.0, 1.0, 20.0, 0.0))
    end
    @test eltype(typeof(first(values(library)))) === Float64
    @test length(collect(values(library))) == 7

    summary=sprint(show, MIME"text/plain"(), library)
    @test occursin("MaterialsLibrary with 7 materials", summary)
    @test occursin("material-1", summary)
    @test occursin("... and 2 more", summary)

    reverse_library=MaterialsLibrary(add_defaults = false)
    for index in 7:-1:1
        add!(reverse_library, "material-$index", Material(index, 1.0, 1.0, 20.0, 0.0))
    end
    @test sprint(show, MIME"text/plain"(), reverse_library) == summary

    dictionary_summary=sprint(show, MIME"text/plain"(), library.data)
    @test occursin("Dict{String, Material} with 7 materials", dictionary_summary)
    @test occursin("... and 2 more", dictionary_summary)
    @test sprint(show, MIME"text/plain"(), reverse_library.data) == dictionary_summary
    singleton_summary=sprint(show, MIME"text/plain"(),
        Dict{String, Material}(
            "only"=>Material(1.0, 1.0, 1.0, 20.0, 0.0),
        ))
    @test occursin("with 1 material", singleton_summary)
    @test sprint(show, MIME"text/plain"(), Dict{String, Material}()) ==
          "Dict{String, Material} with 0 materials"
end
