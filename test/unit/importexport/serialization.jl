@testitem "ImportExport / schema values / lossless supported tags" tags=[:unit] setup=[
    ImportExportTestSupport,
    UseImportExportSupport
] begin
    import LineCableModels.ImportExport as IE

    values=Any[
        nothing, true, 7,
        Float16(1.25), Float32(1.25), 1.25, BigFloat("1.25"),
        Inf, -Inf, NaN, 2.0 - 3.0im, "cable", :phase,
        [1, 2.5, -1.0im], (1.0, 2.0), Dict(:count=>2, :scale=>0.5)
    ]
    for original in values
        encoded=IE.serialize_value(original)
        decoded=IE.deserialize_value(encoded)
        if original isa AbstractFloat&&isnan(original)
            @test isnan(decoded)
        elseif original isa Tuple
            @test decoded == collect(original)
        elseif original isa AbstractDict
            @test decoded == Dict("count"=>2, "scale"=>0.5)
        else
            @test decoded == original
        end
    end

    @test_throws ArgumentError IE.deserialize_value(Dict("__type__"=>"Complex"))
    @test_throws ArgumentError IE.deserialize_value(Dict(
        "__type__"=>"Float64", "special"=>"huge"
    ))
    @test_throws ArgumentError IE.deserialize_value(Dict(
        "__type__"=>"FutureScalar", "value"=>2
    ))
    @test_throws ArgumentError IE.deserialize_value(Dict(
        "type"=>"FutureObject", "value"=>2
    ))
end

@testitem "ImportExport / schema objects / dispatched type tags round trip" tags=[:unit] setup=[
    ImportExportTestSupport,
    UseImportExportSupport,
    TestFixtures
] begin
    import LineCableModels.ImportExport as IE

    design=TestFixtures.mv_cable_design()
    encoded=IE.serialize_value(design)
    @test encoded["type"] == "CableDesign"
    @test encoded["cable_id"] == design.cable_id

    function contains_key(value, key)
        value isa AbstractDict&&(
            haskey(value, key)||any(item->contains_key(item, key), values(value))
        )&&return true
        value isa AbstractVector&&return any(item->contains_key(item, key), value)
        return false
    end
    @test !contains_key(encoded, "temperature")
    @test !contains_key(encoded, "__julia_type__")

    restored=IE.deserialize_value(encoded)
    @test restored isa CableDesign
    @test restored !== design
    @test IE.serialize_value(restored) == encoded

    material=Material(
        :conductor,
        Float32(1), Float32(2), Float32(3), Float32(20), Float32(0.01)
    )
    material_roundtrip=IE.deserialize_value(IE.serialize_value(material))
    @test material_roundtrip == material
    @test eltype(material_roundtrip) === Float32

    missing_radius=deepcopy(encoded)
    delete!(
        missing_radius["components"][1]["conductor_group"]["layers"][1],
        "r_ex"
    )
    @test_throws ArgumentError IE.deserialize_value(missing_radius)

    empty_group=deepcopy(encoded)
    empty!(empty_group["components"][1]["conductor_group"]["layers"])
    @test_throws ArgumentError IE.deserialize_value(empty_group)
end

@testitem "ImportExport / documents / explicit version and schema" tags=[:unit] setup=[
    ImportExportTestSupport,
    UseImportExportSupport,
    TestFixtures
] begin
    import LineCableModels.ImportExport as IE

    materials=MaterialsLibrary(add_defaults = false)
    add!(materials, "copper", TestFixtures.copper_material())
    materials_document=IE._json_document(materials)
    @test materials_document["schema"] == IE.MATERIALS_SCHEMA
    @test materials_document["version"] == IE.JSON_SCHEMA_VERSION
    @test haskey(materials_document, "materials")

    cables=CablesLibrary()
    add!(cables, TestFixtures.mv_cable_design())
    cables_document=IE._json_document(cables)
    @test cables_document["schema"] == IE.CABLES_SCHEMA
    @test cables_document["version"] == IE.JSON_SCHEMA_VERSION
    @test haskey(cables_document, "cables")

    @test_throws ArgumentError IE._json_path("library.archive")
    @test endswith(IE._json_path("library"), "library.json")
end
