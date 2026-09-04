@testitem "ImportExport / v1 scalar and Grid encoding" tags=[:unit] setup=[
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
        decoded=IE.deserialize_value(IE.serialize_value(original))
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

    deterministic=Grid((:conductor, :semicon))
    relative=Grid((1.0, 2.0), (1.0, 2.0))
    absolute=Grid((10.0,), AbsoluteError((0.5,)))
    for grid in (deterministic, relative, absolute)
        restored=IE.deserialize_value(IE.serialize_value(grid))
        @test typeof(restored) === typeof(grid)
        @test collect(restored) == collect(grid)
    end

    material_record=IE._material_record(Material(kind = :conductor, rho = 1.7e-8))
    material_record["kind"]=Dict("grid"=>["conductor", "semicon"])
    kinds=IE._decode_material_record(material_record)
    @test kinds isa Gridspace{Material}
    @test getproperty.(collect(kinds), :kind) == [:conductor, :semicon]

    @test_throws ArgumentError IE.deserialize_value(Dict("__type__"=>"Complex"))
    @test_throws ArgumentError IE.deserialize_value(Dict(
        "__type__"=>"Float64", "special"=>"huge"
    ))
    @test_throws ArgumentError IE.deserialize_value(Dict(
        "type"=>"FutureObject", "value"=>2
    ))
end

@testitem "ImportExport / v1 declaration round trip excludes derived state" tags=[:unit] setup=[
    ImportExportTestSupport,
    UseImportExportSupport,
    TestFixtures
] begin
    import LineCableModels.ImportExport as IE

    design=TestFixtures.mv_cable_design()
    encoded=IE.serialize_value(design)
    @test encoded["kind"] == "cable_design"
    @test encoded["cable_id"] == design.cable_id
    @test !haskey(encoded, "reference_frequency")

    function contains_key(value, key)
        value isa AbstractDict&&(
            haskey(value, key)||any(item->contains_key(item, key), values(value))
        )&&return true
        value isa AbstractVector&&return any(item->contains_key(item, key), value)
        return false
    end
    for derived in (
        "geometry", "terminal_order", "terminal_map", "effective",
        "connection_order", "engine", "grids"
    )
        @test !contains_key(encoded, derived)
    end

    restored=IE.deserialize_value(encoded)
    @test restored isa CableDesign
    @test restored !== design
    @test IE.serialize_value(restored) == encoded
    @test restored.terminal_order == design.terminal_order
    @test getproperty.(restored.geometry.regions, :primitive) ==
          getproperty.(design.geometry.regions, :primitive)
    @test getproperty.(restored.geometry.regions, :terminal) ==
          getproperty.(design.geometry.regions, :terminal)
    @test restored.terminal_map == design.terminal_map

    malformed=deepcopy(encoded)
    empty!(malformed["root"]["items"])
    @test_throws ArgumentError IE.deserialize_value(malformed)

    earth=EarthModel(100.0, 10.0, 1.0)
    connections=Dict(
        terminal=>(index==1 ? 1 : 0)
    for (index, terminal) in enumerate(design.terminal_order)
    )
    system=build(
        LineCableSystem,
        design,
        Pose2(0.0, -1.0);
        connections,
        environment = earth,
        system_id = "round-trip",
        line_length = 1000.0
    )
    encoded_system=IE.serialize_value(system)
    for derived in ("geometry", "terminal_order", "terminal_map", "connection_order")
        @test !contains_key(encoded_system, derived)
    end
    restored_system=IE.deserialize_value(encoded_system)
    @test IE.serialize_value(restored_system) == encoded_system
    @test restored_system.terminal_order == system.terminal_order
    @test restored_system.connection_order == system.connection_order
    @test restored_system.environment.vertical_layers == earth.vertical_layers
    @test restored_system.environment.layers == earth.layers

    problem=LineParametersProblem(
        system;
        temperature = 35.0,
        earth_props = earth,
        frequencies = [0.1, 50.0, 1.0e6],
        Γ = ComplexF64[0.0, 1.0e-5im, 2.0e-4 + 3.0e-4im]
    )
    mktempdir() do directory
        path=joinpath(directory, "problem.json")
        @test export_data(:json, problem; file_name = path) == abspath(path)
        restored_problem=import_data(
            :json,
            LineParametersProblem;
            file_name = path
        )
        @test IE.serialize_value(restored_problem) == IE.serialize_value(problem)
        @test restored_problem.system.terminal_order == system.terminal_order
        @test restored_problem.system.connection_order == system.connection_order
    end

    parametric_system_record=deepcopy(encoded_system)
    parametric_system_record["positions"]=IE.serialize_value(Grid((
        Pose2(-0.1, -1.0),
        Pose2(0.1, -1.0)
    )))
    system_space=IE.deserialize_value(parametric_system_record)
    @test system_space isa Gridspace{LineCableSystem}
    @test eltype(system_space) === Any
    @test Base.IteratorEltype(typeof(system_space)) isa Base.EltypeUnknown
    @test getproperty.(collect(system_space), :positions) == [
        [Pose2(-0.1, -1.0)],
        [Pose2(0.1, -1.0)]
    ]

    copper=Material(kind = :conductor, rho = 1.7241e-8)
    compacted=build(
        CableDesign,
        "compacted-round-trip",
        terminal(
            :core,
            stranded(
                copper;
                shape = Disk(0.35e-3),
                layers = 2,
                n = (6, 12),
                compact = FillFactor(1),
                boundary = Disk(sqrt(19)*0.35e-3)
            )
        )
    )
    compacted_record=IE.serialize_value(compacted)
    restored_compacted=IE.deserialize_value(compacted_record)
    @test IE.serialize_value(restored_compacted) == compacted_record
    @test getproperty.(restored_compacted.geometry.regions, :terminal) ==
          getproperty.(compacted.geometry.regions, :terminal)
    @test area.(getproperty.(restored_compacted.geometry.regions, :primitive)) ==
          area.(getproperty.(compacted.geometry.regions, :primitive))
    @test all(
        region -> region.source.primitive isa Disk &&
                  region.primitive isa LineCableModels.DataModel.Polygon,
        restored_compacted.geometry.regions
    )

    packed=build(
        CableDesign,
        "bounded-round-trip",
        terminal(
            :core,
            stranded(
                copper;
                shape = Disk(0.35e-3),
                layers = 2,
                boundary = Disk(2e-3)
            )
        )
    )
    packed_record=IE.serialize_value(packed)
    @test packed_record["root"]["item"]["items"][1]["items"][1]["boundary"]["kind"] ==
          "disk"
    restored_packed=IE.deserialize_value(packed_record)
    @test restored_packed == packed
    @test IE.serialize_value(restored_packed) == packed_record

    scheduled=build(
        CableDesign,
        "scheduled-round-trip",
        terminal(
            :core,
            stranded(
                copper;
                shape = Disk(0.35e-3),
                layers = 2,
                n = (6, 12),
                lay = (LayRatio(12), Pitch(0.1)),
                dir = (1, -1),
                φ0 = (0.0, 0.2),
                boundary = Disk(2e-3)
            )
        )
    )
    scheduled_record=IE.serialize_value(scheduled)
    restored_scheduled=IE.deserialize_value(scheduled_record)
    @test restored_scheduled == scheduled
    @test IE.serialize_value(restored_scheduled) == scheduled_record

    sector_boundary=Sector(
        span = deg2rad(119.0),
        r_base = 1.10e-3,
        r_back = 10.24e-3,
        fillet = 1.02e-3
    )
    bounded=build(
        CableDesign,
        "bounded-stranded-round-trip",
        terminal(
            :core,
            stranded(
                copper;
                shape = Disk(0.35e-3),
                layers = 1,
                n = 6,
                boundary = sector_boundary
            )
        )
    )
    bounded_record=IE.serialize_value(bounded)
    restored_bounded=IE.deserialize_value(bounded_record)
    @test restored_bounded == bounded
    @test IE.serialize_value(restored_bounded) == bounded_record
    bounded_group=only(only(bounded.root.item.items).items)
    @test bounded_group isa Group
    @test bounded_group.boundary == sector_boundary

    declarations=(
        Ring(capacity(); r = nothing, φ0 = 0.2, gap_frac = 0.03),
        Polar(nr = 2, nφ = 6, r0 = 0.0, dr = 2e-3),
        Fill(r = 10e-3, φ = pi/6),
        Lattice(nx = 2, ny = 3, dx = 10e-3, dy = 12e-3),
        Helix(LayRatio(11); dir = -1, φ0 = 0.1),
        FillFactor(0.9)
    )
    for declaration in declarations
        @test IE.deserialize_value(IE.serialize_value(declaration)) == declaration
    end

    left=terminal(:left, core(copper; r = 0.8e-3))
    right=terminal(:right, core(copper; r = 1.0e-3))
    explicit=assembly(at(left, -2e-3, 0), at(right, 2e-3, 0; φ = 0.1))
    explicit_record=IE.serialize_value(explicit)
    @test explicit_record["kind"] == "assembly"
    @test haskey(explicit_record, "members")
    @test IE.deserialize_value(explicit_record) == explicit

    repeated=cores(
        terminal(:core, core(copper; r = 0.5e-3));
        n = 3,
        r = 2e-3,
        names = (:a, :b, :c)
    )
    repeated_record=IE.serialize_value(repeated)
    @test repeated_record["kind"] == "assembly"
    @test haskey(repeated_record, "pattern")
    @test !haskey(repeated_record, "members")
    @test IE.deserialize_value(repeated_record) == repeated

    descriptive=build(
        CableDesign,
        "descriptive-data",
        left;
        nominal_data = (rated_voltage = 220.0, designation = "test")
    )
    restored_descriptive=IE.deserialize_value(IE.serialize_value(descriptive))
    @test restored_descriptive.nominal_data == descriptive.nominal_data

    library=CablesLibrary()
    catalogue_record=DatasheetInfo(designation_code = "MV cable", voltage_kv = 30.0)
    add!(library, design; catalogue = catalogue_record)
    mktempdir() do directory
        for extension in ("json", "jls")
            path=joinpath(directory, "library.$extension")
            save(library; file_name = path)
            restored_library=CablesLibrary()
            @test load!(restored_library; file_name = path) === restored_library
            @test only(keys(restored_library)) == design.cable_id
            @test IE.serialize_value(only(values(restored_library))) == encoded
            restored_catalogue=catalogue(restored_library, design.cable_id)
            @test restored_catalogue.designation_code == catalogue_record.designation_code
            @test restored_catalogue.voltage_kv == catalogue_record.voltage_kv
        end
    end
end

@testitem "ImportExport / Draft 2020-12 cable document" tags=[:unit] setup=[
    ImportExportTestSupport,
    UseImportExportSupport,
    TestFixtures
] begin
    using JSON3
    using JSONSchema
    import LineCableModels.ImportExport as IE

    materials=MaterialsLibrary(add_defaults = false)
    add!(materials, "copper", TestFixtures.copper_material())
    materials_document=IE._json_document(materials)
    @test materials_document["\$schema"] == IE.JSON_SCHEMA_DIALECT
    @test materials_document["format"] == IE.MATERIALS_SCHEMA
    @test materials_document["version"] == "1.0.0"

    cables=CablesLibrary()
    add!(
        cables,
        TestFixtures.mv_cable_design();
        catalogue = DatasheetInfo(designation_code = "NA2XS(FL)2Y", U = 30.0)
    )
    cables_document=IE._json_document(cables)
    @test cables_document["\$schema"] == IE.JSON_SCHEMA_DIALECT
    @test cables_document["format"] == IE.CABLES_SCHEMA
    @test cables_document["version"] == "1.0.0"
    @test cables_document["root"]["kind"] == "cable_library"
    @test only(values(cables_document["root"]["cables"]))["catalogue"] == Dict(
        "designation_code" => "NA2XS(FL)2Y",
        "U" => IE.serialize_value(30.0)
    )
    @test length(cables_document["materials"]) <
          length(first(values(cables)).geometry.regions)
    schema_path=joinpath(
        pkgdir(LineCableModels), "schemas", "cable-design.schema.json"
    )
    @test isfile(schema_path)
    schema=JSONSchema.Schema(JSON3.read(read(schema_path, String), Dict{String, Any}))
    @test JSONSchema.validate(schema, cables_document) === nothing
    @test JSONSchema.validate(schema, materials_document) === nothing

    packed_cables=CablesLibrary()
    add!(packed_cables,
        build(
            CableDesign,
            "bounded-schema",
            terminal(
                :core,
                stranded(
                    TestFixtures.copper_material();
                    shape = Disk(0.5e-3),
                    layers = 2,
                    boundary = Disk(3e-3)
                )
            )
        ))
    @test JSONSchema.validate(
        schema,
        IE._json_document(packed_cables)
    ) === nothing

    earth=EarthModel(100.0, 10.0, 1.0)
    design=TestFixtures.mv_cable_design()
    connections=Dict(
        terminal=>(index==1 ? 1 : 0)
    for (index, terminal) in enumerate(design.terminal_order)
    )
    system=build(
        LineCableSystem,
        design,
        Pose2(0.0, -1.0);
        connections,
        environment = earth
    )
    system_document=IE._json_document(system)
    @test JSONSchema.validate(schema, system_document) === nothing

    parametric_document=deepcopy(cables_document)
    raw_design=only(values(parametric_document["root"]["cables"]))
    function referenced_tags(node, material_name)
        node isa AbstractDict||return String[]
        tags=String[]
        get(node, "kind", nothing)=="region"&&
            get(node, "material", nothing)==material_name&&
            push!(tags, node["tag"])
        for child_name in ("items", "item", "fill", "wall", "root")
            child=get(node, child_name, nothing)
            child isa AbstractVector&&foreach(
                item->append!(tags, referenced_tags(item, material_name)), child
            )
            child isa AbstractDict&&append!(tags, referenced_tags(child, material_name))
        end
        return tags
    end
    material_names=collect(keys(parametric_document["materials"]))
    material_name=first(sort(
        material_names;
        by = name->length(referenced_tags(raw_design, name)),
        rev = true
    ))
    varied_tags=referenced_tags(raw_design, material_name)
    @test length(varied_tags) > 1
    parametric_document["materials"][material_name]["rho"]=IE.serialize_value(Grid((
        1.7e-8, 2.0e-8)))
    @test JSONSchema.validate(schema, parametric_document) === nothing
    decoded_materials=IE._document_materials(parametric_document)
    design_space=IE._decode_design(raw_design, decoded_materials)
    @test design_space isa Gridspace{CableDesign}
    @test length(design_space) == 2
    decoded_designs=collect(design_space)
    @test all(design->design isa CableDesign, decoded_designs)
    for design in decoded_designs
        selected_rhos=unique(
            region.source.material.rho
        for region in design.geometry.regions
        if String(region.source.tag) in varied_tags
        )
        @test length(selected_rhos) == 1
    end
    @test only(unique(
        region.source.material.rho
    for region in decoded_designs[1].geometry.regions
    if String(region.source.tag) in varied_tags
    )) != only(unique(
        region.source.material.rho
    for region in decoded_designs[2].geometry.regions
    if String(region.source.tag) in varied_tags
    ))

    invalid=deepcopy(cables_document)
    delete!(invalid["root"]["cables"]["test_cable"], "root")
    @test JSONSchema.validate(schema, invalid) !== nothing

    @test_throws ArgumentError IE._json_path("library.archive")
    @test endswith(IE._json_path("library"), "library.json")
end
