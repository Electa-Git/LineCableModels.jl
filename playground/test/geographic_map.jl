using Bonito
using EzXML
using JSON3
using UUIDs
using ZipFile

const BASIC_KML = """
<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
  <Document>
    <Placemark>
      <name>Route A</name>
      <LineString><coordinates>4.3,50.8,0 4.4,50.9,0</coordinates></LineString>
    </Placemark>
    <Placemark>
      <name>Known point</name>
      <Point><coordinates>4.35,50.85,0</coordinates></Point>
    </Placemark>
  </Document>
</kml>
"""

function write_kmz(path, entries)
    writer = ZipFile.Writer(path)
    try
        for (name, contents) in entries
            file = ZipFile.addfile(writer, name)
            write(file, contents)
        end
    finally
        close(writer)
    end
    return path
end

@testset "KML and KMZ validation" begin
    mktempdir() do directory
        kml_path = joinpath(directory, "route.kml")
        write(kml_path, BASIC_KML)
        @test isnothing(validate_kml_upload(kml_path))

        invalid_path = joinpath(directory, "invalid.kml")
        write(invalid_path, "<kml><Placemark><Point></kml>")
        @test_throws EzXML.XMLError validate_kml_upload(invalid_path)

        unsafe_path = joinpath(directory, "unsafe.kml")
        write(
            unsafe_path,
            "<!DOCTYPE kml [<!ENTITY x SYSTEM \"file:///etc/passwd\">]>" * BASIC_KML
        )
        @test_throws ArgumentError validate_kml_upload(unsafe_path)

        kmz_path = joinpath(directory, "route.kmz")
        write_kmz(kmz_path, ["doc.kml" => BASIC_KML])
        @test isnothing(validate_kml_upload(kmz_path))

        ambiguous_path = joinpath(directory, "ambiguous.kmz")
        write_kmz(
            ambiguous_path,
            ["first.kml" => BASIC_KML, "second.kml" => BASIC_KML]
        )
        @test_throws ArgumentError validate_kml_upload(ambiguous_path)

        traversal_path = joinpath(directory, "traversal.kmz")
        write_kmz(traversal_path, ["../doc.kml" => BASIC_KML])
        @test_throws ArgumentError validate_kml_upload(traversal_path)
    end
end

@testset "geographic map component contract" begin
    mktempdir() do directory
        target = UploadTarget(
            DirectoryUploadStore(directory),
            "route-source";
            extensions=[".kml", ".kmz"],
            validate=validate_kml_upload,
            retention=:session
        )
        upload = UploadField(target)
        component = GeographicMap(
            KMLUploadSource(upload);
            tags=["terminal" => "Terminal", "joint" => "Cable joint"],
            basemap=:none
        )
        @test isempty(component.references[])
        @test component.basemap == :none
        @test component.tags == ["terminal" => "Terminal", "joint" => "Cable joint"]

        uid = uuid4()
        component.reference_payload[] = JSON3.write([(
            uid=string(uid),
            source_feature="line-1",
            longitude=4.35,
            latitude=50.85,
            identifier="B4",
            tag="terminal",
        )])
        @test component.references[] == [MapReference(
            uid,
            "line-1",
            4.35,
            50.85,
            "B4",
            "terminal"
        )]

        @test_throws ArgumentError GeographicMap(
            KMLUploadSource(upload);
            tags=Pair{String,String}[]
        )
        @test_throws ArgumentError GeographicMap(
            KMLUploadSource(upload);
            tags=["terminal" => "Terminal", "terminal" => "Duplicate"]
        )
        @test_throws ArgumentError GeographicMap(
            KMLUploadSource(upload);
            tags=["terminal" => "Terminal"],
            basemap=:satellite
        )

        session = Bonito.Session(
            Bonito.NoConnection();
            asset_server=Bonito.NoServer()
        )
        @test !isnothing(Bonito.jsrender(session, component))
        close(session)
    end

    root = normpath(joinpath(@__DIR__, ".."))
    browser_source = read(
        joinpath(root, "assets", "vendor", "geographic-map.entry.js"),
        String
    )
    styles = read(joinpath(root, "assets", "geographic-map.css"), String)
    @test occursin("getClosestPoint", browser_source)
    @test occursin("extractStyles: false", browser_source)
    @test occursin("ResizeObserver", browser_source)
    @test occursin("lcm:upload-committed", browser_source)
    @test occursin("var\"data-map-role\"=\"editor\"", read(
        joinpath(root, "src", "geographic_map.jl"),
        String
    ))
    @test occursin("grid-template-columns", styles)
end
