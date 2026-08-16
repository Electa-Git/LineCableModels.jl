@testitem "ImportExport / PSCAD / system structure and physical values" tags=[:integration] setup=[
    ImportExportTestSupport,
    UseImportExportSupport,
    TestFixtures
] begin
    using EzXML

    system=TestFixtures.three_phase_system()
    frequencies=[50.0, 500.0]
    earth=EarthModel(frequencies, 100.0, 10.0, 1.0)

    mktempdir() do directory
        requested=joinpath(directory, "system.pscx")
        output=export_data(:pscad, system, earth; file_name = requested)
        @test output == joinpath(directory, "$(system.system_id)_system.pscx")
        @test isfile(output)
        @test filesize(output) > 200

        document=readxml(output)
        project=root(document)
        @test nodename(project) == "project"
        @test project["name"] == system.system_id
        @test !isempty(findall("//User", document))
        @test occursin("LineCableModels", read(output, String))
    end
end
