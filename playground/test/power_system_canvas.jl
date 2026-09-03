using Bonito
using JSON3

@testset "power-system canvas contract" begin
    diagram = LineCableModelsPlayground.power_system_specimen_diagram()
    component = PowerSystemCanvas(diagram)

    @test isnothing(snapshot(component))
    encoded = JSON3.write(diagram)
    @test validate_power_system_diagram(encoded) == encoded
    @test length(JSON3.read(encoded).elements) == 5
    @test length(JSON3.read(encoded).buses) == 3

    compiled = JSON3.write((
        schema="lcm.power-system-topology/1",
        diagram=diagram,
        topology=(
            nodes=[(id="node-1", terminals=["GRID.t_bottom", "B-HV"])],
            terminal_to_node=[(terminal="GRID.t_bottom", node="node-1")],
            diagnostics=[],
        ),
    ))
    component.bridge_payload[] = compiled
    accepted = snapshot(component)
    @test accepted["schema"] == "lcm.power-system-topology/1"
    @test accepted["diagram"]["meta"]["title"] == "Industrial feeder specimen"

    invalid_version = JSON3.write(merge(diagram, (version="2",)))
    @test_throws ArgumentError validate_power_system_diagram(invalid_version)
    @test_throws ArgumentError validate_power_system_diagram(JSON3.write((version="1",)))
    @test_throws ArgumentError validate_power_system_topology(JSON3.write((
        schema="unknown",
    )))

    session = Bonito.Session(
        Bonito.NoConnection();
        asset_server=Bonito.NoServer()
    )
    @test !isnothing(Bonito.jsrender(session, component))
    close(session)

    root = normpath(joinpath(@__DIR__, ".."))
    entry = read(
        joinpath(root, "assets", "vendor", "power-system-canvas.entry.jsx"),
        String
    )
    bundle = read(
        joinpath(root, "assets", "vendor", "power-system-canvas.bundle.js"),
        String
    )
    styles = read(joinpath(root, "assets", "power-system-canvas.css"), String)
    @test occursin("compile(diagram)", entry)
    @test occursin("useEditorStore.persist.setOptions", entry)
    @test occursin("MutationObserver", entry)
    @test occursin("LCMPowerSystemCanvas", bundle)
    @test occursin("--lc-panel-bg", styles)
end
