using Bonito
using Observables

@testset "application ribbon composition" begin
    received = ToolbarEvent[]
    binding = ToolbarBinding(event -> push!(received, event); namespace=:test_ribbon)

    normal = ToolbarButton(:open_project; icon=:open, label="Open")
    active = ToolbarButton(:place_bus; icon=:busbar, label="Bus", active=true)
    busy = ToolbarButton(:compile_model; icon=:chart, label="Compile", busy=true)
    disabled = ToolbarButton(
        :remote_run;
        icon=:play,
        label="Remote",
        disabled=true
    )
    toggle = ToolbarToggle(:snap_to_grid; icon=:grid, label="Snap", checked=true)
    number = ToolbarNumber(
        :line_length;
        label="Length",
        value=12,
        minimum=0.1,
        maximum=500,
        step=0.1,
        unit="km"
    )
    choice = ToolbarDropdown(
        :display_layer,
        [:nominal => "Nominal", :results => "Results"];
        icon=:layers,
        label="Layer"
    )

    @test active.active
    @test busy.busy
    @test disabled.disabled
    @test toggle.checked
    @test number.value == 12.0
    @test choice.selected == :nominal

    project = RibbonGroup("Project", normal, active; size=:medium)
    parameters = RibbonGroup("Parameters", toggle, number, choice; size=:small)
    home = RibbonTab(:home, "Home", project, parameters)
    remote = RibbonTab(
        :remote,
        "Remote",
        RibbonGroup("Execution", disabled);
        disabled=true
    )
    component = Ribbon(
        home,
        remote;
        binding,
        quick_access=(ToolbarButton(:quick_save; icon=:save),)
    )
    @test component.events[] === nothing

    session = Bonito.Session(
        Bonito.NoConnection();
        asset_server=Bonito.NoServer()
    )
    rendered = Bonito.jsrender(session, component)
    @test !isnothing(rendered)
    close(session)

    binding(ToolbarEvent(:test_ribbon, :open_project, nothing))
    binding(ToolbarEvent(:test_ribbon, :snap_to_grid, true))
    binding(ToolbarEvent(:test_ribbon, :line_length, 15.5))
    @test getfield.(received, :value) == [nothing, true, 15.5]

    @test_throws ArgumentError RibbonGroup("Empty")
    @test_throws ArgumentError RibbonGroup(
        "Duplicate",
        normal,
        ToolbarButton(:open_project; icon=:open)
    )
    @test_throws ArgumentError RibbonTab(:empty, "Empty")
    @test_throws ArgumentError Ribbon(; binding)
    @test_throws ArgumentError Ribbon(home; binding, active=:missing)
    @test_throws ArgumentError Ribbon(remote; binding)
    @test_throws ArgumentError ToolbarNumber(
        :invalid;
        label="Invalid",
        value=2,
        minimum=0,
        maximum=1
    )

    root = normpath(joinpath(@__DIR__, ".."))
    source = read(joinpath(root, "src", "ribbon.jl"), String)
    styles = read(joinpath(root, "assets", "ribbon.css"), String)
    @test occursin("ResizeObserver", source)
    @test occursin("closeMenus", source)
    @test occursin("role=\"tablist\"", source)
    @test occursin("data-collapsed", styles)
    @test occursin("lc-ribbon-overflow-menu", styles)
end
