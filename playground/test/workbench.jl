const WorkbenchUI = LineCableModelsPlayground.WorkbenchUI
const TemplateWorkbench = LineCableModelsPlayground.TemplateWorkbench

@testset "workbench foundation" begin
    application = TemplateWorkbench.Application()
    state = WorkbenchUI.initialize(application, nothing)
    workbench = WorkbenchUI.compose(application, state)

    @test workbench isa WorkbenchUI.Workbench
    @test workbench.namespace == :template_workbench
    @test workbench.navigation.active === state.active_view
    @test workbench.workspace.active === state.active_view
    @test workbench.output.active === state.active_dock
    @test length(workbench.workspace.views) == 4
    @test length(workbench.output.tabs) == 4
    @test allunique(view.id for view in workbench.workspace.views)
    @test allunique(tab.id for tab in workbench.output.tabs)

    split = WorkbenchUI.SplitPane(
        :scene,
        :inputs;
        orientation=:horizontal,
        ratio=0.7
    )
    @test split.orientation == :horizontal
    @test split.ratio == 0.7
    @test split.resizable
    @test_throws ArgumentError WorkbenchUI.SplitPane(:a, :b; orientation=:diagonal)
    @test_throws ArgumentError WorkbenchUI.SplitPane(:a, :b; ratio=0.95)

    WorkbenchUI.handle!(
        application,
        state,
        TemplateWorkbench.SelectView(:results)
    )
    @test state.active_view[] == :results

    WorkbenchUI.handle!(
        application,
        state,
        TemplateWorkbench.SetLayer(:uncertainty)
    )
    @test state.layer[] == :uncertainty
    @test last(state.messages[]).text == "Display layer changed to uncertainty"

    WorkbenchUI.handle!(application, state, TemplateWorkbench.ResetView())
    @test state.layer[] == :nominal
    @test state.depth[] == 1.0
    @test state.spacing[] == 6.0
    @test state.scenario[] == "Baseline corridor"
    @test state.earth_model[] == :default
    @test !state.loading[]

    WorkbenchUI.handle!(application, state, TemplateWorkbench.ToggleLoading())
    @test state.loading[]
    @test state.active_dock[] == :jobs

    @test_throws ArgumentError WorkbenchUI.NavItem(
        Symbol("Not valid"),
        "Invalid"
    )
    @test_throws ArgumentError WorkbenchUI.Workbench(
        ;
        namespace=Symbol("not-valid"),
        identity=workbench.identity,
        navigation=workbench.navigation,
        workspace=workbench.workspace
    )

    @test haskey(
        Dict(LineCableModelsPlayground.WORKBENCH_ROUTES),
        "/workbenches/template"
    )
    @test TemplateWorkbench.app() isa LineCableModelsPlayground.Bonito.App
end
