using Bonito
using Observables

const Toolkit = LineCableModelsPlayground.Toolkit

@testset "reusable toolkit" begin
    name = TextInput(:scenario; value="Baseline", required=true)
    secret = SecretInput(:credential; required=true)
    notes = TextAreaInput(:notes; rows=3)
    radio = RadioGroup(:mode, [:planning => "Planning", :operations => "Operations"])
    segmented = SegmentedControl(:phase, [:a => "A", :b => "B"])
    combo = ComboBox(:earth_model, [:default => "Default", :carson => "Carson"])
    multiple = MultiSelect(:outputs, [:r => "R", :l => "L"]; selected=[:r])
    number = UnitNumberInput(:length; value=1, minimum=0, maximum=10, unit="km")
    range = RangeInput(:harmonics; lower=1, upper=50, minimum=1, maximum=100)

    @test name.value[] == "Baseline"
    @test !hasfield(SecretInput, :value)
    @test isempty(ComponentXRay.inspection(secret).bindings)
    @test occursin("absent", only(ComponentXRay.inspection(secret).notes))
    @test radio.value[] == "planning"
    @test segmented.value[] == "a"
    @test combo.value[] == "default"
    @test multiple.value[] == ["r"]
    @test number.value[] == 1.0
    @test range.lower[] == 1.0
    @test_throws ArgumentError TextInput(Symbol("Bad name"))
    @test_throws ArgumentError RadioGroup(:mode, String[])
    @test_throws ArgumentError UnitNumberInput(:value; value=11, maximum=10)
    @test_throws ArgumentError RangeInput(:range; lower=5, upper=2)

    form = Form(
        Field("Scenario", name),
        Field("Credential", secret),
        Field("Length", number);
        name=:study_form,
        validator=payload -> isempty(get(payload, "scenario", "")) ?
            Dict("scenario" => "Required") : nothing,
        cancel_label="Cancel"
    )
    @test form.status[] == :pristine
    @test !Toolkit.accept_form_payload!(form, Dict{String,Any}(
        "scenario" => "", "credential" => "not retained", "length" => "1"))
    @test form.status[] == :invalid
    @test form.errors[]["scenario"] == "Required"
    @test Toolkit.accept_form_payload!(form, Dict{String,Any}(
        "scenario" => "Study", "credential" => "not retained", "length" => "1"))
    @test submit_payload(form)["scenario"] == "Study"
    number.value[] = 4.0
    reset_form!(form)
    @test number.value[] == 1.0
    @test isnothing(submit_payload(form))

    callback_count = Ref(0)
    dialog = ConfirmDialog(:confirm, "Confirm", "Continue?";
        onconfirm=() -> (callback_count[] += 1))
    open_dialog!(dialog)
    @test dialog.open[]
    Toolkit.activate_dialog!(dialog, :confirm)
    @test callback_count[] == 1
    @test !dialog.open[]
    @test dialog.event[].action == :confirm
    @test_throws KeyError Toolkit.activate_dialog!(dialog, :missing)

    notice = InlineNotice(:warning, "Warning", "Check input"; dismissible=true)
    dismiss_notice!(notice)
    @test !notice.visible[]

    center = ToastCenter(; capacity=2)
    first_toast = push_toast!(center, :info, "One", "First"; timeout_ms=0)
    push_toast!(center, :success, "Two", "Second"; timeout_ms=0)
    third_toast = push_toast!(center, :warning, "Three", "Third"; timeout_ms=0)
    @test length(center.entries[]) == 2
    @test first_toast.id ∉ getfield.(center.entries[], :id)
    dismiss_toast!(center, third_toast.id)
    @test length(center.entries[]) == 1
    clear_toasts!(center)
    @test isempty(center.entries[])

    properties = PropertyGrid(
        PropertyItem("R", 0.1; unit="Ω/km"),
        PropertyItem("State", "ready"; state=:success);
        columns=2
    )
    table = DataTable(
        TableColumn(:name, "Name"),
        TableColumn(:value, "Value"; align=:right);
        rows=[(id="r", name="Resistance", value=0.1),
            (id="l", name="Inductance", value=0.4)],
        row_key=(row, _) -> row.id
    )
    @test length(Toolkit.table_rows(table)) == 2
    table.filter[] = "res"
    @test only(Toolkit.table_rows(table))[2].id == "r"
    table.filter[] = ""
    table.sort_by[] = :value
    table.sort_direction[] = :descending
    @test first(Toolkit.table_rows(table))[2].id == "l"

    viewport = ViewportFrame("Results", table)
    @test viewport.sizing == :viewport
    @test ViewportFrame("Inputs", form; sizing=:content).sizing == :content
    fitted = ViewportFrame("Fitted drawing", table; sizing=:fill, max_height="32rem")
    @test fitted.sizing == :fill && fitted.max_height == "32rem"
    @test viewport.max_height === nothing
    for length in ("100%", "480px", ".5em", "50vh", "50svh", "50dvh", "50vw")
        @test ViewportFrame("Capped", table; max_height=length).max_height == length
    end
    for length in ("", "none", "-1px", "0px", "NaNpx", "12", "20px; overflow:hidden", "calc(100% - 2rem)")
        @test_throws ArgumentError ViewportFrame("Invalid cap", table; max_height=length)
    end
    @test_throws ArgumentError ViewportFrame("Inputs", form; sizing=:unknown)
    set_viewport_state!(viewport, :loading; message="Preparing")
    @test viewport.state[] == :loading
    @test viewport.message[] == "Preparing"
    @test_throws ArgumentError set_viewport_state!(viewport, :unknown)
    disclosure = Disclosure("Metadata", properties; open=true)

    session = Bonito.Session(Bonito.NoConnection(); asset_server=Bonito.NoServer())
    status = StatusIndicator("Online"; tone=:success)
    action = ActionButton("Refresh status"; busy_label="Refreshing…")
    navigation = NavigationButton("Home"; href="/")
    @test navigation.href == "/" && navigation.target == "_top"
    @test NavigationButton("Dashboard"; href="/dashboard/?view=jobs#recent", target="_self").label == "Dashboard"
    for href in ("//other.example", "https://other.example", "javascript:alert(1)", "/\\other.example", "/\n/other.example")
        @test_throws ArgumentError NavigationButton("Unsafe"; href)
    end
    @test_throws ArgumentError NavigationButton(" "; href="/")
    @test_throws ArgumentError NavigationButton("Home"; href="/", target="popup")
    @test ComponentXRay.inspection(navigation).name == "NavigationButton"
    @test status.tone[] == :success
    @test !status.busy[]
    @test_throws ArgumentError StatusIndicator("Unknown"; tone=:invalid)
    @test !action.busy[] && action.clicks[] == 0
    @test ComponentXRay.inspection(status).name == "StatusIndicator"
    @test ComponentXRay.inspection(action).name == "ActionButton"
    for component in (name, secret, notes, radio, segmented, combo, multiple,
            number, range, form, dialog, notice, center, properties, table,
            viewport, fitted, disclosure, status, action, navigation)
        @test !isnothing(Bonito.jsrender(session, component))
    end
    split = WorkbenchUI.SplitPane(viewport, form; ratio=0.62)
    @test split.scroll == :panes
    canvas_split = WorkbenchUI.SplitPane(fitted, form; scroll=:parent)
    @test canvas_split.scroll == :parent
    @test !isnothing(Bonito.jsrender(session, canvas_split))
    @test_throws ArgumentError WorkbenchUI.SplitPane(viewport, form; scroll=:hidden)
    @test !isnothing(Bonito.jsrender(session, split))
    page = Toolkit.WorkspacePage("Engineering view", split; eyebrow="WORKSPACE", fill=true)
    @test page.fill
    @test !isnothing(Bonito.jsrender(session, page))
    @test ComponentXRay.inspection(page).name == "WorkspacePage"
    close(session)
end
