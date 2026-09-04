@testset "visual contract coverage" begin
    root = normpath(joinpath(@__DIR__, ".."))

    function read_joined(paths)
        return join((read(path, String) for path in paths), '\n')
    end

    widget_sources = read_joined([
        joinpath(root, "src", "widgets.jl"),
        joinpath(root, "src", "console.jl"),
        joinpath(root, "src", "toolbar.jl"),
        joinpath(root, "src", "ribbon.jl"),
        joinpath(root, "src", "repeater.jl"),
        joinpath(root, "src", "uploads.jl"),
        joinpath(root, "src", "geographic_map.jl"),
        joinpath(root, "src", "power_system_canvas.jl"),
        joinpath(root, "src", "toolkit", "Forms.jl"),
        joinpath(root, "src", "toolkit", "Overlays.jl"),
        joinpath(root, "src", "toolkit", "DataViews.jl"),
        joinpath(root, "src", "widgets", "JobControls.jl"),
        joinpath(root, "assets", "geographic-map.css"),
        joinpath(root, "assets", "power-system-canvas.css"),
        joinpath(root, "assets", "repeater.css"),
        joinpath(root, "assets", "ribbon.css"),
        joinpath(root, "assets", "upload.css"),
        joinpath(root, "assets", "forms.css"),
        joinpath(root, "assets", "overlays.css"),
        joinpath(root, "assets", "data-views.css"),
    ])
    workbench_sources = read_joined([
        joinpath(root, "src", "workbench", "WorkbenchUI.jl"),
        joinpath(root, "src", "workbench", "workbench.css"),
        joinpath(root, "src", "workbenches", "TemplateWorkbench.jl"),
        joinpath(root, "src", "workbenches", "template_workbench.css"),
    ])

    # Every gallery specimen is accounted for explicitly. `:analogue` means a
    # workbench has the same visual primitive under its own structural class;
    # `:portable` means the gallery and workbench consume the same component
    # implementation, so a second CSS copy would itself be a defect.
    audit = [
        (route="/widgets/slider", gallery=".lc-native-range",
            parity=:analogue, workbench=".lc-wb-demo-range"),
        (route="/widgets/toggle", gallery=".lc-toggle-line",
            parity=:portable, workbench=nothing),
        (route="/widgets/dropdown", gallery=".lc-native-select",
            parity=:analogue, workbench=".lc-wb-choice-select"),
        (route="/widgets/text-input", gallery=".lc-native-field",
            parity=:analogue, workbench=".lc-wb-demo-fields input"),
        (route="/widgets/number-spinner", gallery=".lc-native-number",
            parity=:analogue, workbench=".lc-wb-demo-fields input"),
        (route="/widgets/actions", gallery=".lc-widget-actions",
            parity=:analogue, workbench=".lc-wb-command"),
        (route="/widgets/progress", gallery=".lc-progress",
            parity=:analogue, workbench=".lc-wb-demo-job progress"),
        (route="/widgets/console", gallery=".lc-console",
            parity=:analogue, workbench=".lc-wb-demo-output"),
        (route="/widgets/tabs", gallery=".lc-tabs-demo",
            parity=:analogue, workbench=".lc-wb-dock-tab"),
        (route="/widgets/toolbar", gallery=".lc-toolbar",
            parity=:analogue, workbench=".lc-wb-toolbar"),
        (route="/widgets/ribbon", gallery=".lc-ribbon",
            parity=:portable, workbench=nothing),
        (route="/widgets/control-panel", gallery=".lc-summary-grid",
            parity=:analogue, workbench=".lc-wb-demo-input-panel"),
        (route="/widgets/form-toolkit", gallery=".lc-form",
            parity=:portable, workbench=nothing),
        (route="/widgets/overlay-toolkit", gallery=".lc-dialog",
            parity=:portable, workbench=nothing),
        (route="/widgets/data-view-toolkit", gallery=".lc-viewport-frame",
            parity=:portable, workbench=nothing),
        (route="/widgets/repeater", gallery=".lc-repeater",
            parity=:portable, workbench=nothing),
        (route="/widgets/job-panel", gallery=".lc-job-panel",
            parity=:analogue, workbench=".lc-wb-demo-job"),
        (route="/widgets/file-upload", gallery=".lc-upload-field",
            parity=:portable, workbench=nothing),
        (route="/widgets/geographic-map", gallery=".lc-map-component",
            parity=:portable, workbench=nothing),
        (route="/widgets/power-system-canvas", gallery=".lc-power-system-canvas",
            parity=:portable, workbench=nothing),
    ]

    registered_routes = Set(first.(LineCableModelsPlayground.WIDGET_ROUTES))
    push!(registered_routes, "/widgets/job-panel")
    @test Set(spec.route for spec in audit) == registered_routes
    @test allunique(spec.route for spec in audit)
    @test all(spec.parity in (:analogue, :portable) for spec in audit)

    gallery_document = read(joinpath(root, "widgets", "index.qmd"), String)
    documented_routes = Set(
        match.captures[1]
        for match in eachmatch(r"route=\"(/widgets/[^\"]+)\"", gallery_document)
    )
    @test documented_routes == registered_routes

    for specimen in audit
        @test occursin(specimen.gallery, widget_sources)
        if specimen.parity == :analogue
            @test !isnothing(specimen.workbench)
            @test occursin(specimen.workbench, workbench_sources)
        else
            @test isnothing(specimen.workbench)
        end
    end

    contract = read(joinpath(root, "assets", "control-contract.css"), String)
    brand = read(joinpath(root, "assets", "brand.css"), String)
    quarto = read(joinpath(root, "_quarto.yml"), String)
    theme_selector = read(joinpath(root, "assets", "theme-selector.html"), String)
    workbench_ui = read(
        joinpath(root, "src", "workbench", "WorkbenchUI.jl"),
        String
    )
    @test occursin("assets/control-contract.css", quarto)
    @test occursin("DOM.style(CONTROL_CONTRACT)", widget_sources)
    @test occursin("DOM.style(CONTROL_CONTRACT)", workbench_ui)
    @test occursin(".lc-control-select option", contract)
    @test occursin(".lc-control-select option:checked", contract)
    @test occursin("lc-control-select lc-publisher-theme-select", theme_selector)

    # Palette values belong exclusively to brand.css. Other stylesheets may
    # consume these tokens but may not silently redefine their meaning.
    palette_tokens = (
        "lc-text",
        "lc-heading",
        "lc-strong-text",
        "lc-focus",
        "lc-hover-bg",
        "lc-active-bg",
        "lc-control-bg",
        "lc-option-bg",
        "lc-option-selected-bg",
        "lc-code-text",
        "lc-code-keyword",
        "lc-code-string",
        "lc-color-scheme",
    )
    for token in palette_tokens
        @test occursin(Regex("--$(token)\\s*:"), brand)
    end

    owned_styles = [
        joinpath(root, "assets", "theme.scss"),
        joinpath(root, "assets", "geographic-map.css"),
        joinpath(root, "assets", "power-system-canvas.css"),
        joinpath(root, "assets", "repeater.css"),
        joinpath(root, "assets", "ribbon.css"),
        joinpath(root, "assets", "upload.css"),
        joinpath(root, "assets", "forms.css"),
        joinpath(root, "assets", "overlays.css"),
        joinpath(root, "assets", "data-views.css"),
        joinpath(root, "src", "diagnostics", "component_xray.css"),
        joinpath(root, "src", "workbench", "workbench.css"),
        joinpath(root, "src", "workbenches", "template_workbench.css"),
        joinpath(root, "templates", "components", "collapsible-sidebar.css"),
    ]
    for path in owned_styles
        source = read(path, String)
        for token in palette_tokens
            @test !occursin(Regex("--$(token)\\s*:"), source)
        end
        @test !occursin(r"(?m)^\s*[^/\n]*\boption(?::checked)?[^\{\n]*\{", source)
    end

    # Literal color values are implementation-local palettes in disguise.
    # All first-party CSS consumes brand.css tokens; only vendored styles and
    # brand.css itself are allowed to own color literals.
    first_party_css = String[]
    for subtree in ("assets", "src", "templates")
        for (directory, directories, files) in walkdir(joinpath(root, subtree))
            filter!(name -> name != "vendor", directories)
            for file in files
                endswith(file, ".css") || continue
                path = joinpath(directory, file)
                path == joinpath(root, "assets", "brand.css") && continue
                push!(first_party_css, path)
            end
        end
    end
    @test !isempty(first_party_css)
    for path in first_party_css
        source = read(path, String)
        @test !occursin(r"(?i)(#[0-9a-f]{3,8}\b|rgb\()", source)
    end

    theme = read(joinpath(root, "assets", "theme.scss"), String)
    collapsible = read(
        joinpath(root, "templates", "components", "collapsible-sidebar.css"),
        String
    )
    @test occursin("code.sourceCode span.kw", theme)
    @test occursin(".quarto-title-block .code-tools-button:hover", theme)
    @test occursin("background: var(--lc-cs-active)", collapsible)
    @test occursin("background: var(--lc-cs-hover)", collapsible)
    @test occursin("background: var(--lc-scene-bg)", collapsible)

    # Every first-party native select has to opt into the semantic contract.
    # Looking forward from each constructor is deliberately simple and strict:
    # new constructors fail this test until their contract enrollment is made
    # explicit in the component source.
    julia_sources = String[]
    for (directory, _, files) in walkdir(joinpath(root, "src"))
        for file in files
            endswith(file, ".jl") || continue
            push!(julia_sources, read(joinpath(directory, file), String))
        end
    end
    all_julia = join(julia_sources, '\n')
    for found in eachmatch(r"DOM\.select\(", all_julia)
        stop = min(lastindex(all_julia), found.offset + 1200)
        @test occursin("lc-control-select", all_julia[found.offset:stop])
    end
    for found in eachmatch(r"(?<![A-Za-z])Dropdown\(", all_julia)
        stop = min(lastindex(all_julia), found.offset + 600)
        @test occursin("lc-control-select", all_julia[found.offset:stop])
    end
end
