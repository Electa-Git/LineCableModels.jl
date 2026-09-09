@testset "presentation foundation" begin
    root = normpath(joinpath(@__DIR__, ".."))
    source = joinpath(root, "presentations", "specimen.qmd")
    output = LineCableModelsPlayground.presentation_output(source)

    @test LineCableModelsPlayground.presentation_source() == source
    @test LineCableModelsPlayground.presentation_source(source) == source
    @test output == joinpath(root, "_site", "presentations", "specimen.html")
    @test LineCableModelsPlayground.presentation_probe_widget() isa
        LineCableModelsPlayground.Bonito.App

    document = read(source, String)
    @test occursin("lcm-deck-revealjs", document)
    @test occursin("lcm-layout-full-canvas", document)
    starter = read(joinpath(root, "presentations", "starter.qmd"), String)
    @test count("::: {.incremental}", starter) == 1
    @test count("::: {.incremental}", document) == 2
    @test occursin(r"(?m)^  - ", starter)
    @test occursin(r"(?m)^  - ", document)
    @test occursin("1. Open the presentation.", document)
    @test occursin("public-url=", document)
    @test !occursin(r"(?im)^\s*#\s+", document)
    @test !occursin(r"(?is)<\s*(script|style|section)(?:\s|>)", document)
    @test !occursin(r"(?is)<\s*iframe(?:\s|>)", document)

    missing_output = tempname() * ".html"
    errors = LineCableModelsPlayground.presentation_contract_errors(source, missing_output)
    @test errors == ["rendered presentation is missing: $missing_output"]

    extension = read(joinpath(root, "_extensions", "lcm-deck", "_extension.yml"), String)
    styles = read(joinpath(root, "_extensions", "lcm-deck", "deck.scss"), String)
    controller = read(joinpath(
        root, "_extensions", "lcm-deck", "deck-controller.js"
    ), String)
    shortcode = read(joinpath(root, "_extensions", "bonito", "bonito.lua"), String)
    @test occursin("disable-layout: true", extension)
    @test occursin("transition: none", extension)
    @test occursin("preload-iframes: false", extension)
    @test !occursin("scale(", styles)
    @test occursin("query.has('receiver')", controller)
    @test occursin("query.has('lcm-print')", controller)
    @test occursin("lcm:viewport-settled", controller)
    @test occursin("overviewshown", controller)
    @test occursin("beforeprint", controller)
    @test !occursin("deck.layout()", controller)
    @test occursin(".reveal.overview .slides > section", styles)
    @test occursin("html.lcm-print-view", styles)
    @test occursin("theme-init.html", extension)
    @test occursin("code-theme.css", extension)
    @test occursin("stopImmediatePropagation", controller)
    @test occursin("lcm-deck-status", controller)
    @test occursin("link.dataset.action = 'home'", controller)
    @test occursin("status.append(playgroundHomeLink())", controller)
    @test !occursin("playgroundHomeLink(true)", controller)
    @test occursin("?.closest('li')?.remove()", controller)
    math_notes = read(joinpath(root, "_extensions", "lcm-deck", "math-notes.js"), String)
    math_styles = read(joinpath(root, "_extensions", "lcm-deck", "math-notes.css"), String)
    authoring_filter = read(joinpath(root, "_extensions", "lcm-deck", "deck.lua"), String)
    @test occursin("registerPlugin('lcm-math-notes'", controller)
    @test occursin("math-notes.js", authoring_filter)
    @test occursin("math-notes.css", authoring_filter)
    @test occursin("data-lcm-math-target", authoring_filter)
    @test occursin("missing \\\\cssId anchor", authoring_filter)
    @test occursin("MessageHook('End Math', bind)", math_notes)
    @test occursin("hub.signal.RemoveHook(hook)", math_notes)
    @test occursin("window.Popper.createPopper", math_notes)
    @test occursin("beforeprint", math_notes)
    @test occursin("overviewshown", math_notes)
    @test occursin("abort.abort(); observer.disconnect()", math_notes)
    @test !occursin("contentDocument", math_notes)
    @test !occursin("deck.next", math_notes)
    @test occursin("var(--lc-panel-bg)", math_styles)
    @test occursin("var(--lc-text)", math_styles)
    @test !occursin("caret-color:", math_styles) # Inherit published-text.css, never duplicate its policy.
    @test occursin("user-select: text", math_styles)
    @test occursin("\\cssId{impedance-imag}", document)
    @test occursin(".lcm-math-note target=\"impedance-imag\"", document)
    initializer = read(joinpath(root, "_extensions", "lcm-deck", "deck-init.js"), String)
    @test occursin("url.searchParams.get('view') === 'print'", initializer)
    @test occursin("url.searchParams.delete('print-pdf')", initializer)
    @test occursin("data-lcm-src", shortcode)
    @test occursin("lcm-live-placeholder", shortcode)

    gallery = read(joinpath(root, "dev", "presentations.qmd"), String)
    recipes = read(joinpath(root, "presentations", "layouts.qmd"), String)
    snippets = collect(eachmatch(r"```\{\.markdown shortcodes=false\}\n(.*?)```"s, recipes))
    @test length(snippets) == length(LineCableModelsPlayground.PRESENTATION_LAYOUTS)
    for (layout, slots) in LineCableModelsPlayground.PRESENTATION_LAYOUTS
        @test occursin("layouts.html#$layout", gallery)
        @test occursin("{#$layout}", recipes)
        matching = filter(snippets) do snippet
            occursin("{.lcm-layout-$layout}", snippet.captures[1])
        end
        @test length(matching) == 1
        if length(matching) == 1
            @test count(r"::: \{\.lcm-slot(?:\s|\})", only(matching).captures[1]) == slots
        end
    end
    @test count("[Preview layout", gallery) == length(snippets)
    brand = read(joinpath(root, "assets", "brand.css"), String)
    @test count(r"--lc-presentation-laser\s*:", brand) == 1
    @test occursin("--lc-presentation-laser: #ff0000", brand)

    export_options = LineCableModelsPlayground.parse_presentation_export_options([
        "--pdf", "--output", "deck.pdf", "--quiet"
    ])
    @test export_options.output == "deck.pdf"
    @test export_options.quiet
    @test_throws ArgumentError LineCableModelsPlayground.parse_presentation_export_options(
        String[]
    )
end
