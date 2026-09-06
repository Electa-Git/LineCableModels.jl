const PRESENTATION_ROOT = joinpath(PLAYGROUND_ROOT, "presentations")
const DEFAULT_PRESENTATION = joinpath(PRESENTATION_ROOT, "specimen.qmd")
const PRESENTATION_LAYOUTS = Dict(
    "full-canvas" => 1,
    "balanced" => 2,
    "feature-sidebar" => 2,
    "top-split" => 3,
    "dashboard-grid" => 4,
    "media-story" => 2,
)

function presentation_probe_widget()
    return App(; title="Presentation canvas probe · LineCableModels playground") do session
        identity = string(uuid4())
        canvas = DOM.canvas(
            tabindex="0",
            role="application",
            var"aria-label"="Presentation coordinate probe"
        )
        status = DOM.output("Select a point inside the canvas")
        probe = DOM.div(
            canvas,
            status;
            class="lc-presentation-probe",
            var"data-session-id"=identity
        )
        root = widget_shell(
            "PRESENTATION PROBE",
            "Persistent canvas and exact pointer coordinates",
            probe
        )
        Bonito.onload(session, probe, js"""
        (element) => {
            const canvas = element.querySelector('canvas');
            const status = element.querySelector('output');
            const context = canvas.getContext('2d');
            let lastPoint = null;

            const draw = () => {
                const box = canvas.getBoundingClientRect();
                const ratio = window.devicePixelRatio || 1;
                canvas.width = Math.max(1, Math.round(box.width * ratio));
                canvas.height = Math.max(1, Math.round(box.height * ratio));
                context.setTransform(ratio, 0, 0, ratio, 0, 0);
                context.clearRect(0, 0, box.width, box.height);
                const style = getComputedStyle(document.documentElement);
                context.strokeStyle = style.getPropertyValue('--lc-border-soft').trim();
                context.lineWidth = 1;
                for (let x = 24; x < box.width; x += 24) {
                    context.beginPath(); context.moveTo(x, 0); context.lineTo(x, box.height); context.stroke();
                }
                for (let y = 24; y < box.height; y += 24) {
                    context.beginPath(); context.moveTo(0, y); context.lineTo(box.width, y); context.stroke();
                }
                context.strokeStyle = style.getPropertyValue('--lc-focus').trim();
                context.strokeRect(0.5, 0.5, Math.max(0, box.width - 1), Math.max(0, box.height - 1));
                if (lastPoint) {
                    context.fillStyle = style.getPropertyValue('--lc-link').trim();
                    context.beginPath();
                    context.arc(lastPoint.x, lastPoint.y, 6, 0, Math.PI * 2);
                    context.fill();
                }
                element.dataset.backingWidth = String(canvas.width);
                element.dataset.backingHeight = String(canvas.height);
            };

            canvas.addEventListener('pointerdown', event => {
                const box = canvas.getBoundingClientRect();
                lastPoint = {x: event.clientX - box.left, y: event.clientY - box.top};
                element.dataset.lastX = lastPoint.x.toFixed(2);
                element.dataset.lastY = lastPoint.y.toFixed(2);
                status.textContent = `x = ${lastPoint.x.toFixed(1)} px · y = ${lastPoint.y.toFixed(1)} px`;
                draw();
            });
            const observer = new ResizeObserver(draw);
            observer.observe(canvas);
            window.addEventListener('lcm:viewport-settled', draw);
            window.addEventListener('lcm:theme-changed', draw);
            draw();
        }
        """)
        return Bonito.jsrender(session, root)
    end
end

function register_presentation_routes!(server)
    Bonito.route!(server, "/presentations/probe" => presentation_probe_widget())
    return server
end

function presentation_source(argument::Union{Nothing,AbstractString}=nothing)
    if isnothing(argument)
        return DEFAULT_PRESENTATION
    end
    candidate = isabspath(argument) ? normpath(argument) : abspath(argument)
    if !isfile(candidate)
        candidate = normpath(joinpath(PLAYGROUND_ROOT, argument))
    end
    isfile(candidate) || throw(ArgumentError("presentation source not found: $argument"))
    endswith(lowercase(candidate), ".qmd") || throw(ArgumentError(
        "presentation source must be a .qmd file: $argument"
    ))
    relative = relpath(realpath(candidate), realpath(PLAYGROUND_ROOT))
    (relative == ".." || startswith(relative, "../") || startswith(relative, "..\\")) &&
        throw(ArgumentError("presentation source must live below $PLAYGROUND_ROOT"))
    return candidate
end

function presentation_output(source)
    relative = relpath(source, PLAYGROUND_ROOT)
    stem, _ = splitext(relative)
    return joinpath(SITE_DIR, "$stem.html")
end

function render_presentation(source=DEFAULT_PRESENTATION; quiet=false)
    source = presentation_source(source)
    executable = require_quarto()
    println("Rendering presentation with Quarto ($executable)")
    cd(PLAYGROUND_ROOT) do
        Quarto.render(source; execute=false, quiet)
    end
    output = presentation_output(source)
    isfile(output) || error("Quarto did not produce $output")
    println("Rendered $output")
    return output
end

function presentation_contract_errors(source, output=presentation_output(source))
    errors = String[]
    document = read(source, String)
    lowercase_document = lowercase(document)

    occursin("lcm-deck-revealjs", document) || push!(errors,
        "format must be lcm-deck-revealjs")
    occursin(r"(?im)^\s*#\s+", document) && push!(errors,
        "level-one Markdown headings are outside the flat-slide contract")
    occursin(r"(?is)<\s*(script|style|section)(?:\s|>)", document) && push!(errors,
        "raw script, style, and section tags are outside the authoring contract")
    occursin(r"(?is)<\s*iframe(?:\s|>)", document) && push!(errors,
        "use the bonito shortcode instead of authoring iframe markup")

    for match in eachmatch(r"\{\{<\s*bonito\b(.*?)>\}\}"s, document)
        body = match.captures[1]
        occursin(r"\broute\s*=\s*\"/", body) || push!(errors,
            "every bonito shortcode requires an absolute same-origin route")
        occursin(r"\bpublic-url\s*=\s*\"(?:/|https?://)", body) || push!(errors,
            "every presentation bonito shortcode requires public-url")
    end

    for match in eachmatch(r"\blcm-layout-([a-z0-9-]+)\b", lowercase_document)
        haskey(PRESENTATION_LAYOUTS, match.captures[1]) || push!(errors,
            "unknown presentation layout: $(match.captures[1])")
    end

    if !isfile(output)
        push!(errors, "rendered presentation is missing: $output")
        return errors
    end

    html = read(output, String)
    occursin("disableLayout: true", html) || push!(errors,
        "generated deck did not disable Reveal layout")
    occursin("lcm-deck-controller", html) || push!(errors,
        "generated deck is missing the LCM controller")
    live_count = count("class=\"lcm-live-viewport\"", html)
    frame_count = count("data-lcm-src=", html)
    placeholder_count = count("class=\"lcm-live-placeholder\"", html)
    live_count == frame_count == placeholder_count || push!(errors,
        "live frames and static placeholders are not one-to-one")
    occursin(r"<iframe[^>]*\ssrc\s*="i, html) && push!(errors,
        "presentation live frames must not load before the controller selects audience mode")
    occursin(r"(?is)<div[^>]*class=\"[^\"]*lcm-layout[^\"]*\"[^>]*>\s*<section", html) &&
        push!(errors, "a layout slot compiled into a nested section")
    parsed = EzXML.parsehtml(html)
    for layout in findall("//*[@data-lcm-layout]", root(parsed))
        name = layout["data-lcm-layout"]
        expected = get(PRESENTATION_LAYOUTS, name, 0)
        slots = findall(
            "./*[contains(concat(' ', normalize-space(@class), ' '), ' lcm-slot ')]",
            layout
        )
        length(slots) == expected || push!(errors,
            "rendered layout '$name' requires $expected direct slots; found $(length(slots))")
        isempty(findall(".//section", layout)) || push!(errors,
            "rendered layout '$name' contains a nested section")
    end
    occursin("scale(", read(joinpath(
        PLAYGROUND_ROOT, "_extensions", "lcm-deck", "deck.scss"
    ), String)) && push!(errors, "deck-owned CSS must not scale presentation geometry")
    return unique(errors)
end

function check_presentation(source=DEFAULT_PRESENTATION; quiet=false, render=true)
    source = presentation_source(source)
    output = render ? render_presentation(source; quiet) : presentation_output(source)
    errors = presentation_contract_errors(source, output)
    isempty(errors) || throw(ArgumentError(
        "presentation contract failed:\n  - " * join(errors, "\n  - ")
    ))
    println("Presentation contract passed: $source")
    return output
end

function presentation_browser()
    configured = get(ENV, "LCM_BROWSER", "")
    if !isempty(configured)
        executable = Sys.which(configured)
        isnothing(executable) && isfile(configured) && return configured
        isnothing(executable) || return executable
        throw(ArgumentError("LCM_BROWSER is not executable: $configured"))
    end
    for name in ("google-chrome", "chromium", "chromium-browser", "google-chrome-stable")
        executable = Sys.which(name)
        isnothing(executable) || return executable
    end
    throw(ArgumentError(
        "Chromium browser not found. Set LCM_BROWSER to a Chrome/Chromium executable."
    ))
end

function export_presentation(source=DEFAULT_PRESENTATION; output=nothing, quiet=false)
    source = presentation_source(source)
    html = check_presentation(source; quiet)
    pdf = isnothing(output) ? first(splitext(html)) * ".pdf" : abspath(output)
    mkpath(dirname(pdf))

    server = Bonito.Server(DEFAULT_HOST, 0)
    try
        register_static_site_routes!(server)
        route = "/" * replace(relpath(html, SITE_DIR), '\\' => '/')
        url = Bonito.online_url(server, route) * "?lcm-print"
        browser = presentation_browser()
        mktempdir() do profile
            command = `$browser --headless=new --disable-gpu --no-first-run --no-default-browser-check --run-all-compositor-stages-before-draw --virtual-time-budget=3500 --no-pdf-header-footer --user-data-dir=$profile --print-to-pdf=$pdf $url`
            run(command)
        end
    finally
        close(server)
    end
    isfile(pdf) || error("Chromium did not produce $pdf")
    println("Exported $pdf")
    return pdf
end

function split_presentation_arguments(arguments)
    if !isempty(arguments) && !startswith(arguments[1], "-")
        return presentation_source(arguments[1]), arguments[2:end]
    end
    return DEFAULT_PRESENTATION, arguments
end
