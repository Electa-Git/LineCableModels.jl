"""Session-local scientific views shared unchanged by the showcase and CableStudy."""
module ScientificViews

using Bonito, Observables, UUIDs
using ..Toolkit
using ..ComponentXRay
using ..LineCableModelsPlayground: RuntimeClient, WorkerSelector, PreparationStatus,
    ScientificJob, ScientificResult, WorkerDiagnostics, JuliaTerminal

include("StudyCases.jl")
using .StudyCases
export StudyCases, ScientificView, CableGeometry, StudyRuntime

const STYLE_PATH = normpath(joinpath(@__DIR__, "..", "..", "assets", "scientific-views.css"))
include_dependency(STYLE_PATH)
const STYLES = read(STYLE_PATH, String)
styles() = DOM.style(STYLES; var"data-lcm-css-source"="assets/scientific-views.css")

number(id, name, label, value, low, high, step, unit) = Field(label,
    UnitNumberInput(name; id="$id-$name", value, minimum=low, maximum=high, step, unit, required=true))

function frequency_fields(id, low, high, points)
    return (minimum_frequency_hz=number(id, :minimum_frequency_hz, "From", low, 1, 1e6, 1, "Hz"),
        maximum_frequency_hz=number(id, :maximum_frequency_hz, "To", high, 1, 1e6, 1, "Hz"),
        frequency_points=number(id, :frequency_points, "Samples", points, 2, 200, 1, ""))
end

function case_fields(case::LineParameters, id)
    spec = inputs(case)
    frequencies = spec["frequencies_hz"]
    return (; separation_m=number(id, :separation_m, "Cable separation", spec["separation_m"], 0.05, 100, 0.05, "m"),
        depth_m=number(id, :depth_m, "Burial depth", spec["depth_m"], 0.05, 100, 0.05, "m"),
        earth_resistivity_ohm_m=number(id, :earth_resistivity_ohm_m, "Earth resistivity", spec["earth_resistivity_ohm_m"], 0.01, 1e6, 0.01, "Ω m"),
        frequency_fields(id, first(frequencies), last(frequencies), length(frequencies))...)
end

function case_fields(case::CorridorImpedance, id)
    spec = inputs(case)
    return (; ugc_share=number(id, :ugc_share, "Underground share", spec["ugc_share"], 0.001, 0.999, 0.001, "p.u."),
        corridor_length_m=number(id, :corridor_length_m, "Corridor length", spec["corridor_length_m"], 1000, 1e6, 1000, "m"),
        length_error_percent=number(id, :length_error_percent, "Length variation ±", spec["length_error_percent"], 0, 50, 1, "%"),
        frequency_fields(id, spec["minimum_frequency_hz"], spec["maximum_frequency_hz"], spec["frequency_points"])...)
end

quantities(::LineParameters) = ["resistance"=>"Re Z · resistance", "reactance"=>"Im Z · reactance",
    "conductance"=>"Re Y · conductance", "susceptance"=>"Im Y · susceptance"]
quantities(::CorridorImpedance) = ["impedance"=>"Driving-point impedance at B5"]
curve_labels(::LineParameters) = ("Self (1,1)", "Mutual (1,2)")
curve_labels(::CorridorImpedance) = ("Shorter corridor", "Nominal corridor", "Longer corridor")
profile(::LineParameters) = "line-parameters"
profile(::CorridorImpedance) = "power-flow"
assumptions(::LineParameters) = "Identical coaxial cables: 10 mm core, 5 mm insulation, 1 mm sheath; 20 °C. Grounded sheaths, homogeneous earth. Values are per metre, not statistical envelopes."
assumptions(::CorridorImpedance) = "Reference case ohl_ugc_transition_v1; earth 100 Ω m. Reuses the prepared active-device linearization. The three curves vary passive corridor length only; they are not confidence intervals."
preparation_note(::LineParameters) = "Prepare runs representative endpoint frequencies. It warms the executor; Run still evaluates your requested inputs."
preparation_note(::CorridorImpedance) = "Prepare solves and linearizes the reference network. Changing corridor inputs reuses that model; it does not silently change active-device setpoints."

"""
    ScientificView(session, case, client)

Compose shared typed fields, explicit scientific-job controls and a persistent
SVG/table result view for `LineParameters()` or `CorridorImpedance()`. Construction
and field changes never prepare or submit work. Bindings belong to `session`;
the shared `ScientificJob` retains provenance and fences outdated results.
"""
struct ScientificView{C<:AbstractStudyCase,F,J}
    "Approved scientific case; dispatch selects the operation and projection."
    case::C
    "Owned run context, not an execution authority."
    client::RuntimeClient
    "Named collection of the existing typed fields."
    fields::F
    "Browser validity state; incomplete edits disable submission."
    valid::Observable{Bool}
    "Display-only selection; changing it does not submit work."
    quantity::ComboBox
    "Shared fenced scientific job and last-good result binding."
    job::J
end

function ScientificView(session::Session, case::AbstractStudyCase, client::RuntimeClient)
    id = "study-$(uuid4())"
    fields = case_fields(case, id)
    valid = Observable(true)
    parameters = Observable{Any}(inputs(case))
    map!(session, parameters, valid, (f.control.value for f in values(fields))...) do accepted, vals...
        accepted || return nothing
        try
            return inputs(case; NamedTuple{keys(fields)}(vals)...)
        catch exception
            exception isa ArgumentError || rethrow()
            return nothing
        end
    end
    quantity = ComboBox(:quantity, quantities(case); id="$id-quantity")
    job = ScientificJob(client, role(case), operation(case); parameters)
    return ScientificView(case, client, fields, valid, quantity, job)
end

"""Keep the empty and populated scientific display states in one Observable type."""
struct ScientificDisplay
    "Optional bounded series with its case-declared frequency and ordinate units."
    series::Union{Nothing,NamedTuple}
    "Plain-language current, outdated or unavailable display status."
    status::String
end

function display_state(view, result, quantity)
    result === nothing && return ScientificDisplay(nothing, "No calculation yet. Select a worker, prepare, then Run.")
    series = try
        result_series(view.case, result.value, quantity)
    catch
        return ScientificDisplay(nothing, "Result cannot be displayed: unexpected shape or nonfinite values. The job receipt remains available below.")
    end
    status = result.current ? "Current result · worker provenance below" : "Outdated result · retained for comparison; Run again for current inputs and assignment"
    return ScientificDisplay(series, status)
end

tick(value) = string(round(value; sigdigits=4))

# A bounded display projection, not a chart framework: fixed axes and persistent
# SVG nodes. Normalize before subtraction to avoid overflow on finite wire data.
function plot_coordinates(series)
    series === nothing && return (paths=String[], xticks=String[], yticks=String[], unit="")
    x = log10.(series.frequency)
    all_values = reduce(vcat, (curve.values for curve in series.curves))
    scale = max(maximum(abs, all_values), eps(Float64))
    low, high = extrema(all_values ./ scale)
    padding = max((high-low)*0.06, 0.02)
    low = max(-1.0, low-padding); high = min(1.0, high+padding)
    xp(v) = 78 + 600*(v-first(x))/(last(x)-first(x))
    yp(v) = 282 - 236*(v/scale-low)/(high-low)
    paths = [join(("$(xp(x[i])),$(yp(curve.values[i]))" for i in eachindex(x)), ' ') for curve in series.curves]
    return (paths=paths, xticks=tick.(exp10.(range(first(x), last(x); length=5))),
        yticks=tick.(range(low, high; length=5) .* scale), unit=series.unit)
end

function result_plot(session, view, state)
    coordinates = map(s -> plot_coordinates(s.series), session, state)
    lines = [SVG.polyline(; points="", var"data-study-curve"=i,
        class="lc-study-curve lc-study-series-$i") for i in eachindex(curve_labels(view.case))]
    # Observable children render HTML wrappers in Bonito. Inside SVG those
    # wrappers close the foreign-content context, displacing following labels.
    # Keep text nodes native SVG and update their textContent in place instead.
    xticks = [SVG.text(""; var"data-study-xtick"=i,
        x=78+150*(i-1), y=305, var"text-anchor"="middle") for i in 1:5]
    yticks = [SVG.text(""; var"data-study-ytick"=i,
        x=69, y=286-59*(i-1), var"text-anchor"="end") for i in 1:5]
    svg = SVG.svg(SVG.title("Frequency response; logarithmic frequency axis"),
        SVG.path(; d="M78 46 V282 H678", class="lc-study-axis"), lines..., xticks..., yticks...,
        SVG.text("Frequency (Hz) · logarithmic"; x=378, y=336, var"text-anchor"="middle"),
        SVG.text(""; var"data-study-unit"="", x=78, y=24);
        viewBox="0 0 720 350", role="img", class="lc-study-plot")
    Bonito.onload(session, svg, js"""
        element => {
            const coordinates = $(coordinates);
            const update = value => {
                // SVG geometry uses attributes, not the read-only animated DOM
                // properties used by the generic HTML Observable binding.
                element.querySelectorAll('[data-study-curve]').forEach((node, i) => {
                    node.setAttribute('points', value.paths[i] || '');
                });
                for (const axis of ['x', 'y']) element.querySelectorAll(`[data-study-${axis}tick]`).forEach((node, i) => {
                    node.textContent = value[axis + 'ticks'][i] || '';
                });
                element.querySelector('[data-study-unit]').textContent = value.unit;
            };
            coordinates.on(update);
            update(coordinates.value);
        }
    """)
    legend = DOM.div((DOM.span(label; class="lc-study-series-$i")
        for (i, label) in enumerate(curve_labels(view.case)))...; class="lc-study-legend")
    table = DataTable(TableColumn(:frequency, "Frequency (Hz)"; format=tick, align=:right),
        (TableColumn(Symbol("series$i"), label; format=tick, align=:right)
            for (i,label) in enumerate(curve_labels(view.case)))...; filterable=false, label="Scientific samples")
    map!(session, table.rows, state) do snapshot
        s = snapshot.series
        s === nothing && return Any[]
        return Any[Dict{Symbol,Any}(:frequency=>s.frequency[i],
            (Symbol("series$j")=>curve.values[i] for (j,curve) in enumerate(s.curves))...)
            for i in eachindex(s.frequency)]
    end
    return DOM.div(ViewportFrame("Scientific response", svg; footer=legend),
        DOM.p(map(s -> s.status, session, state); class="lc-study-result-status", role="status"),
        Disclosure("Samples · values in the plotted unit", table); class="lc-study-result")
end

function Bonito.jsrender(session::Session, view::ScientificView)
    state = map((result, quantity) -> display_state(view, result, quantity), session, view.job.result, view.quantity.value)
    validation = map(session, view.job.parameters) do parameters
        parameters === nothing ? "Complete all fields within their limits. From must be less than To; Samples must be an integer." : "Input edits do not run a calculation."
    end
    fields = DOM.div(values(view.fields)...; class="lc-study-fields", var"data-runtime-input-fields"="",
        oninput=js"event => $(view.valid).notify([...event.currentTarget.querySelectorAll('input')].every(input => input.checkValidity()))")
    setup = Disclosure("Worker selection and preparation", DOM.div(
        DOM.p(preparation_note(view.case); class="lc-study-note"),
        WorkerSelector(view.client, role(view.case); profiles=(profile(view.case),)),
        PreparationStatus(view.client, role(view.case); parameters=preparation_inputs(view.case))))
    node = DOM.section(DOM.header(DOM.h2(case_title(view.case)), DOM.p(assumptions(view.case); class="lc-study-note")),
        setup, DOM.div(DOM.div(view.job, fields, DOM.p(validation; class="lc-study-note", role="status"),
            Field("Display quantity", view.quantity); class="lc-study-inputs"), result_plot(session, view, state);
            class="lc-study-layout"); class="lc-study-view", var"data-runtime-input-scope"="")
    return Bonito.jsrender(session, DOM.div(styles(), ComponentXRay.instrument(session, node, view); style="display: contents;"))
end

function ComponentXRay.inspection(view::ScientificView)
    return ComponentXRay.ComponentInspection(view; name="ScientificView",
        source=ComponentXRay.source_reference(@__MODULE__, @__FILE__, @__LINE__),
        parameters=[ComponentXRay.PropertyInspection(:case, nameof(typeof(view.case))),
            ComponentXRay.PropertyInspection(:role, role(view.case)),
            ComponentXRay.PropertyInspection(:operation, operation(view.case))],
        css_scopes=[".lc-study-view", ".lc-study-layout", ".lc-study-fields", ".lc-study-inputs",
            ".lc-study-note", ".lc-study-result", ".lc-study-result-status", ".lc-study-plot", ".lc-study-axis",
            ".lc-study-curve", ".lc-study-legend", ".lc-study-series-1", ".lc-study-series-2", ".lc-study-series-3"],
        notes=["Parameters, result payloads and run identities are not diagnostic metadata.",
            "Typed fields and ScientificJob own their actual callbacks and bindings."])
end

include("CableGeometry.jl")
include("StudyRuntime.jl")

end
