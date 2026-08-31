module ICHQP2026LineParametersDeck

using Bonito
import LineCableModels as LCM
using Makie
using Markdown
using Measurements: Measurement, measurement, uncertainty, value
using Observables: off
using Printf
using WGLMakie
using ..PageAuthoring

const FREQUENCIES = collect(
    10.0 .^ range(log10(0.5), log10(50.0 * 50); length = 72)
)
const CONDUCTOR_DIAMETER_RANGE_MM = 10.0:0.5:50.0
const INSULATION_THICKNESS_RANGE_MM = 2.0:0.25:20.0
const LATERAL_DISTANCE_RANGE_M = 1.0:0.25:10.0
const EARTH_RESISTIVITY_LOG_RANGE = 1.0:0.025:3.0
const LINE_LENGTH_RANGE_M = 100.0:100.0:5000.0
const UNCERTAINTY_RANGE_PCT = 0.0:0.25:10.0
const BURIAL_DEPTH_M = 1.0
const DEFAULT_INPUT = (
    conductor_diameter = 25.0e-3,
    conductor_uncertainty = 1.0,
    insulation_thickness = 8.0e-3,
    insulation_uncertainty = 1.0,
    lateral_distance = 4.0,
    distance_uncertainty = 1.0,
    earth_resistivity = 100.0,
    earth_uncertainty = 1.0,
    line_length = 1000.0,
    length_uncertainty = 1.0
)
const COPPER = LCM.Material(
    :conductor,
    1.7241e-8,
    1.0,
    0.999994,
    20.0,
    0.00393
)
const XLPE = LCM.Material(
    :insulator,
    1.0e14,
    2.3,
    1.0,
    20.0,
    0.0;
    tan_delta = 2.0e-4
)
const FIXED_CABLE = (
    conductor_diameter = 25.0e-3,
    insulation_thickness = 8.0e-3
)
const LINE_FORMULATION = LCM.Formulation(
    options = (ideal_transposition = false,)
)
const PLOT_BACKGROUND = RGBf(0.9, 0.9, 0.9)
const PLOT_TEXT = RGBf(0.16, 0.18, 0.18)
const CURVE_COLOR = RGBf(0.13, 0.38, 0.88)
const BAND_COLOR = RGBAf(0.13, 0.38, 0.88, 0.22)
const QUANTITIES = (:R, :L, :C, :G)
const QUANTITY_FUNCTIONS = (
    R = LCM.R,
    L = LCM.L,
    C = LCM.C,
    G = LCM.G
)
const SELF_UNITS = (
    R = LCM.display_unit(LCM.R, :total),
    L = LCM.display_unit(LCM.L, :total),
    C = LCM.display_unit(LCM.C, :total),
    G = LCM.display_unit(LCM.G, :total; prefix = :micro)
)
const MUTUAL_UNITS = (
    R = LCM.display_unit(LCM.R, :total),
    L = LCM.display_unit(LCM.L, :total),
    C = LCM.display_unit(LCM.C, :total; prefix = :nano),
    G = LCM.display_unit(LCM.G, :total; prefix = :micro)
)
const PREFLIGHT_STATE = preflight_state(
    label = "Line-parameter dashboards not activated"
)
const WARMUP_TASK = Ref{Union{Nothing, Task}}(nothing)
const WARMUP_RESULT = Ref{Any}(nothing)
const WARMUP_LOCK = ReentrantLock()

function positive_finite(name::AbstractString, number::Real)
    isfinite(number) && number > zero(number) || throw(DomainError(
        number,
        "$name must be positive and finite"
    ))
    return number
end

function uncertainty_percent(name::AbstractString, number::Real)
    isfinite(number) && zero(number) <= number <= oftype(number, 100) ||
        throw(DomainError(number, "$name must lie between 0 and 100 percent"))
    return number
end

function measured(nominal::Real, uncertainty_pct::Real)
    positive_finite("nominal value", nominal)
    uncertainty_percent("standard uncertainty", uncertainty_pct)
    return measurement(nominal, nominal * uncertainty_pct / 100)
end

"""Build one package-owned solid-core cable with an XLPE outer shell."""
function build_cable(
        cable_id::AbstractString,
        conductor_diameter_m::Real,
        insulation_thickness_m::Real
)
    positive_finite("conductor diameter", conductor_diameter_m)
    positive_finite("insulation thickness", insulation_thickness_m)
    core = LCM.Group(
        :core,
        LCM.Region(
            :core_metal,
            LCM.Disk(conductor_diameter_m / 2),
            COPPER
        )
    )
    insulation = LCM.Region(
        :core_insulation,
        LCM.Shell(insulation_thickness_m),
        XLPE
    )
    return LCM.build(
        LCM.CableDesign,
        cable_id,
        LCM.Stack(core, insulation)
    )
end

"""
Build the uncertain two-cable line problem entirely through LineCableModels.

Only cable 1 is geometrically uncertain. Both cable centres remain at 1 m
depth and move symmetrically as their uncertain lateral separation changes.
The second cable is the fixed reference design used by both matrix views.
"""
function build_problem(input; frequencies = FREQUENCIES)
    cable_1 = build_cable(
        "ichqp2026-line-cable-1",
        measured(input.conductor_diameter, input.conductor_uncertainty),
        measured(input.insulation_thickness, input.insulation_uncertainty)
    )
    cable_2 = build_cable(
        "ichqp2026-line-cable-2",
        FIXED_CABLE.conductor_diameter,
        FIXED_CABLE.insulation_thickness
    )
    distance = measured(
        input.lateral_distance,
        input.distance_uncertainty
    )
    line_length = measured(input.line_length, input.length_uncertainty)
    earth = LCM.Earth(
        rho = measured(input.earth_resistivity, input.earth_uncertainty),
        eps_r = 10.0,
        mu_r = 1.0
    )
    system = LCM.build(
        LCM.LineCableSystem,
        [cable_1, cable_2],
        [
            LCM.Pose2(-distance / 2, -BURIAL_DEPTH_M),
            LCM.Pose2(distance / 2, -BURIAL_DEPTH_M)
        ];
        connections = ((core = 1,), (core = 2,)),
        environment = earth,
        system_id = "ichqp2026-two-cable-line",
        line_length
    )
    problem = LCM.LineParametersProblem(
        system;
        temperature = 20.0,
        earth_props = earth,
        frequencies
    )
    return (; cable_1, cable_2, earth, system, problem)
end

function curve_units(element::Symbol)
    element === :self && return SELF_UNITS
    element === :mutual && return MUTUAL_UNITS
    throw(ArgumentError("unknown line-parameter matrix element: $element"))
end

function curve_payload(parameters::LCM.LineParameters, element::Symbol)
    indices = element === :self ? (1, 1) :
              element === :mutual ? (1, 2) :
              throw(ArgumentError("unknown matrix element: $element"))
    units = curve_units(element)
    curves = map(QUANTITIES) do quantity
        observable = getproperty(QUANTITY_FUNCTIONS, quantity)
        unit = getproperty(units, quantity)
        factor = LCM.scale_factor(observable, :total, unit)
        measured_values = LCM.observe(parameters, observable, indices...)
        nominal_values = Float64.(value.(measured_values)) .* factor
        sigma_values = Float64.(uncertainty.(measured_values)) .* abs(factor)
        return (
            nominal = nominal_values,
            lower = nominal_values .- sigma_values,
            upper = nominal_values .+ sigma_values,
            unit = LCM.label(unit)
        )
    end
    return NamedTuple{QUANTITIES}(curves)
end

function compute_line_parameters(input; frequencies = FREQUENCIES)
    started = time_ns()
    model = build_problem(input; frequencies)
    parameters = LCM.compute(
        model.problem,
        LINE_FORMULATION;
        options = (output_basis = :total,)
    )
    elapsed_ms = (time_ns() - started) / 1.0e6
    frequency_values = Float64.(value.(LCM.frequencies(parameters)))
    return (;
        input,
        model,
        parameters,
        frequencies = frequency_values,
        curves = (
            self = curve_payload(parameters, :self),
            mutual = curve_payload(parameters, :mutual)
        ),
        elapsed_ms,
        finished_at = time()
    )
end

function quantity_title(quantity::Symbol, element::Symbol)
    suffix = element === :self ? "₁₁" : "₁₂"
    return "$(quantity)$suffix(f)"
end

function line_parameter_figure(result, element::Symbol)
    payload = getproperty(result.curves, element)
    frequencies = Observable(copy(result.frequencies))
    figure = Figure(
        size = (1100, 720),
        backgroundcolor = PLOT_BACKGROUND,
        figure_padding = 12
    )
    axes = Axis[]
    curves = map(enumerate(QUANTITIES)) do (index, quantity)
        row = index <= 2 ? 1 : 2
        column = isodd(index) ? 1 : 2
        initial = getproperty(payload, quantity)
        nominal = Observable(copy(initial.nominal))
        lower = Observable(copy(initial.lower))
        upper = Observable(copy(initial.upper))
        axis = Axis(
            figure[row, column];
            title = quantity_title(quantity, element),
            xlabel = row == 2 ? "Frequency (Hz)" : "",
            ylabel = "$(quantity) ($(initial.unit))",
            xscale = log10,
            backgroundcolor = PLOT_BACKGROUND,
            titlecolor = PLOT_TEXT,
            xlabelcolor = PLOT_TEXT,
            ylabelcolor = PLOT_TEXT,
            xticklabelcolor = PLOT_TEXT,
            yticklabelcolor = PLOT_TEXT,
            leftspinecolor = RGBf(0.35, 0.38, 0.38),
            rightspinecolor = RGBf(0.35, 0.38, 0.38),
            topspinecolor = RGBf(0.35, 0.38, 0.38),
            bottomspinecolor = RGBf(0.35, 0.38, 0.38),
            xgridcolor = RGBAf(0.35, 0.38, 0.38, 0.28),
            ygridcolor = RGBAf(0.35, 0.38, 0.38, 0.28),
            titlesize = 18,
            xlabelsize = 16,
            ylabelsize = 16,
            xticklabelsize = 13,
            yticklabelsize = 13
        )
        band_plot = band!(
            axis,
            frequencies,
            lower,
            upper;
            color = BAND_COLOR
        )
        line_plot = lines!(
            axis,
            frequencies,
            nominal;
            color = CURVE_COLOR,
            linewidth = 2.4
        )
        push!(axes, axis)
        return (; axis, nominal, lower, upper, band_plot, line_plot)
    end
    named_curves = NamedTuple{QUANTITIES}(Tuple(curves))
    linkxaxes!(axes...)
    colgap!(figure.layout, 14)
    rowgap!(figure.layout, 12)
    inspector = DataInspector(figure)
    return (; figure, axes, frequencies, curves = named_curves, inspector, element)
end

function update_line_parameter_figure!(handle, result)
    payload = getproperty(result.curves, handle.element)
    handle.frequencies[] = copy(result.frequencies)
    for quantity in QUANTITIES
        target = getproperty(handle.curves, quantity)
        source = getproperty(payload, quantity)
        target.lower[] = copy(source.lower)
        target.upper[] = copy(source.upper)
        target.nominal[] = copy(source.nominal)
    end
    foreach(autolimits!, handle.axes)
    return handle
end

function reset_line_parameter_views!(handles...)
    for handle in handles
        foreach(autolimits!, handle.axes)
    end
    return nothing
end

mutable struct SweepWorker
    queue::Channel{Any}
    generation::Int
    closed::Bool
    lock::ReentrantLock
    task::Union{Nothing, Task}
end

SweepWorker() = SweepWorker(Channel{Any}(1), 0, false, ReentrantLock(), nothing)

function enqueue!(worker::SweepWorker, input)
    return lock(worker.lock) do
        worker.closed && return worker.generation
        worker.generation += 1
        while isready(worker.queue)
            take!(worker.queue)
        end
        put!(worker.queue, (generation = worker.generation, input))
        return worker.generation
    end
end

function worker_is_current(worker::SweepWorker, generation::Integer)
    return lock(worker.lock) do
        !worker.closed && worker.generation == generation
    end
end

function close!(worker::SweepWorker)
    lock(worker.lock) do
        worker.closed && return nothing
        worker.closed = true
        close(worker.queue)
    end
    return nothing
end

function run_worker!(worker::SweepWorker, result_state, status_state)
    while !worker.closed
        item = try
            take!(worker.queue)
        catch error
            error isa InvalidStateException && break
            rethrow()
        end

        sleep(0.22)
        while isready(worker.queue)
            item = take!(worker.queue)
        end
        worker_is_current(worker, item.generation) || continue
        status_state[] = (
            phase = :computing,
            label = "Running one shared uncertain frequency sweep…"
        )
        try
            calculation = Threads.@spawn compute_line_parameters(item.input)
            result = fetch(calculation)
            worker_is_current(worker, item.generation) || continue
            result_state[] = result
            status_state[] = (
                phase = :ready,
                label = @sprintf(
                    "Package sweep ready · %d points · %.0f ms",
                    length(result.frequencies),
                    result.elapsed_ms
                )
            )
        catch error
            worker_is_current(worker, item.generation) || continue
            status_state[] = (
                phase = :error,
                label = "Line-parameter update failed: $(sprint(showerror, error))"
            )
        end
    end
    return nothing
end

function make_slider(values, default, id, label)
    return Bonito.Slider(
        collect(values);
        value = default,
        id,
        ariaLabel = label
    )
end

const CONTROL_NAMES = (
    :conductor_diameter,
    :conductor_uncertainty,
    :insulation_thickness,
    :insulation_uncertainty,
    :lateral_distance,
    :distance_uncertainty,
    :earth_resistivity,
    :earth_uncertainty,
    :line_length,
    :length_uncertainty
)
const CONTROL_DEFAULTS = (
    conductor_diameter = 1.0e3 * DEFAULT_INPUT.conductor_diameter,
    conductor_uncertainty = DEFAULT_INPUT.conductor_uncertainty,
    insulation_thickness = 1.0e3 * DEFAULT_INPUT.insulation_thickness,
    insulation_uncertainty = DEFAULT_INPUT.insulation_uncertainty,
    lateral_distance = DEFAULT_INPUT.lateral_distance,
    distance_uncertainty = DEFAULT_INPUT.distance_uncertainty,
    earth_resistivity = log10(DEFAULT_INPUT.earth_resistivity),
    earth_uncertainty = DEFAULT_INPUT.earth_uncertainty,
    line_length = DEFAULT_INPUT.line_length,
    length_uncertainty = DEFAULT_INPUT.length_uncertainty
)

function reset_control_bundle!(controls)
    reset_selectors!((
        getproperty(controls, name) => getproperty(CONTROL_DEFAULTS, name)
    for name in CONTROL_NAMES)...)
    return nothing
end

function make_control_bundle(prefix::AbstractString)
    return (
        conductor_diameter = make_slider(
            CONDUCTOR_DIAMETER_RANGE_MM,
            1.0e3 * DEFAULT_INPUT.conductor_diameter,
            "$prefix-conductor-diameter-input",
            "Cable 1 conductor diameter in millimetres"
        ),
        conductor_uncertainty = make_slider(
            UNCERTAINTY_RANGE_PCT,
            DEFAULT_INPUT.conductor_uncertainty,
            "$prefix-conductor-uncertainty-input",
            "Cable 1 conductor-diameter standard uncertainty in percent"
        ),
        insulation_thickness = make_slider(
            INSULATION_THICKNESS_RANGE_MM,
            1.0e3 * DEFAULT_INPUT.insulation_thickness,
            "$prefix-insulation-thickness-input",
            "Cable 1 insulation thickness in millimetres"
        ),
        insulation_uncertainty = make_slider(
            UNCERTAINTY_RANGE_PCT,
            DEFAULT_INPUT.insulation_uncertainty,
            "$prefix-insulation-uncertainty-input",
            "Cable 1 insulation-thickness standard uncertainty in percent"
        ),
        lateral_distance = make_slider(
            LATERAL_DISTANCE_RANGE_M,
            DEFAULT_INPUT.lateral_distance,
            "$prefix-lateral-distance-input",
            "Cable-centre lateral distance in metres"
        ),
        distance_uncertainty = make_slider(
            UNCERTAINTY_RANGE_PCT,
            DEFAULT_INPUT.distance_uncertainty,
            "$prefix-distance-uncertainty-input",
            "Lateral-distance standard uncertainty in percent"
        ),
        earth_resistivity = make_slider(
            EARTH_RESISTIVITY_LOG_RANGE,
            log10(DEFAULT_INPUT.earth_resistivity),
            "$prefix-earth-resistivity-input",
            "Earth resistivity on a logarithmic scale"
        ),
        earth_uncertainty = make_slider(
            UNCERTAINTY_RANGE_PCT,
            DEFAULT_INPUT.earth_uncertainty,
            "$prefix-earth-uncertainty-input",
            "Earth-resistivity standard uncertainty in percent"
        ),
        line_length = make_slider(
            LINE_LENGTH_RANGE_M,
            DEFAULT_INPUT.line_length,
            "$prefix-line-length-input",
            "Line length in metres"
        ),
        length_uncertainty = make_slider(
            UNCERTAINTY_RANGE_PCT,
            DEFAULT_INPUT.length_uncertainty,
            "$prefix-length-uncertainty-input",
            "Line-length standard uncertainty in percent"
        )
    )
end

function link_slider_values!(session::Session, left, right)
    left_to_right = on(session, left.value) do number
        isequal(right.value[], number) || (right.value[] = number)
        return nothing
    end
    right_to_left = on(session, right.value) do number
        isequal(left.value[], number) || (left.value[] = number)
        return nothing
    end
    return (; left_to_right, right_to_left)
end

function setup(session::Session)
    self_controls = make_control_bundle("line-self")
    mutual_controls = make_control_bundle("line-mutual")
    linked_controls = map(CONTROL_NAMES) do name
        link_slider_values!(
            session,
            getproperty(self_controls, name),
            getproperty(mutual_controls, name)
        )
    end

    input = map(
        (
            conductor_diameter,
            conductor_uncertainty,
            insulation_thickness,
            insulation_uncertainty,
            lateral_distance,
            distance_uncertainty,
            earth_resistivity,
            earth_uncertainty,
            line_length,
            length_uncertainty
        ) -> (;
            conductor_diameter = 1.0e-3 * Float64(conductor_diameter),
            conductor_uncertainty = Float64(conductor_uncertainty),
            insulation_thickness = 1.0e-3 * Float64(insulation_thickness),
            insulation_uncertainty = Float64(insulation_uncertainty),
            lateral_distance = Float64(lateral_distance),
            distance_uncertainty = Float64(distance_uncertainty),
            earth_resistivity = 10.0^Float64(earth_resistivity),
            earth_uncertainty = Float64(earth_uncertainty),
            line_length = Float64(line_length),
            length_uncertainty = Float64(length_uncertainty)
        ),
        session,
        (getproperty(self_controls, name).value for name in CONTROL_NAMES)...
    )

    cached = WARMUP_RESULT[]
    initial_result = isnothing(cached) ? compute_line_parameters(input[]) : cached
    result_state = Observable(initial_result)
    status_state = Observable((
        phase = :ready,
        label = @sprintf(
            "Package sweep ready · %d points · %.0f ms",
            length(initial_result.frequencies),
            initial_result.elapsed_ms
        )
    ))
    self_plot = line_parameter_figure(initial_result, :self)
    mutual_plot = line_parameter_figure(initial_result, :mutual)
    result_observer = on(result_state) do result
        update_line_parameter_figure!(self_plot, result)
        update_line_parameter_figure!(mutual_plot, result)
        return nothing
    end

    worker = SweepWorker()
    worker.task = @async run_worker!(worker, result_state, status_state)
    errormonitor(worker.task)
    input_observer = on(input) do latest
        status_state[] = (
            phase = :queued,
            label = "Waiting for controls to settle…"
        )
        enqueue!(worker, latest)
        return nothing
    end

    self_reset_button = Bonito.Button(
        "Reset all plot views";
        style = nothing,
        id = "line-self-reset-view-control",
        class = "lc-action-button",
        type = "button",
        ariaLabel = "Reset self and mutual line-parameter plot views"
    )
    mutual_reset_button = Bonito.Button(
        "Reset all plot views";
        style = nothing,
        id = "line-mutual-reset-view-control",
        class = "lc-action-button",
        type = "button",
        ariaLabel = "Reset self and mutual line-parameter plot views"
    )
    self_selector_reset_button = Bonito.Button(
        "Reset selectors";
        style = nothing,
        id = "line-self-reset-selectors-control",
        class = "lc-action-button",
        type = "button",
        ariaLabel = "Reset all self line-parameter selectors to their initial values"
    )
    mutual_selector_reset_button = Bonito.Button(
        "Reset selectors";
        style = nothing,
        id = "line-mutual-reset-selectors-control",
        class = "lc-action-button",
        type = "button",
        ariaLabel = "Reset all mutual line-parameter selectors to their initial values"
    )
    self_reset_observer = on(self_reset_button.value) do clicked
        clicked && reset_line_parameter_views!(self_plot, mutual_plot)
        return nothing
    end
    mutual_reset_observer = on(mutual_reset_button.value) do clicked
        clicked && reset_line_parameter_views!(self_plot, mutual_plot)
        return nothing
    end
    self_selector_reset_observer = on(self_selector_reset_button.value) do clicked
        clicked && reset_control_bundle!(self_controls)
        return nothing
    end
    mutual_selector_reset_observer = on(
        mutual_selector_reset_button.value
    ) do clicked
        clicked && reset_control_bundle!(mutual_controls)
        return nothing
    end

    on(session.on_close) do _
        off(result_observer)
        off(input_observer)
        off(self_reset_observer)
        off(mutual_reset_observer)
        off(self_selector_reset_observer)
        off(mutual_selector_reset_observer)
        close!(worker)
        return nothing
    end
    return (;
        input,
        result_state,
        status_state,
        self_controls,
        mutual_controls,
        linked_controls,
        self_plot,
        mutual_plot,
        worker,
        self_reset_button,
        mutual_reset_button,
        self_selector_reset_button,
        mutual_selector_reset_button
    )
end

function control_values(session::Session, input)
    return (
        conductor_diameter = map(
            state -> @sprintf("%.1f mm", 1.0e3 * state.conductor_diameter),
            session,
            input
        ),
        conductor_uncertainty = map(
            state -> @sprintf("±%.2g %%", state.conductor_uncertainty),
            session,
            input
        ),
        insulation_thickness = map(
            state -> @sprintf("%.2f mm", 1.0e3 * state.insulation_thickness),
            session,
            input
        ),
        insulation_uncertainty = map(
            state -> @sprintf("±%.2g %%", state.insulation_uncertainty),
            session,
            input
        ),
        lateral_distance = map(
            state -> @sprintf("%.2f m", state.lateral_distance),
            session,
            input
        ),
        distance_uncertainty = map(
            state -> @sprintf("±%.2g %%", state.distance_uncertainty),
            session,
            input
        ),
        earth_resistivity = map(
            state -> @sprintf("%.1f Ω·m", state.earth_resistivity),
            session,
            input
        ),
        earth_uncertainty = map(
            state -> @sprintf("±%.2g %%", state.earth_uncertainty),
            session,
            input
        ),
        line_length = map(
            state -> @sprintf("%.0f m", state.line_length),
            session,
            input
        ),
        length_uncertainty = map(
            state -> @sprintf("±%.2g %%", state.length_uncertainty),
            session,
            input
        )
    )
end

function parameter_pair(
        label::AbstractString,
        nominal_widget,
        uncertainty_widget,
        nominal_value,
        uncertainty_value;
        id::AbstractString
)
    return DOM.div(
        DOM.div(
            DOM.strong(label),
            DOM.span(
                DOM.output(nominal_value),
                DOM.output(uncertainty_value);
                class = "lc-line-parameter-pair-values"
            );
            class = "lc-control-heading"
        ),
        DOM.div(
            DOM.label(
                DOM.span("nominal"),
                nominal_widget;
                class = "lc-line-parameter-slider"
            ),
            DOM.label(
                DOM.span("standard uncertainty"),
                uncertainty_widget;
                class = "lc-line-parameter-slider"
            );
            class = "lc-line-parameter-slider-pair"
        );
        id,
        class = "lc-control lc-line-parameter-pair"
    )
end

function controls_panel(
        session::Session,
        state,
        widgets,
        reset_button,
        selector_reset_button;
        prefix::AbstractString
)
    values = control_values(session, state.input)
    status = map(snapshot -> snapshot.label, session, state.status_state)
    return webpart(
        parameter_pair(
            "Cable 1 conductor diameter",
            widgets.conductor_diameter,
            widgets.conductor_uncertainty,
            values.conductor_diameter,
            values.conductor_uncertainty;
            id = "$prefix-conductor-control"
        ),
        parameter_pair(
            "Cable 1 insulation thickness",
            widgets.insulation_thickness,
            widgets.insulation_uncertainty,
            values.insulation_thickness,
            values.insulation_uncertainty;
            id = "$prefix-insulation-control"
        ),
        parameter_pair(
            "Lateral centre distance",
            widgets.lateral_distance,
            widgets.distance_uncertainty,
            values.lateral_distance,
            values.distance_uncertainty;
            id = "$prefix-distance-control"
        ),
        parameter_pair(
            "Earth resistivity",
            widgets.earth_resistivity,
            widgets.earth_uncertainty,
            values.earth_resistivity,
            values.earth_uncertainty;
            id = "$prefix-earth-control"
        ),
        parameter_pair(
            "Line length",
            widgets.line_length,
            widgets.length_uncertainty,
            values.line_length,
            values.length_uncertainty;
            id = "$prefix-length-control"
        ),
        value_list(
            "Burial depth" => "$(BURIAL_DEPTH_M) m",
            "Frequency scan" => "0.5 Hz–2.5 kHz",
            "Cable 2" => "25 mm core · 8 mm XLPE"
        ),
        DOM.p(
            "Curves show the package result and its ±1σ envelope. One debounced solve updates both matrix views.";
            class = "lc-line-parameter-caption"
        ),
        action_row(reset_button, selector_reset_button),
        status_line(status; class = "lc-cable-constants-status"),
        diagnostic(
            "CableDesign → LineCableSystem → LineParametersProblem → compute → observe";
            suffix = " · default formula bundle · total basis"
        );
        kind = :controls,
        class = "lc-cable-constants-controls lc-line-parameters-controls"
    )
end

function self_page(session::Session, state)
    plot = webpart(
        wgl_figure(state.self_plot.figure);
        kind = :plot,
        title = "Self elements",
        meta = "Z[1,1] · Y[1,1] · ±1σ",
        id = "line-parameters-self-plot",
        body_class = "lc-transition-plot-host lc-line-parameters-plot-host"
    )
    controls = controls_panel(
        session,
        state,
        state.self_controls,
        state.self_reset_button,
        state.self_selector_reset_button;
        prefix = "line-self"
    )
    canvas = layout_two_columns(
        plot,
        controls;
        ratio = :wide_narrow,
        height = "100%",
        class = "lc-transition-grid lc-line-parameters-grid"
    )
    return (body = slide(
        "Self line parameters",
        canvas;
        lede = md"""
        Total line parameters of cable 1 from 0.5 Hz to the 50th harmonic. The physical system is rebuilt after the controls settle; persistent Makie plots receive only the new package observations.
        """
    ),)
end

function mutual_page(session::Session, state)
    plot = webpart(
        wgl_figure(state.mutual_plot.figure);
        kind = :plot,
        title = "Mutual elements",
        meta = "Z[1,2] · Y[1,2] · ±1σ",
        id = "line-parameters-mutual-plot",
        body_class = "lc-transition-plot-host lc-line-parameters-plot-host"
    )
    controls = controls_panel(
        session,
        state,
        state.mutual_controls,
        state.mutual_reset_button,
        state.mutual_selector_reset_button;
        prefix = "line-mutual"
    )
    canvas = layout_two_columns(
        plot,
        controls;
        ratio = :wide_narrow,
        height = "100%",
        class = "lc-transition-grid lc-line-parameters-grid"
    )
    return (body = slide(
        "Mutual line parameters",
        canvas;
        lede = md"""
        The same computed result exposes the mutual coupling terms. Controls on this page are synchronized with the self view; changing either page schedules the same latest-only sweep.
        """
    ),)
end

function warm_line_parameters!()
    try
        set_preflight!(
            PREFLIGHT_STATE,
            :model,
            0.12,
            "Building the two-cable uncertain line system…"
        )
        yield()
        model = build_problem(DEFAULT_INPUT; frequencies = FREQUENCIES[1:3])
        model.system.terminal_order == [
            (cable = 1, terminal = :core),
            (cable = 2, terminal = :core)
        ] || throw(ArgumentError("unexpected two-cable terminal order"))

        set_preflight!(
            PREFLIGHT_STATE,
            :solver,
            0.32,
            "Compiling the default line-parameter formulations…"
        )
        calculation = Threads.@spawn compute_line_parameters(DEFAULT_INPUT)
        result = fetch(calculation)
        WARMUP_RESULT[] = result
        yield()

        set_preflight!(
            PREFLIGHT_STATE,
            :curves,
            0.76,
            "Preparing persistent self and mutual plots…"
        )
        WGLMakie.activate!()
        warm_app = App() do session
            state = setup(session)
            return DOM.div(
                self_page(session, state).body...,
                mutual_page(session, state).body...
            )
        end
        parent = Session(NoConnection(); asset_server = NoServer())
        rendered = nothing
        try
            rendered = Bonito.show_html(IOBuffer(), warm_app; parent)
        finally
            isnothing(rendered) || close(rendered)
            close(parent)
        end

        set_preflight!(
            PREFLIGHT_STATE,
            :hot,
            1.0,
            "Line-parameter dashboards hot"
        )
    catch error
        set_preflight!(
            PREFLIGHT_STATE,
            :error,
            PREFLIGHT_STATE[].progress,
            "Line-parameter warmup failed: $(sprint(showerror, error))"
        )
        @error "ICHQP2026 line-parameter warmup failed" exception = (
            error,
            catch_backtrace()
        )
    end
    return nothing
end

function start_warmup!()
    return lock(WARMUP_LOCK) do
        if preflight_ready(PREFLIGHT_STATE)
            return WARMUP_TASK[]
        end
        task = WARMUP_TASK[]
        !isnothing(task) && !istaskdone(task) && return task
        set_preflight!(
            PREFLIGHT_STATE,
            :queued,
            0.02,
            "Starting line-parameter warmup…"
        )
        WARMUP_TASK[] = @async warm_line_parameters!()
        errormonitor(WARMUP_TASK[])
        return WARMUP_TASK[]
    end
end

const PREFLIGHT_RESOURCE = preflight_resource(
    id = "line-parameters-dashboard",
    title = "Two-cable line parameters under uncertainty",
    activate = start_warmup!,
    state = PREFLIGHT_STATE
)

const DECK = deck_descriptor(
    id = "line-parameters",
    group = "ICHQP2026",
    title = "Line parameters",
    order = 40,
    render = true,
    setup = setup,
    resources = (PREFLIGHT_RESOURCE,),
    pages = (
        deck_page(
            "Self line parameters";
            id = "self-line-parameters",
            class = "lc-application-slide lc-fill-page",
            build = self_page
        ),
        deck_page(
            "Mutual line parameters";
            id = "mutual-line-parameters",
            class = "lc-application-slide lc-fill-page",
            build = mutual_page
        )
    )
)

end


ICHQP2026LineParametersDeck.DECK
