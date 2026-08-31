isdefined(@__MODULE__, :OHLUGCTransitionPage) || include(
    joinpath(
        @__DIR__,
        "..",
        "..",
        "live-manual",
        "content",
        "020_ohl_ugc_transition.jl"
    )
)

module ICHQP2026ApplicationCaseDeck

using Bonito
using Makie
using Markdown
using Observables: async_latest, off
using Printf
using WGLMakie
import ..OHLUGCTransitionPage
using ..PageAuthoring

const Case = OHLUGCTransitionPage
const SHARE_RANGE = 0.0:0.01:1.0
const ERROR_RANGE_PCT = 0.0:0.25:10.0
const DEFAULT_INPUT = (
    share = Case.DEFAULT_SHARE,
    error_pct = 5.0
)
const LOWER_COLOR = RGBf(0.80, 0.24, 0.20)
const NOMINAL_COLOR = RGBf(0.13, 0.38, 0.84)
const UPPER_COLOR = RGBf(0.68, 0.27, 0.72)
const CURVE_KEYS = (:lower, :nominal, :upper)
const CURVE_COLORS = (
    lower = LOWER_COLOR,
    nominal = NOMINAL_COLOR,
    upper = UPPER_COLOR
)
const PREFLIGHT_STATE = preflight_state(
    label = "Application-case power flow not activated"
)
const WARMUP_TASK = Ref{Union{Nothing, Task}}(nothing)
const WARMUP_BUNDLE = Ref{Any}(nothing)
const WARMUP_LOCK = ReentrantLock()

function error_percent(number::Real)
    isfinite(number) && zero(number) <= number <= oftype(number, 100) ||
        throw(DomainError(number, "length error must lie between 0 and 100 percent"))
    return Float64(number)
end

"""Return deterministic lower, nominal, and upper UGC transition shares."""
function transition_bounds(share::Real, error_pct::Real)
    nominal = Case.normalized_share(share)
    error = error_percent(error_pct)
    offset = nominal * error / 100
    return (
        lower = clamp(nominal - offset, 0.0, 1.0),
        nominal,
        upper = clamp(nominal + offset, 0.0, 1.0)
    )
end

curve_key(share::Real) = round(Int, 1.0e6 * Case.normalized_share(share))

function curve_at!(network, network_model, cache, share::Real)
    normalized = Case.normalized_share(share)
    return get!(cache, curve_key(normalized)) do
        Case.set_corridor_lengths!(network, normalized)
        Case.harmonic_curve(Case.harmonic_response(network_model))
    end
end

function curve_bundle!(network, network_model, cache, input)
    bounds = transition_bounds(input.share, input.error_pct)
    curves = map(CURVE_KEYS) do key
        curve_at!(network, network_model, cache, getproperty(bounds, key))
    end
    Case.set_corridor_lengths!(network, input.share)
    return (; bounds, curves = NamedTuple{CURVE_KEYS}(curves))
end

function curve_title(key::Symbol, share::Real, error_pct::Real)
    length_km = Case.CORRIDOR_LENGTH * share / 1.0e3
    qualifier = key === :lower ? "Base − $(error_pct)%" :
                key === :upper ? "Base + $(error_pct)%" : "Base"
    return @sprintf("%s · %.1f km UGC", qualifier, length_km)
end

function impedance_bounds_figure(bundle, error_pct::Real)
    figure = Figure(size = (780, 860), backgroundcolor = Case.CANVAS_COLOR)
    axes = Any[]
    lines = Any[]
    for (row, key) in enumerate(CURVE_KEYS)
        curve = getproperty(bundle.curves, key)
        axis = Axis(
            figure[row, 1];
            title = curve_title(
                key,
                getproperty(bundle.bounds, key),
                error_pct
            ),
            xlabel = row == length(CURVE_KEYS) ? "Frequency (Hz)" : "",
            ylabel = "|Z| (dBΩ)",
            xscale = log10,
            backgroundcolor = Case.CANVAS_COLOR,
            xgridcolor = RGBf(0.74, 0.74, 0.74),
            ygridcolor = RGBf(0.74, 0.74, 0.74),
            titlesize = 14,
            xlabelsize = 14,
            ylabelsize = 14,
            xticklabelsize = 11,
            yticklabelsize = 11
        )
        line = lines!(
            axis,
            curve.frequencies_hz,
            curve.magnitude_db;
            color = getproperty(CURVE_COLORS, key),
            linewidth = 2.8
        )
        row < length(CURVE_KEYS) && hidexdecorations!(axis; grid = false)
        push!(axes, axis)
        push!(lines, line)
    end
    linkxaxes!(axes...)
    rowgap!(figure.layout, 5)
    plot = (;
        figure,
        axes = NamedTuple{CURVE_KEYS}(Tuple(axes)),
        lines = NamedTuple{CURVE_KEYS}(Tuple(lines))
    )
    update_impedance_bounds!(plot, bundle, error_pct)
    return plot
end

function update_impedance_bounds!(plot, bundle, error_pct::Real)
    magnitudes = Float64[]
    for key in CURVE_KEYS
        curve = getproperty(bundle.curves, key)
        axis = getproperty(plot.axes, key)
        line = getproperty(plot.lines, key)
        Makie.update!(line; arg2 = curve.magnitude_db)
        axis.title[] = curve_title(
            key,
            getproperty(bundle.bounds, key),
            error_pct
        )
        xlims!(axis, extrema(curve.frequencies_hz))
        append!(magnitudes, curve.magnitude_db)
    end
    lower, upper = extrema(magnitudes)
    padding = max(1.0, 0.06 * (upper - lower))
    for axis in values(plot.axes)
        ylims!(axis, lower - padding, upper + padding)
    end
    return bundle
end

mutable struct ApplicationRuntime
    case_runtime::Any
    impedance_plot::Any
    cache::Dict{Int, NamedTuple}
    input::NamedTuple
    lock::ReentrantLock
end

function initial_bundle!(case_runtime, input)
    cache = Dict{Int, NamedTuple}()
    cache[curve_key(case_runtime.share)] =
        case_runtime.cache[Case.share_key(case_runtime.share)]
    warmed = WARMUP_BUNDLE[]
    if !isnothing(warmed) && warmed.input == input
        for key in CURVE_KEYS
            share = getproperty(warmed.bundle.bounds, key)
            cache[curve_key(share)] = getproperty(warmed.bundle.curves, key)
        end
        return warmed.bundle, cache
    end
    return curve_bundle!(
        case_runtime.network,
        case_runtime.network_model,
        cache,
        input
    ), cache
end

function build_runtime(input = DEFAULT_INPUT)
    case_runtime = Case.build_runtime(initial_share = input.share)
    Case.update_corridor_diagram!(case_runtime.diagram, input.share)
    bundle, cache = initial_bundle!(case_runtime, input)
    plot = impedance_bounds_figure(bundle, input.error_pct)
    return ApplicationRuntime(
        case_runtime,
        plot,
        cache,
        input,
        ReentrantLock()
    )
end

function update!(runtime::ApplicationRuntime, input)
    return lock(runtime.lock) do
        bundle = curve_bundle!(
            runtime.case_runtime.network,
            runtime.case_runtime.network_model,
            runtime.cache,
            input
        )
        update_impedance_bounds!(runtime.impedance_plot, bundle, input.error_pct)
        runtime.input = input
        return bundle
    end
end

function setup(session::Session)
    share_slider = Bonito.Slider(
        collect(SHARE_RANGE);
        value = DEFAULT_INPUT.share,
        id = "application-transition-length-input",
        ariaLabel = "Underground cable share of the B4 to B5 corridor"
    )
    error_slider = Bonito.Slider(
        collect(ERROR_RANGE_PCT);
        value = DEFAULT_INPUT.error_pct,
        id = "application-transition-error-input",
        ariaLabel = "Deterministic transition length error in percent"
    )
    input = map(
        (share, error_pct) -> (;
            share = Float64(share),
            error_pct = Float64(error_pct)
        ),
        session,
        share_slider.value,
        error_slider.value
    )
    runtime = build_runtime(input[])
    status = Observable("Three deterministic harmonic responses ready")

    share_label = map(
        state -> @sprintf("%.1f km UGC", Case.CORRIDOR_LENGTH * state.share / 1.0e3),
        session,
        input
    )
    error_label = map(
        state -> @sprintf("±%.2g %%", state.error_pct),
        session,
        input
    )
    lower_label = map(
        state -> @sprintf(
            "%.1f km",
            Case.CORRIDOR_LENGTH * transition_bounds(
                state.share,
                state.error_pct
            ).lower / 1.0e3
        ),
        session,
        input
    )
    nominal_label = map(
        state -> @sprintf("%.1f km", Case.CORRIDOR_LENGTH * state.share / 1.0e3),
        session,
        input
    )
    upper_label = map(
        state -> @sprintf(
            "%.1f km",
            Case.CORRIDOR_LENGTH * transition_bounds(
                state.share,
                state.error_pct
            ).upper / 1.0e3
        ),
        session,
        input
    )

    diagram_observer = on(input) do state
        Case.update_corridor_diagram!(runtime.case_runtime.diagram, state.share)
        return nothing
    end
    latest_input = async_latest(input, 1)
    response_observer = on(latest_input) do state
        state == runtime.input && return nothing
        status[] = "Updating deterministic length bounds…"
        started = time_ns()
        try
            update!(runtime, state)
            elapsed_ms = (time_ns() - started) / 1.0e6
            status[] = @sprintf(
                "Three harmonic responses ready · %.0f ms",
                elapsed_ms
            )
        catch error
            status[] = "Update failed: $(sprint(showerror, error))"
            @error "ICHQP2026 application-case update failed" exception = (
                error,
                catch_backtrace()
            )
        end
        return nothing
    end
    on(session.on_close) do _
        off(diagram_observer)
        off(response_observer)
        return nothing
    end

    controls = webpart(
        control(
            "UGC transition length",
            share_slider;
            value = share_label,
            id = "application-transition-length-control"
        ),
        control(
            "Deterministic length error",
            error_slider;
            value = error_label,
            id = "application-transition-error-control"
        ),
        color_key(
            "OHL segment" => "is-ohl",
            "UGC segment" => "is-ugc";
            class = "lc-corridor-key"
        ),
        value_list(
            "Base − error" => lower_label,
            "Base" => nominal_label,
            "Base + error" => upper_label,
            "Observation node" => "B5",
            "Frequency sweep" => "100 Hz–5 kHz"
        ),
        status_line(status; class = "lc-transition-status"),
        diagnostic(
            "Float64 lower / nominal / upper";
            suffix = " · one frozen power-flow linearization"
        );
        kind = :controls,
        class = "lc-application-case-controls"
    )
    diagram = webpart(
        wgl_figure(runtime.case_runtime.diagram.figure);
        kind = :plot,
        title = "Network diagram",
        meta = "B4–B5 transition location",
        id = "application-network-diagram",
        body_class = "lc-transition-plot-host lc-application-case-diagram-host"
    )
    impedance = webpart(
        wgl_figure(runtime.impedance_plot.figure);
        kind = :plot,
        title = "B5 harmonic impedance",
        meta = "deterministic length bounds",
        id = "application-impedance-bounds",
        body_class = "lc-transition-plot-host lc-application-case-impedance-host"
    )
    return (;
        runtime,
        input,
        share_slider,
        error_slider,
        status,
        controls,
        diagram,
        impedance,
        diagram_observer,
        response_observer
    )
end

function application_page(::Session, state)
    side = webgrid(
        reshape([:controls, :diagram], 2, 1);
        rows = "minmax(13rem, 2fr) minmax(0, 3fr)",
        height = "100%",
        stack = false,
        class = "lc-application-case-side",
        controls = state.controls,
        diagram = state.diagram
    )
    canvas = webgrid(
        [:impedance :side];
        columns = "minmax(0, 11fr) minmax(24rem, 9fr)",
        compact_columns = "minmax(0, 1fr) minmax(21rem, 1fr)",
        rows = "minmax(0, 1fr)",
        height = "100%",
        stack_rows = "42rem 36rem",
        class = "lc-application-case-grid",
        impedance = state.impedance,
        side
    )
    return (body = slide(
        "Application case – parametric OHL/UGC transition",
        canvas;
        lede = md"""
        Compare deterministic lower, nominal, and upper transition-length
        boundaries at ``B_5``. No probabilistic number type enters this run;
        the cached PowerImpedance operating point remains fixed.
        """
    ),)
end

function warm_application_case!()
    try
        set_preflight!(
            PREFLIGHT_STATE,
            :network,
            0.05,
            "Preparing the PowerImpedance operating point…"
        )
        preparation = Case.start_preparation!()
        while !isnothing(preparation) && !istaskdone(preparation)
            source = Case.PREFLIGHT_STATE[]
            set_preflight!(
                PREFLIGHT_STATE,
                source.phase,
                min(0.76, 0.76 * source.progress),
                source.label
            )
            sleep(0.1)
        end
        isnothing(preparation) || wait(preparation)
        Case.is_prepared() || throw(ErrorException(Case.PREFLIGHT_STATE[].label))

        set_preflight!(
            PREFLIGHT_STATE,
            :bounds,
            0.80,
            "Computing deterministic transition-length bounds…"
        )
        prepared = Case.prepared_case()
        network = deepcopy(prepared.network)
        network_model = Case.linearized_network(network, prepared.powerflow)
        cache = Dict{Int, NamedTuple}(
            curve_key(DEFAULT_INPUT.share) => prepared.curve
        )
        calculation = Threads.@spawn curve_bundle!(
            network,
            network_model,
            cache,
            DEFAULT_INPUT
        )
        bundle = fetch(calculation)
        WARMUP_BUNDLE[] = (; input = DEFAULT_INPUT, bundle)

        set_preflight!(
            PREFLIGHT_STATE,
            :render,
            0.92,
            "Preparing persistent plots and network glyphs…"
        )
        WGLMakie.activate!()
        warm_app = App() do session
            state = setup(session)
            return DOM.div(application_page(session, state).body...)
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
            "Application-case power flow and deterministic bounds hot"
        )
    catch error
        set_preflight!(
            PREFLIGHT_STATE,
            :error,
            PREFLIGHT_STATE[].progress,
            "Application-case warmup failed: $(sprint(showerror, error))"
        )
        @error "ICHQP2026 application-case warmup failed" exception = (
            error,
            catch_backtrace()
        )
    end
    return nothing
end

function start_warmup!()
    return lock(WARMUP_LOCK) do
        preflight_ready(PREFLIGHT_STATE) && return WARMUP_TASK[]
        task = WARMUP_TASK[]
        !isnothing(task) && !istaskdone(task) && return task
        set_preflight!(
            PREFLIGHT_STATE,
            :queued,
            0.02,
            "Starting application-case warmup…"
        )
        WARMUP_TASK[] = @async warm_application_case!()
        errormonitor(WARMUP_TASK[])
        return WARMUP_TASK[]
    end
end

const PREFLIGHT_RESOURCE = preflight_resource(
    id = "ichqp2026-application-case",
    title = "Parametric OHL/UGC harmonic impedance",
    activate = start_warmup!,
    state = PREFLIGHT_STATE
)

const DECK = deck_descriptor(
    id = "application-case",
    group = "ICHQP2026",
    title = "Application case",
    order = 50,
    render = true,
    setup = setup,
    resources = (PREFLIGHT_RESOURCE,),
    pages = (
        deck_page(
            "Parametric OHL/UGC transition";
            id = "parametric-ohl-ugc-transition",
            class = "lc-application-slide lc-fill-page",
            build = application_page
        ),
    )
)

end

ICHQP2026ApplicationCaseDeck.DECK
