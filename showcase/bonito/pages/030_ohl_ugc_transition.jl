module OHLUGCTransitionPage

using Bonito
using Distributed
using GraphMakie
using Graphs
using Makie
using Markdown
using NetworkLayout
using Observables: async_latest, off
using PowerImpedance
using PowerImpedance.NetworkBuilder: Grid, NetworkState
using WGLMakie
using ..PageAuthoring
using ..PowerImpedanceDiagramExt

const DEFAULT_SHARE = 0.50
const DEFAULT_FREQUENCY_RANGE = (100.0, 5.0e3, 400)
const TRANSMISSION_VOLTAGE = 380 / sqrt(3)
const TRANSMISSION_POWER = 100
const EARTH_RESISTIVITY = 100.0
const CORRIDOR_LENGTH = 100.0e3
const CANVAS_COLOR = RGBf(0.9, 0.9, 0.9)
const OHL_COLOR = RGBf(0.796, 0.235, 0.2)
const UGC_COLOR = RGBf(0.251, 0.388, 0.847)
const CORRIDOR_ELEMENTS = (:ohl, :ugc)
const REFERENCE_POWERFLOW = Ref{Any}(nothing)
const REFERENCE_POWERFLOW_GENERATION = Ref(0)
const REFERENCE_POWERFLOW_LOCK = ReentrantLock()
const PREPARED_CASE = Ref{Any}(nothing)
const PREPARATION_TASK = Ref{Union{Nothing, Task}}(nothing)
const PREPARATION_WORKER = Ref{Union{Nothing, Int}}(nothing)
const PREPARATION_LOCK = ReentrantLock()

const CONNECTIONS = (
    (node = :B3d, element = :tl1, side = 2, terminal = 1),
    (node = :B3d, element = :c1, side = 2, terminal = 1),
    (node = :B3q, element = :tl1, side = 2, terminal = 2),
    (node = :B3q, element = :c1, side = 2, terminal = 2),
    (node = :B2d, element = :g4, side = 1, terminal = 1),
    (node = :B2d, element = :tl1, side = 1, terminal = 1),
    (node = :B2q, element = :g4, side = 1, terminal = 2),
    (node = :B2q, element = :tl1, side = 1, terminal = 2),
    (node = :B4, element = :ugc, side = 1, terminal = 1),
    (node = :B4, element = :c1, side = 1, terminal = 1),
    (node = :BX, element = :ugc, side = 2, terminal = 1),
    (node = :BX, element = :ohl, side = 1, terminal = 1),
    (node = :B5, element = :ohl, side = 2, terminal = 1),
    (node = :B5, element = :c2, side = 1, terminal = 1),
    (node = :B6d, element = :tl78, side = 1, terminal = 1),
    (node = :B6d, element = :c2, side = 2, terminal = 1),
    (node = :B6q, element = :tl78, side = 1, terminal = 2),
    (node = :B6q, element = :c2, side = 2, terminal = 2),
    (node = :B7d, element = :g1, side = 1, terminal = 1),
    (node = :B7d, element = :tl78, side = 2, terminal = 1),
    (node = :B7q, element = :g1, side = 1, terminal = 2),
    (node = :B7q, element = :tl78, side = 2, terminal = 2)
)

const BUILDER_OPTIONS = (;
    voltageBase = TRANSMISSION_VOLTAGE,
    power_flow = (; is_bounded = (; bus_voltage = true))
)

function ac_terminal_line(length)
    return overhead_line(
        length = length,
        conductors = Conductors(
            organization = :flat,
            nᵇ = 3,
            nˢᵇ = 1,
            Rᵈᶜ = 0.063,
            rᶜ = 0.015,
            yᵇᶜ = 30,
            Δyᵇᶜ = 0,
            Δxᵇᶜ = 10,
            Δ̃xᵇᶜ = 0,
            dˢᵇ = 0,
            dˢᵃᵍ = 10
        ),
        groundwires = Groundwires(
            nᵍ = 2,
            Rᵍᵈᶜ = 0.92,
            rᵍ = 0.0062,
            Δxᵍ = 6.5,
            Δyᵍ = 7.5,
            dᵍˢᵃᵍ = 10
        ),
        earth_parameters = (1, 1, EARTH_RESISTIVITY),
        transformation = true
    )
end

function converter_one()
    return only(mmc(
        Grid;
        Vᵈᶜ = 640,
        vDCbase = 640,
        Vₘ = TRANSMISSION_VOLTAGE,
        P_max = 1500,
        P_min = -1500,
        P = -TRANSMISSION_POWER,
        Q = 100,
        Q_max = 500,
        Q_min = -500,
        occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
        ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
        pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
        q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159),
        dc = PI_control(Kₚ = 6, Kᵢ = 15),
        timeDelay = 200.0e-6,
        padeOrderNum = 5,
        padeOrderDen = 5
    ))
end

function converter_two()
    return only(mmc(
        Grid;
        Vᵈᶜ = 640,
        vDCbase = 640,
        Vₘ = TRANSMISSION_VOLTAGE,
        P_max = 1000,
        P_min = -1000,
        P = TRANSMISSION_POWER,
        Q = 100,
        Q_max = 1000,
        Q_min = -1000,
        vACbase_LL_RMS = 333,
        turnsRatio = 333 / 380,
        Lᵣ = 0.0461,
        Rᵣ = 0.4103,
        Lₐᵣₘ = 30.0e-3,
        occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
        ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
        pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
        p = PI_control(Kₚ = 0.1, Kᵢ = 31.4159),
        q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159),
        timeDelay = 200.0e-6,
        padeOrderNum = 5,
        padeOrderDen = 5
    ))
end

function corridor_elements(share)
    regularized_share = clamp(share, 1.0e-3, 1 - 1.0e-3)
    ohl = overhead_line(
        length = CORRIDOR_LENGTH * (1 - regularized_share),
        conductors = Conductors(
            organization = :flat,
            nᵇ = 2,
            nˢᵇ = 1,
            Rᵈᶜ = 0.0266,
            rᶜ = 44.8e-3 / 2,
            yᵇᶜ = 18.0,
            Δyᵇᶜ = 0.0,
            Δxᵇᶜ = 7.3,
            Δ̃xᵇᶜ = 0.0,
            dˢᵇ = 0.0,
            dˢᵃᵍ = 6.0
        ),
        groundwires = Groundwires(
            nᵍ = 2,
            Rᵍᵈᶜ = 0.92,
            rᵍ = 0.0062,
            Δxᵍ = 7.3,
            Δyᵍ = 7.0,
            dᵍˢᵃᵍ = 6.0
        ),
        earth_parameters = (1, 1, EARTH_RESISTIVITY),
        transformation = true
    )
    ugc = cable(
        length = CORRIDOR_LENGTH * regularized_share,
        positions = [(-0.5, 1), (0.5, 1)],
        C1 = Conductor(rₒ = 0.02622, ρ = 2.354e-8, μᵣ = 1.035),
        I1 = Insulator(rᵢ = 0.02622, rₒ = 0.06006, ϵᵣ = 2.67, μᵣ = 1.469),
        C2 = Conductor(rᵢ = 0.06006, rₒ = 0.06336, ρ = 2.14e-7, μᵣ = 1.0),
        I2 = Insulator(rᵢ = 0.06336, rₒ = 0.06636, ϵᵣ = 2.3, μᵣ = 1.0),
        C3 = Conductor(rᵢ = 0.06636, rₒ = 0.06651, ρ = 2.826e-8, μᵣ = 1.0),
        I3 = Insulator(rᵢ = 0.06651, rₒ = 0.07256, ϵᵣ = 2.3, μᵣ = 1.0),
        earth_parameters = (1, 1, EARTH_RESISTIVITY),
        transformation = true
    )
    return (; ohl, ugc)
end

function normalized_share(share::Real)
    isfinite(share) || throw(DomainError(share, "UGC share must be finite"))
    return clamp(Float64(share), 0.0, 1.0)
end

function transition_network(share::Real)
    share = normalized_share(share)
    corridor = corridor_elements(share)
    elements = (;
        g1 = ac_source(
            V = TRANSMISSION_VOLTAGE,
            P = TRANSMISSION_POWER,
            P_min = -2000,
            P_max = 2000,
            Q_max = 1000,
            Q_min = -1000,
            pins = 3,
            transformation = true
        ),
        c1 = converter_one(),
        c2 = converter_two(),
        corridor.ugc,
        corridor.ohl,
        g4 = ac_source(
            V = TRANSMISSION_VOLTAGE,
            P = TRANSMISSION_POWER,
            P_min = -2000,
            P_max = 2000,
            Q_max = 1000,
            Q_min = -1000,
            pins = 3,
            transformation = true
        ),
        tl1 = ac_terminal_line(25.0e3),
        tl78 = ac_terminal_line(90.0e3)
    )
    return PowerImpedance.NetworkBuilder.define(
        elements,
        CONNECTIONS;
        options = BUILDER_OPTIONS
    )
end

function reference_powerflow(network::NetworkState = transition_network(DEFAULT_SHARE))
    return PowerImpedance.compute(
        PowerImpedance.PowerFlowProblem(network),
        PowerImpedance.ACDCPowerFlow()
    )
end

function cached_reference_powerflow(; force::Bool = false)
    return lock(REFERENCE_POWERFLOW_LOCK) do
        if force || isnothing(REFERENCE_POWERFLOW[])
            result = reference_powerflow()
            REFERENCE_POWERFLOW[] = result
            REFERENCE_POWERFLOW_GENERATION[] += 1
        end
        return (;
            powerflow = REFERENCE_POWERFLOW[],
            generation = REFERENCE_POWERFLOW_GENERATION[]
        )
    end
end

function linearized_network(network::NetworkState, powerflow)
    result = PowerImpedance.compute(
        PowerImpedance.LinearizationProblem(network, powerflow),
        PowerImpedance.AdmittanceLinearization()
    )
    return result.network_model
end

function harmonic_response(
        network_model;
        frequency_range = DEFAULT_FREQUENCY_RANGE
)
    problem = PowerImpedance.PowerImpedanceProblem(
        network_model;
        nodes = [:B5],
        eliminated_elements = [:c2],
        frequency_range
    )
    return PowerImpedance.compute(problem, PowerImpedance.NodalImpedance())
end

function report_preparation!(updates, phase, progress, label)
    state = (; phase, progress, label)
    if isnothing(updates)
        set_preparation!(phase, progress, label)
    else
        put!(updates, state)
    end
    return state
end

function compute_prepared_case(updates = nothing)
    report_preparation!(
        updates,
        :network,
        0.1,
        "Constructing the reference network…"
    )
    network = transition_network(DEFAULT_SHARE)

    report_preparation!(
        updates,
        :powerflow,
        0.25,
        "Solving the reference power flow…"
    )
    reference = cached_reference_powerflow()

    report_preparation!(
        updates,
        :linearization,
        0.55,
        "Linearizing the network model…"
    )
    network_model = linearized_network(network, reference.powerflow)

    report_preparation!(
        updates,
        :response,
        0.82,
        "Computing the initial B5 impedance response…"
    )
    response = harmonic_response(network_model)
    curve = harmonic_curve(response)

    return (;
        network,
        network_model,
        powerflow = reference.powerflow,
        generation = reference.generation,
        response,
        curve
    )
end

function preparation_error!(error)
    set_preparation!(
        :error,
        1.0,
        "Preparation failed: $(sprint(showerror, error))"
    )
    return error
end

function prepare_case!()
    try
        return lock(PREPARATION_LOCK) do
            !isnothing(PREPARED_CASE[]) && return PREPARED_CASE[]
            PREPARED_CASE[] = compute_prepared_case()
            set_preparation!(:ready, 1.0, "Application case ready")
            return PREPARED_CASE[]
        end
    catch error
        preparation_error!(error)
        rethrow()
    end
end

is_prepared() = !isnothing(PREPARED_CASE[])

function start_preparation!()
    is_prepared() && return PREPARATION_TASK[]
    task = PREPARATION_TASK[]
    !isnothing(task) && !istaskdone(task) && return task
    set_preparation!(:queued, 0.02, "Starting numerical preparation…")

    PREPARATION_TASK[] = @async begin
        worker = nothing
        try
            project = Base.active_project()
            flags = ["--project=$(dirname(project))", "--startup-file=no"]
            worker = only(addprocs(1; exeflags = flags))
            PREPARATION_WORKER[] = worker

            imports = quote
                using LineCableModels
                import PowerImpedance
                using Bonito
                using GraphMakie
                using Graphs
                using NetworkLayout
                using WGLMakie
                WGLMakie.activate!()
            end
            remotecall_wait(Core.eval, worker, Main, imports)
            remotecall_wait(
                Base.include,
                worker,
                Main,
                abspath(joinpath(@__DIR__, "..", "app.jl"))
            )

            updates = RemoteChannel(() -> Channel{NamedTuple}(16), 1)
            future = remotecall(compute_prepared_case, worker, updates)
            while !isready(future)
                while isready(updates)
                    state = take!(updates)
                    set_preparation!(state.phase, state.progress, state.label)
                end
                sleep(0.1)
            end
            while isready(updates)
                state = take!(updates)
                set_preparation!(state.phase, state.progress, state.label)
            end
            PREPARED_CASE[] = fetch(future)
            set_preparation!(:ready, 1.0, "Application case ready")
        catch error
            preparation_error!(error)
            @error "OHL/UGC application-case preparation failed" exception = (
                error,
                catch_backtrace()
            )
        finally
            if !isnothing(worker) && worker in workers()
                try
                    rmprocs(worker)
                catch error
                    @warn "Could not stop the numerical-preparation worker" exception = (
                        error,
                        catch_backtrace()
                    )
                end
            end
            PREPARATION_WORKER[] = nothing
        end
    end
    errormonitor(PREPARATION_TASK[])
    return PREPARATION_TASK[]
end

prepared_case() = is_prepared() ? PREPARED_CASE[] : prepare_case!()

function set_corridor_lengths!(network::NetworkState, share::Real)
    share = normalized_share(share)
    regularized_share = clamp(share, 1.0e-3, 1 - 1.0e-3)
    network.elements.ohl.element_model.length = CORRIDOR_LENGTH * (1 - regularized_share)
    network.elements.ugc.element_model.length = CORRIDOR_LENGTH * regularized_share
    return share
end

function harmonic_curve(response; node = :B5)
    nodes = collect(PowerImpedance.response_nodes(response))
    node_index = only(findall(==(node), nodes))
    frequencies_hz = real.(PowerImpedance.angular_frequencies(response)) ./ (2π)
    impedance = @view PowerImpedance.response_values(response)[
        node_index,
        node_index,
        :
    ]
    return (;
        frequencies_hz,
        magnitude_db = 20 .* log10.(abs.(vec(impedance)))
    )
end

function transition_color(share::Real)
    amount = clamp(Float32(share), 0.0f0, 1.0f0)
    return RGBf(
        (1 - amount) * OHL_COLOR.r + amount * UGC_COLOR.r,
        (1 - amount) * OHL_COLOR.g + amount * UGC_COLOR.g,
        (1 - amount) * OHL_COLOR.b + amount * UGC_COLOR.b
    )
end

function corridor_edge_indices(diagram)
    indices = Int[]
    for (index, edge) in enumerate(Graphs.edges(diagram.model.graph))
        left = diagram.model.vertex_keys[Graphs.src(edge)]
        right = diagram.model.vertex_keys[Graphs.dst(edge)]
        component = hasproperty(left, :element) ? left : right
        component.element in CORRIDOR_ELEMENTS && push!(indices, index)
    end
    return indices
end

function color_corridor!(diagram, share::Real)
    color = transition_color(share)
    corridor_edges = corridor_edge_indices(diagram)
    edge_colors = copy(diagram.plots.graph.edge_color[])
    edge_colors[corridor_edges] .= color
    diagram.plots.graph.edge_color[] = edge_colors
    edge_widths = copy(diagram.plots.graph.edge_width[])
    edge_widths[corridor_edges] .= 6.0f0
    diagram.plots.graph.edge_width[] = edge_widths
    for (component, plot) in zip(diagram.model.components, diagram.plots.components)
        component.key.element in CORRIDOR_ELEMENTS || continue
        plot.color[] = color
        plot.strokecolor[] = color
    end
    return color
end

function use_node_labels!(diagram)
    for (bus, label) in zip(diagram.model.buses, diagram.plots.bus_labels)
        label.text[] = [join(string.(bus.nodes), "/")]
    end
    return diagram
end

mutable struct TransitionRuntime
    network::NetworkState
    network_model::Any
    reference_powerflow::Any
    reference_generation::Int
    impedance_plot::Any
    impedance_line::Any
    diagram::Any
    cache::Dict{Int, NamedTuple}
    share::Float64
    lock::ReentrantLock
end

share_key(share::Real) = round(Int, 100normalized_share(share))

function impedance_figure(response)
    curve = harmonic_curve(response)
    figure = Figure(
        size = (960, 420),
        backgroundcolor = CANVAS_COLOR
    )
    axis = Axis(
        figure[1, 1];
        xlabel = "Frequency (Hz)",
        ylabel = "|Z| (dBΩ)",
        xscale = log10,
        backgroundcolor = CANVAS_COLOR,
        xgridcolor = RGBf(0.74, 0.74, 0.74),
        ygridcolor = RGBf(0.74, 0.74, 0.74),
        xlabelsize = 16,
        ylabelsize = 16,
        xticklabelsize = 13,
        yticklabelsize = 13
    )
    line = lines!(
        axis,
        curve.frequencies_hz,
        curve.magnitude_db;
        color = UGC_COLOR,
        linewidth = 3
    )
    xlims!(axis, extrema(curve.frequencies_hz))
    return (; figure, axis, line)
end

function build_runtime(;
        initial_share::Real = DEFAULT_SHARE,
        frequency_range = DEFAULT_FREQUENCY_RANGE
)
    prepared = prepared_case()
    network, network_model = deepcopy((prepared.network, prepared.network_model))
    set_corridor_lengths!(network, initial_share)
    response = if initial_share == DEFAULT_SHARE &&
                  frequency_range == DEFAULT_FREQUENCY_RANGE
        prepared.response
    else
        harmonic_response(network_model; frequency_range)
    end
    curve = initial_share == DEFAULT_SHARE && frequency_range == DEFAULT_FREQUENCY_RANGE ?
            prepared.curve : harmonic_curve(response)
    impedance_plot = impedance_figure(response)
    impedance_line = impedance_plot.line

    diagram = PowerImpedanceDiagramExt.diagram(
        network,
        prepared.powerflow;
        layout = NetworkLayout.Stress(; seed = 7),
        interactive = false,
        style = (;
            figure_size = (960, 300),
            background = CANVAS_COLOR,
            line_color = RGBf(0.35, 0.35, 0.35),
            label_color = RGBf(0.12, 0.12, 0.12)
        )
    )
    use_node_labels!(diagram)
    diagram.axis.aspect[] = nothing
    diagram.axis.backgroundcolor[] = CANVAS_COLOR
    color_corridor!(diagram, initial_share)

    cache = Dict{Int, NamedTuple}(share_key(initial_share) => curve)
    return TransitionRuntime(
        network,
        network_model,
        prepared.powerflow,
        prepared.generation,
        impedance_plot,
        impedance_line,
        diagram,
        cache,
        Float64(initial_share),
        ReentrantLock()
    )
end

function display_curve!(runtime::TransitionRuntime, curve)
    Makie.update!(runtime.impedance_line; arg2 = curve.magnitude_db)
    Makie.autolimits!(runtime.impedance_plot.axis)
    Makie.xlims!(runtime.impedance_plot.axis, extrema(curve.frequencies_hz))
    return curve
end

function update_impedance!(
        runtime::TransitionRuntime,
        share::Real;
        frequency_range = DEFAULT_FREQUENCY_RANGE
)
    key = share_key(share)
    normalized_share = key / 100
    return lock(runtime.lock) do
        set_corridor_lengths!(runtime.network, normalized_share)
        curve = get!(runtime.cache, key) do
            response = harmonic_response(runtime.network_model; frequency_range)
            harmonic_curve(response)
        end
        display_curve!(runtime, curve)
        runtime.share = normalized_share
        return curve
    end
end

function update!(
        runtime::TransitionRuntime,
        share::Real;
        frequency_range = DEFAULT_FREQUENCY_RANGE
)
    normalized_share = share_key(share) / 100
    color_corridor!(runtime.diagram, normalized_share)
    return update_impedance!(runtime, normalized_share; frequency_range)
end

function recache!(
        runtime::TransitionRuntime,
        share::Real = runtime.share;
        frequency_range = DEFAULT_FREQUENCY_RANGE,
        reference = cached_reference_powerflow(; force = true)
)
    normalized_share = share_key(share) / 100
    return lock(runtime.lock) do
        set_corridor_lengths!(runtime.network, normalized_share)
        network_model = linearized_network(runtime.network, reference.powerflow)
        response = harmonic_response(network_model; frequency_range)
        curve = harmonic_curve(response)

        runtime.network_model = network_model
        runtime.reference_powerflow = reference.powerflow
        runtime.reference_generation = reference.generation
        empty!(runtime.cache)
        runtime.cache[share_key(normalized_share)] = curve
        display_curve!(runtime, curve)
        color_corridor!(runtime.diagram, normalized_share)
        runtime.share = normalized_share
        return curve
    end
end

function content(session::Session)
    slider = Bonito.Slider(
        collect(0.0:0.01:1.0);
        value = DEFAULT_SHARE,
        id = "ugc-share-input",
        ariaLabel = "Underground cable share"
    )
    recache_button = Bonito.Button(
        "Re-cache power flow + linearization";
        style = nothing,
        id = "recache-powerflow-control",
        class = "lc-recache-button",
        type = "button",
        ariaLabel = "Recompute the reference power flow and linearization"
    )
    runtime = build_runtime()
    requested_share = slider.value
    color_observer = on(requested_share) do share
        color_corridor!(runtime.diagram, normalized_share(share))
        return nothing
    end
    latest_share = async_latest(requested_share, 1)
    status = Observable(
        "Operating point #$(runtime.reference_generation) · response ready"
    )
    share_label = map(
        share -> "$(round(Int, 100share))% UGC",
        session,
        requested_share
    )
    ugc_length_label = map(
        share -> "$(round(Int, CORRIDOR_LENGTH * share / 1.0e3)) km",
        session,
        requested_share
    )
    ohl_length_label = map(
        share -> "$(round(Int, CORRIDOR_LENGTH * (1 - share) / 1.0e3)) km",
        session,
        requested_share
    )

    response_observer = on(latest_share) do share
        share_key(share) == share_key(runtime.share) && return nothing
        cached = haskey(runtime.cache, share_key(share))
        status[] = cached ? "Loading cached response…" : "Updating B5 impedance…"
        try
            update_impedance!(runtime, share)
            status[] = cached ?
                       "Cached response · operating point #$(runtime.reference_generation)" :
                       "Response ready · operating point #$(runtime.reference_generation)"
        catch error
            status[] = "Update failed: $(sprint(showerror, error))"
            @error "OHL/UGC transition update failed" exception = (
                error,
                catch_backtrace()
            )
        end
        return nothing
    end

    recache_running = Ref(false)
    recache_observer = on(recache_button.value) do clicked
        clicked || return nothing
        recache_running[] && return nothing
        recache_running[] = true
        status[] = "Recomputing power flow + linearization…"
        @async begin
            try
                reference = cached_reference_powerflow(; force = true)
                recache!(runtime, requested_share[]; reference)
                status[] = "Operating point #$(runtime.reference_generation) re-cached · response ready"
            catch error
                status[] = "Re-cache failed: $(sprint(showerror, error))"
                @error "OHL/UGC transition re-cache failed" exception = (
                    error,
                    catch_backtrace()
                )
            finally
                recache_running[] = false
            end
        end
        return nothing
    end

    on(session.on_close) do _
        off(color_observer)
        off(response_observer)
        off(recache_observer)
        return nothing
    end

    controls = webpart(
        control(
            "Underground-cable share",
            slider;
            value = share_label,
            id = "ugc-share-control",
            class = "lc-transition-control"
        ),
        color_key(
            "OHL" => "is-ohl",
            "UGC" => "is-ugc";
            class = "lc-corridor-key"
        ),
        value_list(
            "Overhead line" => ohl_length_label,
            "Underground cable" => ugc_length_label,
            "Observation node" => "B5",
            "Frequency sweep" => "100 Hz–5 kHz";
            class = "lc-transition-values"
        ),
        recache_button,
        status_line(status; class = "lc-transition-status"),
        diagnostic(
            "PowerImpedance.jl";
            suffix = " · fixed operating point",
            class = "lc-transition-diagnostic"
        );
        kind = :controls,
        class = "lc-transition-controls"
    )
    diagram = webpart(
        WGLMakie.WithConfig(runtime.diagram.figure; resize_to = :parent);
        kind = :plot,
        title = "Network diagram",
        meta = "B4–B5 corridor",
        id = "transition-network-diagram",
        body_class = "lc-transition-plot-host"
    )
    impedance = webpart(
        WGLMakie.WithConfig(runtime.impedance_plot.figure; resize_to = :parent);
        kind = :plot,
        title = "Driving-point impedance",
        meta = "Z[B5, B5] · 100 Hz–5 kHz",
        id = "b5-impedance-plot",
        body_class = "lc-transition-plot-host"
    )
    body = webgrid(
        [:controls :diagram; :controls :impedance];
        columns = "minmax(16rem, 3fr) minmax(0, 7fr)",
        compact_columns = "minmax(15rem, 2fr) minmax(0, 5fr)",
        rows = "minmax(0, 2fr) minmax(0, 3fr)",
        gap = "0.75rem 1rem",
        height = "auto",
        stack_rows = "auto 17rem 26rem",
        class = "lc-interactive-grid lc-transition-grid",
        controls,
        diagram,
        impedance
    )
    return (;
        body,
        runtime,
        slider,
        recache_button,
        status,
        observer = response_observer,
        color_observer,
        recache_observer
    )
end

function build(session::Session)
    page = content(session)
    return merge(page,
        (;
            body = slide(
            "Application case - parametric OHL/UGC transition",
            page.body;
            lede = md"""
            Move the corridor share from overhead line to underground cable and inspect the driving-point impedance ``Z_{B_5,B_5}(f)``.
            """
        )
        ))
end

const PAGE = page_descriptor(
    id = "parametric-ohl-ugc-transition",
    group = "Application cases",
    title = "Parametric OHL/UGC transition",
    order = 30,
    render = true,
    class = "lc-application-slide",
    build = build,
    prepare = start_preparation!,
    ready = is_prepared
)

end

OHLUGCTransitionPage.PAGE
