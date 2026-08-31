module ICHQP2026CableConstantsDeck

using Bonito
import LineCableModels as LCM
using LinearAlgebra: norm
using Makie
using Markdown
using Measurements: Measurement, measurement, uncertainty, value
using Observables: off
using Printf
using WGLMakie
using ..PageAuthoring

const STRAND_LAYER_RANGE = 1:1:10
const FREQUENCY_LOG_RANGE = 1.0:0.05:6.0
const WIRE_DIAMETER_RANGE_MM = 2.0:0.1:8.0
const UNCERTAINTY_RANGE_PCT = 0.0:0.1:10.0
const INSULATION_THICKNESS_RANGE_MM = 2.0:0.1:20.0
const PREVIEW_MARGIN = 1.08
const CLONE_LATERAL_OFFSET_MM = 1.0e3
const ENVELOPE_COVERAGE = 1.0
const ENVELOPE_ALPHA = 0.22
const WIRE_ENVELOPE_COLOR = RGBAf(0.10, 0.78, 0.88, ENVELOPE_ALPHA)
const INSULATION_ENVELOPE_COLOR = RGBAf(0.92, 0.24, 0.72, ENVELOPE_ALPHA)
const BOUNDARY_SAMPLES = 128
const DEFAULT_INPUT = (
    frequency = 1.0e3,
    strand_layers = 2,
    wire_diameter = 4.7e-3,
    wire_uncertainty = 1.0,
    insulation_thickness = 8.0e-3,
    insulation_uncertainty = 1.0
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
const CABLE_CONSTANTS_FORMULATION = LCM.CableConstantsFormulation()
const PLOT_BACKGROUND = RGBf(0.9, 0.9, 0.9)
const PLOT_TEXT = RGBf(0.16, 0.18, 0.18)
const RESULT_UNITS = (
    R = LCM.display_unit(LCM.R, :pul),
    L = LCM.display_unit(LCM.L, :pul),
    C = LCM.display_unit(LCM.C, :pul),
    G = LCM.display_unit(LCM.G, :pul; prefix = :micro)
)
const RESULT_FACTORS = (
    R = LCM.scale_factor(LCM.R, :pul, RESULT_UNITS.R),
    L = LCM.scale_factor(LCM.L, :pul, RESULT_UNITS.L),
    C = LCM.scale_factor(LCM.C, :pul, RESULT_UNITS.C),
    G = LCM.scale_factor(LCM.G, :pul, RESULT_UNITS.G)
)
const RESULT_UNIT_LABELS = (
    R = LCM.label(RESULT_UNITS.R),
    L = LCM.label(RESULT_UNITS.L),
    C = LCM.label(RESULT_UNITS.C),
    G = LCM.label(RESULT_UNITS.G)
)
const PREFLIGHT_STATE = preflight_state(
    label = "Cable-constants dashboard not activated"
)
const WARMUP_TASK = Ref{Union{Nothing, Task}}(nothing)
const WARMUP_RESULT = Ref{Any}(nothing)
const WARMUP_LOCK = ReentrantLock()
const SI_PREFIX = Dict(
    -12 => "p",
    -9 => "n",
    -6 => "μ",
    -3 => "m",
    0 => "",
    3 => "k",
    6 => "M",
    9 => "G",
    12 => "T"
)

function positive_finite(name::AbstractString, number::Real)
    isfinite(number) && number > zero(number) || throw(DomainError(
        number,
        "$name must be positive and finite"
    ))
    return number
end

function strand_count(layers::Integer)
    layers >= 0 || throw(DomainError(
        layers,
        "the number of stranded layers must be non-negative"
    ))
    return 1 + 3layers * (layers + 1)
end

"""
Build the physical cable through the current public LineCableModels grammar.

The function is deliberately generic in its physical dimensions. Passing
ordinary floating-point values builds the visual design; passing `Measurement`
values builds the numerically evaluated uncertain design through the same code
path. `strand_layers` remains a discrete topology choice, with `6k` wires in
ring `k` around the centre wire.
"""
function build_design(
        wire_diameter_m::Real,
        insulation_thickness_m::Real,
        strand_layers::Integer
)
    positive_finite("wire diameter", wire_diameter_m)
    positive_finite("insulation thickness", insulation_thickness_m)
    strand_count(strand_layers)

    wire_radius = wire_diameter_m / 2
    wire = LCM.Region(:core_wire, LCM.Disk(wire_radius), COPPER)
    parts = LCM.AbstractCablePart[]
    for layer in 0:strand_layers
        count = iszero(layer) ? 1 : 6layer
        lay_radius = layer * wire_diameter_m
        push!(
            parts,
            LCM.Group(
                :core,
                wire;
                pattern = LCM.Ring(count; r = lay_radius)
            )
        )
    end
    push!(
        parts,
        LCM.Region(
            :core_insulation,
            LCM.Shell(insulation_thickness_m),
            XLPE
        )
    )
    return LCM.build(
        LCM.CableDesign,
        "ichqp2026-cable-constants",
        LCM.Stack(parts)
    )
end

function package_preview_shapes(design::LCM.CableDesign)
    shapes = Any[]
    for region in design.geometry.regions
        context = (
            label = String(region.source.tag),
            group = Symbol("preview_", region.source.tag),
            include_label = false
        )
        append!(
            shapes,
            LCM.DataModel.preview_shapes(region, context)
        )
    end
    return shapes
end

function millimetre_polygon(geometry)
    point(point) = Point2f(1.0e3 * point[1], 1.0e3 * point[2])
    exterior = point.(geometry.exterior)
    interiors = [point.(interior) for interior in geometry.interiors]
    return Makie.GeometryBasics.Polygon(exterior, interiors)
end

"""Project only package-resolved preview geometry into persistent Makie data."""
function preview_projection(design::LCM.CableDesign)
    shapes = package_preview_shapes(design)
    length(shapes) == length(design.geometry.regions) || throw(DimensionMismatch(
        "the package preview projection must return one shape per resolved region"
    ))

    strand_centers = Point2f[]
    strand_polygons = Makie.GeometryBasics.Polygon[]
    strand_diameter_mm = nothing
    strand_color = nothing
    layers = Makie.GeometryBasics.Polygon[]
    layer_colors = Any[]

    for (region, shape) in zip(design.geometry.regions, shapes)
        if region.terminal === :core && region.source.tag === :core_wire
            center = LCM.centroid(region.primitive)
            push!(strand_centers, Point2f(1.0e3 * center[1], 1.0e3 * center[2]))
            push!(strand_polygons, millimetre_polygon(shape.geometry))
            diameter = 2.0e3 * LCM.r_ex(region.primitive)
            if strand_diameter_mm === nothing
                strand_diameter_mm = Float32(diameter)
                strand_color = shape.color
            elseif !isapprox(strand_diameter_mm, diameter)
                throw(ArgumentError(
                    "the cable-constants preview requires one fixed strand diameter"
                ))
            end
        else
            push!(layers, millimetre_polygon(shape.geometry))
            push!(layer_colors, shape.color)
        end
    end

    isempty(strand_centers) && throw(ArgumentError(
        "the package preview contains no resolved conductor strands"
    ))
    isempty(layers) && throw(ArgumentError(
        "the package preview contains no dielectric geometry"
    ))
    diameter_mm = something(strand_diameter_mm)
    outer_diameter_mm = Float64(2.0e3 * LCM.outer_radius(design))
    return (;
        strand_centers,
        strand_polygons,
        strand_diameter_mm = diameter_mm,
        strand_color,
        layers,
        layer_colors,
        outer_diameter_mm
    )
end

function draw_preview!(axis, projection; outlines::Bool = true)
    layer_plot = poly!(
        axis,
        projection.layers;
        color = projection.layer_colors,
        strokecolor = :black,
        strokewidth = outlines ? 1.1 : 0
    )
    strand_plot = poly!(
        axis,
        projection.strand_polygons;
        color = fill(projection.strand_color, length(projection.strand_polygons)),
        strokecolor = :black,
        strokewidth = outlines ? 0.8 : 0
    )
    return (; layer_plot, strand_plot)
end

function update_preview_plots!(layer_plot, strand_plot, projection)
    Makie.update!(
        layer_plot,
        projection.layers;
        color = projection.layer_colors
    )
    Makie.update!(
        strand_plot,
        projection.strand_polygons;
        color = fill(projection.strand_color, length(projection.strand_polygons))
    )
    return projection
end

function translate_polygon(polygon, offset::Point2f)
    exterior = [point + offset for point in polygon.exterior]
    interiors = [
        [point + offset for point in interior] for interior in polygon.interiors
    ]
    return Makie.GeometryBasics.Polygon(exterior, interiors)
end

"""Translate preview geometry for display only; no second domain object is built."""
function translated_preview_projection(projection, offset::Point2f)
    return merge(
        projection,
        (;
            strand_centers = [point + offset for point in projection.strand_centers],
            strand_polygons = [
                translate_polygon(polygon, offset)
                for polygon in projection.strand_polygons
            ],
            layers = [
                translate_polygon(polygon, offset)
                for polygon in projection.layers
            ]
        )
    )
end

function preview_limits!(axis, radius_limit_mm::Real)
    limits!(
        axis,
        -radius_limit_mm,
        CLONE_LATERAL_OFFSET_MM + radius_limit_mm,
        -radius_limit_mm,
        radius_limit_mm
    )
    return axis
end

function millimetre_points(points)
    converted = Point2f[
        Point2f(1.0e3 * point[1], 1.0e3 * point[2]) for point in points
    ]
    if length(converted) > 1 && isapprox(first(converted), last(converted))
        pop!(converted)
    end
    return converted
end

function signed_area(points)
    isempty(points) && return 0.0
    return 0.5 * sum(eachindex(points)) do index
        next_index = index == lastindex(points) ? firstindex(points) : index + 1
        return points[index][1] * points[next_index][2] -
               points[next_index][1] * points[index][2]
    end
end

function boundary_signature(points)
    xs = first.(points)
    ys = last.(points)
    return (
        round(minimum(xs); digits = 4),
        round(maximum(xs); digits = 4),
        round(minimum(ys); digits = 4),
        round(maximum(ys); digits = 4),
        round(abs(signed_area(points)); digits = 4)
    )
end

"""
Project the package's resolved preview shapes into stable physical boundaries.

The IDs follow resolved-region identity rather than drawing order. Coincident
inner/outer contours are retained once, so a shared material interface cannot
acquire a spurious darker uncertainty band.
"""
function boundary_projection(
        design::LCM.CableDesign;
        unique_physical::Bool = true
)
    shapes = package_preview_shapes(design)
    counts = Dict{Tuple{Any, Symbol}, Int}()
    boundaries = NamedTuple[]
    signatures = Set{Any}()

    for (region, shape) in zip(design.geometry.regions, shapes)
        key = (region.terminal, region.source.tag)
        ordinal = get(counts, key, 0) + 1
        counts[key] = ordinal
        contours = ((:outer, 0, shape.geometry.exterior),)
        contours = (
            contours...,
            ((:inner, index, points) for
             (index, points) in enumerate(shape.geometry.interiors))...
        )
        for (side, contour_index, contour) in contours
            points = millimetre_points(contour)
            signature = boundary_signature(points)
            unique_physical && signature in signatures && continue
            push!(signatures, signature)
            push!(boundaries, (;
                id = (key..., ordinal, side, contour_index),
                points,
                signature
            ))
        end
    end
    return boundaries
end

function resample_boundary(points, count::Integer = BOUNDARY_SAMPLES)
    length(points) >= 3 || throw(ArgumentError(
        "a physical boundary requires at least three points"
    ))
    count >= 3 || throw(ArgumentError(
        "a resampled physical boundary requires at least three points"
    ))
    segment_lengths = Float64[]
    for index in eachindex(points)
        next_index = index == lastindex(points) ? firstindex(points) : index + 1
        push!(segment_lengths, norm(points[next_index] - points[index]))
    end
    perimeter = sum(segment_lengths)
    perimeter > 0 || throw(ArgumentError("a physical boundary must have area"))
    cumulative = cumsum(segment_lengths)
    result = Point2f[]
    sizehint!(result, count)
    for distance in range(0.0, perimeter; length = count + 1)[1:end-1]
        segment = searchsortedfirst(cumulative, distance)
        previous = segment == 1 ? 0.0 : cumulative[segment - 1]
        fraction = (distance - previous) / segment_lengths[segment]
        next_segment = segment == lastindex(points) ? firstindex(points) : segment + 1
        push!(result, Point2f(
            (1 - fraction) * points[segment] + fraction * points[next_segment]
        ))
    end
    return result
end

function boundary_strip_mesh(
        nominal_boundaries,
        lower_design::LCM.CableDesign,
        upper_design::LCM.CableDesign
)
    lower = Dict(
        boundary.id => boundary.points for
        boundary in boundary_projection(lower_design; unique_physical = false)
    )
    upper = Dict(
        boundary.id => boundary.points for
        boundary in boundary_projection(upper_design; unique_physical = false)
    )
    vertices = Point2f[]
    faces = Makie.GeometryBasics.GLTriangleFace[]
    affected = Any[]

    for boundary in nominal_boundaries
        haskey(lower, boundary.id) && haskey(upper, boundary.id) || throw(KeyError(
            "perturbed designs do not contain nominal boundary $(boundary.id)"
        ))
        lower_points = resample_boundary(lower[boundary.id])
        upper_points = resample_boundary(upper[boundary.id])
        maximum(norm.(upper_points .- lower_points)) <= 1.0e-5 && continue
        push!(affected, boundary.id)
        offset = length(vertices)
        for index in eachindex(lower_points)
            push!(vertices, lower_points[index], upper_points[index])
        end
        count = length(lower_points)
        for index in 1:count
            next_index = index == count ? 1 : index + 1
            lower_index = offset + 2index - 1
            upper_index = offset + 2index
            next_lower = offset + 2next_index - 1
            next_upper = offset + 2next_index
            push!(
                faces,
                Makie.GeometryBasics.GLTriangleFace(
                    lower_index,
                    upper_index,
                    next_upper
                ),
                Makie.GeometryBasics.GLTriangleFace(
                    lower_index,
                    next_upper,
                    next_lower
                )
            )
        end
    end
    return (;
        mesh = Makie.GeometryBasics.Mesh(vertices, faces),
        affected,
        vertex_count = length(vertices),
        face_count = length(faces)
    )
end

function perturbed_designs(input, source::Symbol)
    wire_sigma = ENVELOPE_COVERAGE * input.wire_diameter *
                 input.wire_uncertainty / 100
    insulation_sigma = ENVELOPE_COVERAGE * input.insulation_thickness *
                       input.insulation_uncertainty / 100
    if source === :wire_diameter
        lower = build_design(
            input.wire_diameter - wire_sigma,
            input.insulation_thickness,
            input.strand_layers
        )
        upper = build_design(
            input.wire_diameter + wire_sigma,
            input.insulation_thickness,
            input.strand_layers
        )
    elseif source === :insulation_thickness
        lower = build_design(
            input.wire_diameter,
            input.insulation_thickness - insulation_sigma,
            input.strand_layers
        )
        upper = build_design(
            input.wire_diameter,
            input.insulation_thickness + insulation_sigma,
            input.strand_layers
        )
    else
        throw(ArgumentError("unknown uncertainty source: $source"))
    end
    return (; lower, upper)
end

function envelope_projection(
        nominal_design::LCM.CableDesign,
        input,
        source::Symbol
)
    boundaries = boundary_projection(nominal_design)
    perturbed = perturbed_designs(input, source)
    return boundary_strip_mesh(
        boundaries,
        perturbed.lower,
        perturbed.upper
    )
end

function outline_segments(boundaries)
    points = Point2f[]
    for boundary in boundaries
        contour = resample_boundary(boundary.points)
        for index in eachindex(contour)
            next_index = index == lastindex(contour) ? firstindex(contour) : index + 1
            push!(points, contour[index], contour[next_index])
        end
    end
    return points
end

function preview_figure(design::LCM.CableDesign)
    projection = preview_projection(design)
    clone_projection = translated_preview_projection(
        projection,
        Point2f(CLONE_LATERAL_OFFSET_MM, 0)
    )
    radius_limit_mm = Ref(Float32(
        PREVIEW_MARGIN * projection.outer_diameter_mm / 2
    ))
    title = Observable(
        "Two $(length(projection.strand_centers))-strand cables · 1.0 m spacing · $(round(projection.outer_diameter_mm; digits = 1)) mm OD"
    )
    figure = Figure(size = (860, 700), backgroundcolor = PLOT_BACKGROUND)
    axis = Axis(
        figure[1, 1];
        title,
        titlecolor = PLOT_TEXT,
        backgroundcolor = PLOT_BACKGROUND,
        aspect = DataAspect()
    )
    hidedecorations!(axis)
    hidespines!(axis)
    preview = draw_preview!(axis, projection)
    clone = draw_preview!(axis, clone_projection)
    preview_limits!(axis, radius_limit_mm[])
    return (;
        figure,
        axis,
        layer_plot = preview.layer_plot,
        strand_plot = preview.strand_plot,
        clone_layer_plot = clone.layer_plot,
        clone_strand_plot = clone.strand_plot,
        title,
        radius_limit_mm,
        projection
    )
end

function update_preview!(handle, design::LCM.CableDesign)
    projection = preview_projection(design)
    clone_projection = translated_preview_projection(
        projection,
        Point2f(CLONE_LATERAL_OFFSET_MM, 0)
    )
    update_preview_plots!(
        handle.layer_plot,
        handle.strand_plot,
        projection
    )
    update_preview_plots!(
        handle.clone_layer_plot,
        handle.clone_strand_plot,
        clone_projection
    )
    handle.title[] = "Two $(length(projection.strand_centers))-strand cables · 1.0 m spacing · $(round(projection.outer_diameter_mm; digits = 1)) mm OD"
    handle.radius_limit_mm[] = Float32(
        PREVIEW_MARGIN * projection.outer_diameter_mm / 2
    )
    preview_limits!(handle.axis, handle.radius_limit_mm[])
    return projection
end

function envelope_axis_limit(strand_layers::Integer)
    maximum_wire = 1.0e-3 * last(WIRE_DIAMETER_RANGE_MM) *
                   (1 + last(UNCERTAINTY_RANGE_PCT) / 100)
    maximum_insulation = 1.0e-3 * last(INSULATION_THICKNESS_RANGE_MM) *
                         (1 + last(UNCERTAINTY_RANGE_PCT) / 100)
    design = build_design(maximum_wire, maximum_insulation, strand_layers)
    return Float32(PREVIEW_MARGIN * 1.0e3 * LCM.outer_radius(design))
end

function uncertainty_figure(design::LCM.CableDesign, input)
    projection = preview_projection(design)
    boundaries = boundary_projection(design)
    wire = envelope_projection(design, input, :wire_diameter)
    insulation = envelope_projection(design, input, :insulation_thickness)
    radius_limit_mm = Ref(envelope_axis_limit(input.strand_layers))
    title = Observable(
        "Package geometry · ±$(Int(ENVELOPE_COVERAGE))σ envelopes"
    )
    figure = Figure(size = (900, 700), backgroundcolor = PLOT_BACKGROUND)
    axis = Axis(
        figure[1, 1];
        title,
        titlecolor = PLOT_TEXT,
        backgroundcolor = PLOT_BACKGROUND,
        aspect = DataAspect()
    )
    hidedecorations!(axis)
    hidespines!(axis)

    preview = draw_preview!(axis, projection; outlines = false)
    wire_envelope_plot = mesh!(
        axis,
        wire.mesh;
        color = WIRE_ENVELOPE_COLOR,
        transparency = true,
        shading = NoShading
    )
    insulation_envelope_plot = mesh!(
        axis,
        insulation.mesh;
        color = INSULATION_ENVELOPE_COLOR,
        transparency = true,
        shading = NoShading
    )
    outline_plot = linesegments!(
        axis,
        outline_segments(boundaries);
        color = RGBAf(0.05, 0.06, 0.07, 0.92),
        linewidth = 1.2
    )
    limits!(
        axis,
        -radius_limit_mm[],
        radius_limit_mm[],
        -radius_limit_mm[],
        radius_limit_mm[]
    )
    return (;
        figure,
        axis,
        layer_plot = preview.layer_plot,
        strand_plot = preview.strand_plot,
        wire_envelope_plot,
        insulation_envelope_plot,
        outline_plot,
        title,
        radius_limit_mm,
        nominal_design = Ref{Any}(design),
        boundaries = Ref(boundaries),
        projection,
        wire,
        insulation
    )
end

function update_uncertainty_nominal!(handle, design::LCM.CableDesign, strand_layers)
    projection = preview_projection(design)
    boundaries = boundary_projection(design)
    update_preview_plots!(
        handle.layer_plot,
        handle.strand_plot,
        projection
    )
    Makie.update!(
        handle.outline_plot;
        arg1 = outline_segments(boundaries)
    )
    handle.nominal_design[] = design
    handle.boundaries[] = boundaries
    handle.radius_limit_mm[] = envelope_axis_limit(strand_layers)
    limits!(
        handle.axis,
        -handle.radius_limit_mm[],
        handle.radius_limit_mm[],
        -handle.radius_limit_mm[],
        handle.radius_limit_mm[]
    )
    return projection
end

function update_uncertainty_source!(handle, input, source::Symbol)
    envelope = envelope_projection(handle.nominal_design[], input, source)
    plot = source === :wire_diameter ? handle.wire_envelope_plot :
           source === :insulation_thickness ? handle.insulation_envelope_plot :
           throw(ArgumentError("unknown uncertainty source: $source"))
    Makie.update!(plot; arg1 = envelope.mesh)
    return envelope
end

function engineering(number::Real, unit::AbstractString)
    numeric = Float64(number)
    isfinite(numeric) || return "$numeric $unit"
    iszero(numeric) && return "0 $unit"
    exponent = clamp(3floor(Int, log10(abs(numeric)) / 3), -12, 12)
    scaled = numeric / 10.0^exponent
    return @sprintf("%.4g %s%s", scaled, SI_PREFIX[exponent], unit)
end

function display_measurement(number::Measurement, quantity::Symbol)
    factor = getproperty(RESULT_FACTORS, quantity)
    unit = getproperty(RESULT_UNIT_LABELS, quantity)
    displayed = number * factor
    return @sprintf("%.4g ± %.2g %s",
        value(displayed),
        uncertainty(displayed),
        unit)
end

function display_physical_measurement(
        number::Measurement,
        factor::Real,
        unit::AbstractString
)
    displayed = number * factor
    return @sprintf("%.4g ± %.2g %s",
        value(displayed),
        uncertainty(displayed),
        unit)
end

function measured_inputs(input)
    wire_nominal = input.wire_diameter
    insulation_nominal = input.insulation_thickness
    wire = measurement(
        wire_nominal,
        wire_nominal * input.wire_uncertainty / 100
    )
    insulation = measurement(
        insulation_nominal,
        insulation_nominal * input.insulation_uncertainty / 100
    )
    return wire, insulation
end

function conductor_diagnostics(design::LCM.CableDesign)
    strands = filter(design.geometry.regions) do region
        region.terminal === :core && region.source.tag === :core_wire
    end
    isempty(strands) && throw(ArgumentError(
        "the cable design contains no resolved core strands"
    ))
    conductor_area = sum(
        region -> LCM.area(region.primitive),
        strands;
        init = zero(eltype(design))
    )
    dc_resistance = COPPER.rho / conductor_area
    return (; conductor_area, dc_resistance)
end

function numerical_design_key(input)
    return (
        input.strand_layers,
        input.wire_diameter,
        input.wire_uncertainty,
        input.insulation_thickness,
        input.insulation_uncertainty
    )
end

function build_uncertain_design(input)
    wire, insulation = measured_inputs(input)
    return build_design(wire, insulation, input.strand_layers)
end

function compute_constants(input, design::LCM.CableDesign)
    started = time_ns()
    constants = LCM.CableConstants(
        design;
        frequency = input.frequency,
        formulation = CABLE_CONSTANTS_FORMULATION
    )
    observations = (
        R = only(LCM.observe(constants, LCM.R)),
        L = only(LCM.observe(constants, LCM.L)),
        C = only(LCM.observe(constants, LCM.C)),
        G = only(LCM.observe(constants, LCM.G))
    )
    diagnostics = conductor_diagnostics(design)
    elapsed_ms = (time_ns() - started) / 1.0e6
    return (;
        input,
        design,
        constants,
        observations,
        diagnostics,
        elapsed_ms,
        finished_at = time()
    )
end

compute_constants(input) = compute_constants(input, build_uncertain_design(input))

mutable struct NumericalWorker
    queue::Channel{Any}
    generation::Int
    closed::Bool
    lock::ReentrantLock
    task::Union{Nothing, Task}
    design_key::Any
    design::Any
end

function NumericalWorker(design_key, design)
    NumericalWorker(
        Channel{Any}(1),
        0,
        false,
        ReentrantLock(),
        nothing,
        design_key,
        design
    )
end

function enqueue!(worker::NumericalWorker, input)
    lock(worker.lock) do
        worker.closed && return worker.generation
        worker.generation += 1
        while isready(worker.queue)
            take!(worker.queue)
        end
        put!(worker.queue, (generation = worker.generation, input))
        return worker.generation
    end
end

function close!(worker::NumericalWorker)
    lock(worker.lock) do
        worker.closed && return nothing
        worker.closed = true
        close(worker.queue)
    end
    return nothing
end

function worker_is_current(worker::NumericalWorker, generation::Int)
    return lock(worker.lock) do
        !worker.closed && worker.generation == generation
    end
end

function run_worker!(worker::NumericalWorker, result_state, status_state)
    while !worker.closed
        item = try
            take!(worker.queue)
        catch error
            error isa InvalidStateException && break
            rethrow()
        end

        sleep(0.1)
        while isready(worker.queue)
            item = take!(worker.queue)
        end
        worker_is_current(worker, item.generation) || continue
        status_state[] = (
            phase = :updating,
            label = "Updating package cable constants…"
        )
        try
            design_key = numerical_design_key(item.input)
            if design_key != worker.design_key
                worker.design = build_uncertain_design(item.input)
                worker.design_key = design_key
            end
            calculation = Threads.@spawn compute_constants(
                item.input,
                worker.design
            )
            result = fetch(calculation)
            worker_is_current(worker, item.generation) || continue
            result_state[] = result
            status_state[] = (
                phase = :ready,
                label = @sprintf("Package result ready · %.1f ms",
                    result.elapsed_ms)
            )
        catch error
            worker_is_current(worker, item.generation) || continue
            status_state[] = (
                phase = :error,
                label = "Cable-constants update failed: $(sprint(showerror, error))"
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

function link_slider_values!(session::Session, left, right)
    left_to_right = on(session, left.value) do value
        isequal(right.value[], value) || (right.value[] = value)
        return nothing
    end
    right_to_left = on(session, right.value) do value
        isequal(left.value[], value) || (left.value[] = value)
        return nothing
    end
    return (; left_to_right, right_to_left)
end

function setup(session::Session)
    frequency_slider = make_slider(
        FREQUENCY_LOG_RANGE,
        log10(DEFAULT_INPUT.frequency),
        "cable-constants-frequency-input",
        "Cable-constants frequency on a logarithmic scale"
    )
    strand_layers_slider = make_slider(
        STRAND_LAYER_RANGE,
        DEFAULT_INPUT.strand_layers,
        "cable-strand-layers-input",
        "Number of concentric stranded conductor layers"
    )
    wire_diameter_slider = make_slider(
        WIRE_DIAMETER_RANGE_MM,
        1.0e3 * DEFAULT_INPUT.wire_diameter,
        "cable-wire-diameter-input",
        "Individual conductor wire diameter in millimetres"
    )
    wire_uncertainty_slider = make_slider(
        UNCERTAINTY_RANGE_PCT,
        DEFAULT_INPUT.wire_uncertainty,
        "cable-wire-uncertainty-input",
        "Wire diameter standard uncertainty in percent"
    )
    insulation_thickness_slider = make_slider(
        INSULATION_THICKNESS_RANGE_MM,
        1.0e3 * DEFAULT_INPUT.insulation_thickness,
        "cable-insulation-thickness-input",
        "Insulation thickness in millimetres"
    )
    insulation_uncertainty_slider = make_slider(
        UNCERTAINTY_RANGE_PCT,
        DEFAULT_INPUT.insulation_uncertainty,
        "cable-insulation-uncertainty-input",
        "Insulation thickness standard uncertainty in percent"
    )
    geometry_strand_layers_slider = make_slider(
        STRAND_LAYER_RANGE,
        DEFAULT_INPUT.strand_layers,
        "geometry-envelope-strand-layers-input",
        "Number of concentric stranded conductor layers"
    )
    geometry_wire_diameter_slider = make_slider(
        WIRE_DIAMETER_RANGE_MM,
        1.0e3 * DEFAULT_INPUT.wire_diameter,
        "geometry-envelope-wire-diameter-input",
        "Individual conductor wire diameter in millimetres"
    )
    geometry_wire_uncertainty_slider = make_slider(
        UNCERTAINTY_RANGE_PCT,
        DEFAULT_INPUT.wire_uncertainty,
        "geometry-envelope-wire-uncertainty-input",
        "Wire diameter standard uncertainty in percent"
    )
    geometry_insulation_thickness_slider = make_slider(
        INSULATION_THICKNESS_RANGE_MM,
        1.0e3 * DEFAULT_INPUT.insulation_thickness,
        "geometry-envelope-insulation-thickness-input",
        "Insulation thickness in millimetres"
    )
    geometry_insulation_uncertainty_slider = make_slider(
        UNCERTAINTY_RANGE_PCT,
        DEFAULT_INPUT.insulation_uncertainty,
        "geometry-envelope-insulation-uncertainty-input",
        "Insulation thickness standard uncertainty in percent"
    )
    linked_sliders = (
        link_slider_values!(
            session,
            strand_layers_slider,
            geometry_strand_layers_slider
        ),
        link_slider_values!(
            session,
            wire_diameter_slider,
            geometry_wire_diameter_slider
        ),
        link_slider_values!(
            session,
            wire_uncertainty_slider,
            geometry_wire_uncertainty_slider
        ),
        link_slider_values!(
            session,
            insulation_thickness_slider,
            geometry_insulation_thickness_slider
        ),
        link_slider_values!(
            session,
            insulation_uncertainty_slider,
            geometry_insulation_uncertainty_slider
        )
    )

    input = map(
        (
            log_frequency,
            strand_layers,
            wire_diameter,
            wire_uncertainty,
            insulation_thickness,
            insulation_uncertainty
        ) -> (;
            frequency = 10.0^Float64(log_frequency),
            strand_layers = Int(strand_layers),
            wire_diameter = 1.0e-3 * Float64(wire_diameter),
            wire_uncertainty = Float64(wire_uncertainty),
            insulation_thickness = 1.0e-3 * Float64(insulation_thickness),
            insulation_uncertainty = Float64(insulation_uncertainty)
        ),
        session,
        frequency_slider.value,
        strand_layers_slider.value,
        wire_diameter_slider.value,
        wire_uncertainty_slider.value,
        insulation_thickness_slider.value,
        insulation_uncertainty_slider.value
    )
    geometry_input = map(
        (strand_layers, wire_diameter, insulation_thickness) -> (;
            strand_layers = Int(strand_layers),
            wire_diameter = 1.0e-3 * Float64(wire_diameter),
            insulation_thickness = 1.0e-3 * Float64(insulation_thickness)
        ),
        session,
        strand_layers_slider.value,
        wire_diameter_slider.value,
        insulation_thickness_slider.value
    )

    initial_design = build_design(
        DEFAULT_INPUT.wire_diameter,
        DEFAULT_INPUT.insulation_thickness,
        DEFAULT_INPUT.strand_layers
    )
    plot = preview_figure(initial_design)
    uncertainty_plot = uncertainty_figure(initial_design, DEFAULT_INPUT)
    visual_state = Observable(plot.projection)
    geometry_observer = on(geometry_input) do geometry
        design = build_design(
            geometry.wire_diameter,
            geometry.insulation_thickness,
            geometry.strand_layers
        )
        visual_state[] = update_preview!(plot, design)
        update_uncertainty_nominal!(
            uncertainty_plot,
            design,
            geometry.strand_layers
        )
        latest = merge(input[], geometry)
        update_uncertainty_source!(
            uncertainty_plot,
            latest,
            :wire_diameter
        )
        update_uncertainty_source!(
            uncertainty_plot,
            latest,
            :insulation_thickness
        )
        return nothing
    end
    wire_envelope_observer = on(wire_uncertainty_slider.value) do uncertainty
        latest = merge(input[], (; wire_uncertainty = Float64(uncertainty)))
        update_uncertainty_source!(
            uncertainty_plot,
            latest,
            :wire_diameter
        )
        return nothing
    end
    insulation_envelope_observer = on(
        insulation_uncertainty_slider.value
    ) do uncertainty
        latest = merge(
            input[],
            (; insulation_uncertainty = Float64(uncertainty))
        )
        update_uncertainty_source!(
            uncertainty_plot,
            latest,
            :insulation_thickness
        )
        return nothing
    end

    initial_numerical_design = build_uncertain_design(input[])
    result_state = Observable(compute_constants(input[], initial_numerical_design))
    status_state = Observable((
        phase = :ready,
        label = @sprintf("Package result ready · %.1f ms",
            result_state[].elapsed_ms)
    ))
    worker = NumericalWorker(
        numerical_design_key(input[]),
        initial_numerical_design
    )
    worker.task = @async run_worker!(worker, result_state, status_state)
    errormonitor(worker.task)
    numerical_observer = on(input) do latest
        status_state[] = (
            phase = :queued,
            label = "Updating… waiting for controls to settle"
        )
        enqueue!(worker, latest)
        return nothing
    end

    reset_button = Bonito.Button(
        "Reset plot view";
        style = nothing,
        id = "cable-constants-reset-view-control",
        class = "lc-action-button",
        type = "button",
        ariaLabel = "Reset the cable preview plot view"
    )
    selector_reset_button = Bonito.Button(
        "Reset selectors";
        style = nothing,
        id = "cable-constants-reset-selectors-control",
        class = "lc-action-button",
        type = "button",
        ariaLabel = "Reset all cable-constants selectors to their initial values"
    )
    reset_observer = on(reset_button.value) do clicked
        clicked || return nothing
        preview_limits!(plot.axis, plot.radius_limit_mm[])
        return nothing
    end
    selector_reset_observer = on(selector_reset_button.value) do clicked
        clicked || return nothing
        reset_selectors!(
            frequency_slider => log10(DEFAULT_INPUT.frequency),
            strand_layers_slider => DEFAULT_INPUT.strand_layers,
            wire_diameter_slider => 1.0e3 * DEFAULT_INPUT.wire_diameter,
            wire_uncertainty_slider => DEFAULT_INPUT.wire_uncertainty,
            insulation_thickness_slider =>
                1.0e3 * DEFAULT_INPUT.insulation_thickness,
            insulation_uncertainty_slider =>
                DEFAULT_INPUT.insulation_uncertainty
        )
        return nothing
    end
    geometry_reset_button = Bonito.Button(
        "Reset plot view";
        style = nothing,
        id = "geometry-envelope-reset-view-control",
        class = "lc-action-button",
        type = "button",
        ariaLabel = "Reset the geometry uncertainty plot view"
    )
    geometry_selector_reset_button = Bonito.Button(
        "Reset selectors";
        style = nothing,
        id = "geometry-envelope-reset-selectors-control",
        class = "lc-action-button",
        type = "button",
        ariaLabel = "Reset all geometry-envelope selectors to their initial values"
    )
    geometry_reset_observer = on(geometry_reset_button.value) do clicked
        clicked || return nothing
        limits!(
            uncertainty_plot.axis,
            -uncertainty_plot.radius_limit_mm[],
            uncertainty_plot.radius_limit_mm[],
            -uncertainty_plot.radius_limit_mm[],
            uncertainty_plot.radius_limit_mm[]
        )
        return nothing
    end
    geometry_selector_reset_observer = on(
        geometry_selector_reset_button.value
    ) do clicked
        clicked || return nothing
        reset_selectors!(
            geometry_strand_layers_slider => DEFAULT_INPUT.strand_layers,
            geometry_wire_diameter_slider =>
                1.0e3 * DEFAULT_INPUT.wire_diameter,
            geometry_wire_uncertainty_slider => DEFAULT_INPUT.wire_uncertainty,
            geometry_insulation_thickness_slider =>
                1.0e3 * DEFAULT_INPUT.insulation_thickness,
            geometry_insulation_uncertainty_slider =>
                DEFAULT_INPUT.insulation_uncertainty
        )
        return nothing
    end

    on(session.on_close) do _
        off(geometry_observer)
        off(wire_envelope_observer)
        off(insulation_envelope_observer)
        off(numerical_observer)
        off(reset_observer)
        off(selector_reset_observer)
        off(geometry_reset_observer)
        off(geometry_selector_reset_observer)
        close!(worker)
        return nothing
    end
    return (;
        input,
        geometry_input,
        visual_state,
        result_state,
        status_state,
        plot,
        uncertainty_plot,
        worker,
        frequency_slider,
        strand_layers_slider,
        wire_diameter_slider,
        wire_uncertainty_slider,
        insulation_thickness_slider,
        insulation_uncertainty_slider,
        geometry_strand_layers_slider,
        geometry_wire_diameter_slider,
        geometry_wire_uncertainty_slider,
        geometry_insulation_thickness_slider,
        geometry_insulation_uncertainty_slider,
        linked_sliders,
        reset_button,
        selector_reset_button,
        geometry_reset_button,
        geometry_selector_reset_button
    )
end

function build_geometry_page(session::Session, state)
    strand_layers = map(
        input -> string(input.strand_layers),
        session,
        state.input
    )
    wire_diameter = map(
        input -> @sprintf("%.1f mm", 1.0e3 * input.wire_diameter),
        session,
        state.input
    )
    wire_uncertainty = map(
        input -> @sprintf(
            "%.1f %% · σd = %.3f mm",
            input.wire_uncertainty,
            1.0e3 * input.wire_diameter * input.wire_uncertainty / 100
        ),
        session,
        state.input
    )
    insulation_thickness = map(
        input -> @sprintf("%.1f mm", 1.0e3 * input.insulation_thickness),
        session,
        state.input
    )
    insulation_uncertainty = map(
        input -> @sprintf(
            "%.1f %% · σt = %.3f mm",
            input.insulation_uncertainty,
            1.0e3 * input.insulation_thickness *
            input.insulation_uncertainty / 100
        ),
        session,
        state.input
    )
    outer_diameter = map(
        projection -> @sprintf("%.1f mm", projection.outer_diameter_mm),
        session,
        state.visual_state
    )
    resolved_strands = map(
        projection -> string(length(projection.strand_centers)),
        session,
        state.visual_state
    )

    plot = webpart(
        wgl_figure(state.uncertainty_plot.figure);
        kind = :plot,
        title = "Geometry influence in physical space",
        meta = "nominal fill · source envelopes · nominal outlines",
        id = "cable-geometry-uncertainty-preview",
        body_class = "lc-transition-plot-host"
    )
    controls = webpart(
        control(
            "Stranded layers",
            state.geometry_strand_layers_slider;
            value = strand_layers,
            id = "geometry-envelope-strand-layers-control"
        ),
        control(
            "Wire diameter",
            state.geometry_wire_diameter_slider;
            value = wire_diameter,
            id = "geometry-envelope-wire-diameter-control"
        ),
        control(
            "Wire-diameter standard uncertainty",
            state.geometry_wire_uncertainty_slider;
            value = wire_uncertainty,
            id = "geometry-envelope-wire-uncertainty-control"
        ),
        control(
            "Insulation thickness",
            state.geometry_insulation_thickness_slider;
            value = insulation_thickness,
            id = "geometry-envelope-insulation-thickness-control"
        ),
        control(
            "Insulation-thickness standard uncertainty",
            state.geometry_insulation_uncertainty_slider;
            value = insulation_uncertainty,
            id = "geometry-envelope-insulation-uncertainty-control"
        ),
        color_key(
            "Wire diameter ±1σ" => "is-wire-envelope",
            "Insulation thickness ±1σ" => "is-insulation-envelope";
            separator = "+",
            class = "lc-geometry-envelope-key"
        ),
        value_list(
            "Resolved strands (1 + Σ6k)" => resolved_strands,
            "Nominal cable outer diameter" => outer_diameter,
            "Coverage factor" => "k = $(Int(ENVELOPE_COVERAGE))"
        ),
        DOM.p(
            "Colored regions show ±1 standard-uncertainty geometry envelopes. Overlapping colors indicate overlapping spatial influence, not probability density.";
            class = "lc-geometry-envelope-caption"
        ),
        action_row(
            state.geometry_reset_button,
            state.geometry_selector_reset_button
        ),
        diagnostic(
            "lower/nominal/upper CableDesign → package preview geometry";
            suffix = " · two batched persistent overlay meshes"
        );
        kind = :controls,
        class = "lc-cable-constants-controls lc-geometry-envelope-controls"
    )
    canvas = layout_two_columns(
        plot,
        controls;
        ratio = :wide_narrow,
        height = "100%",
        class = "lc-transition-grid lc-cable-constants-grid lc-geometry-envelope-grid"
    )
    return (body = slide(
        "Geometry uncertainty envelopes",
        canvas;
        lede = md"""
        The nominal cable comes from the package preview. Cyan and magenta overlays show where the same physical boundaries move when the two uncertain dimensions are rebuilt at ±1σ; strand expansion and translation are both retained.
        """
    ),)
end

function build_page(session::Session, state)
    frequency = map(
        input -> engineering(input.frequency, "Hz"),
        session,
        state.input
    )
    strand_layers = map(
        input -> string(input.strand_layers),
        session,
        state.input
    )
    wire_diameter = map(
        input -> @sprintf("%.1f mm", 1.0e3 * input.wire_diameter),
        session,
        state.input
    )
    wire_uncertainty = map(
        input -> @sprintf("%.1f %%", input.wire_uncertainty),
        session,
        state.input
    )
    insulation_thickness = map(
        input -> @sprintf("%.1f mm", 1.0e3 * input.insulation_thickness),
        session,
        state.input
    )
    insulation_uncertainty = map(
        input -> @sprintf("%.1f %%", input.insulation_uncertainty),
        session,
        state.input
    )
    outer_diameter = map(
        projection -> @sprintf("%.1f mm", projection.outer_diameter_mm),
        session,
        state.visual_state
    )
    resolved_strands = map(
        projection -> string(length(projection.strand_centers)),
        session,
        state.visual_state
    )
    conductor_area = map(
        result -> display_physical_measurement(
            result.diagnostics.conductor_area,
            1.0e6,
            "mm²"
        ),
        session,
        state.result_state
    )
    dc_resistance = map(
        result -> display_physical_measurement(
            result.diagnostics.dc_resistance,
            1.0e3,
            "Ω/km"
        ),
        session,
        state.result_state
    )
    resistance = map(
        result -> display_measurement(result.observations.R, :R),
        session,
        state.result_state
    )
    inductance = map(
        result -> display_measurement(result.observations.L, :L),
        session,
        state.result_state
    )
    capacitance = map(
        result -> display_measurement(result.observations.C, :C),
        session,
        state.result_state
    )
    conductance = map(
        result -> display_measurement(result.observations.G, :G),
        session,
        state.result_state
    )
    result_frequency = map(
        result -> engineering(result.input.frequency, "Hz"),
        session,
        state.result_state
    )
    status = map(snapshot -> snapshot.label, session, state.status_state)

    plot = webpart(
        wgl_figure(state.plot.figure);
        kind = :plot,
        title = "Package preview geometry",
        meta = "persistent Makie scene",
        id = "cable-constants-preview",
        body_class = "lc-transition-plot-host"
    )
    controls = webpart(
        control(
            "Frequency",
            state.frequency_slider;
            value = frequency,
            id = "cable-constants-frequency-control"
        ),
        control(
            "Stranded layers",
            state.strand_layers_slider;
            value = strand_layers,
            id = "cable-strand-layers-control"
        ),
        control(
            "Wire diameter",
            state.wire_diameter_slider;
            value = wire_diameter,
            id = "cable-wire-diameter-control"
        ),
        control(
            "Wire-diameter standard uncertainty [%]",
            state.wire_uncertainty_slider;
            value = wire_uncertainty,
            id = "cable-wire-uncertainty-control"
        ),
        control(
            "Insulation thickness",
            state.insulation_thickness_slider;
            value = insulation_thickness,
            id = "cable-insulation-thickness-control"
        ),
        control(
            "Insulation-thickness standard uncertainty [%]",
            state.insulation_uncertainty_slider;
            value = insulation_uncertainty,
            id = "cable-insulation-uncertainty-control"
        ),
        value_list(
            # "Resolved strands (1 + Σ6k)" => resolved_strands,
            "Cable outer diameter" => outer_diameter,
            "Conductor area" => conductor_area,
            "Rdc = ρ/A" => dc_resistance
        ),
        DOM.div(
            DOM.div(
                DOM.strong("One-frequency package result"),
                DOM.span(result_frequency);
                class = "lc-cable-result-heading"
            ),
            value_list(
                "R" => resistance,
                "L" => inductance,
                "C" => capacitance,
                "G" => conductance;
                class = "lc-cable-constants-values"
            );
            class = "lc-cable-result-panel"
        ),
        action_row(
            state.reset_button,
            state.selector_reset_button
        ),
        status_line(status; class = "lc-cable-constants-status"),
        diagnostic(
            "CableDesign → CableConstants → observe(R, L, C, G)";
            suffix = " · latest-only debounced numerical lane"
        );
        kind = :controls,
        class = "lc-cable-constants-controls"
    )
    canvas = layout_two_columns(
        plot,
        controls;
        ratio = :wide_narrow,
        height = "100%",
        class = "lc-transition-grid lc-cable-constants-grid"
    )
    return (body = slide(
        "Cable constants under uncertainty",
        canvas;
        lede = md"""
        Nominal geometry updates the persistent package preview immediately. After the controls settle, the same builder reconstructs a `CableDesign{Measurement}` and the package evaluates one frequency point with correlated uncertainty propagation.
        """
    ),)
end

function warm_cable_constants!()
    try
        set_preflight!(
            PREFLIGHT_STATE,
            :nominal_design,
            0.14,
            "Building the nominal cable through LineCableModels…"
        )
        yield()
        design = build_design(
            DEFAULT_INPUT.wire_diameter,
            DEFAULT_INPUT.insulation_thickness,
            DEFAULT_INPUT.strand_layers
        )

        set_preflight!(
            PREFLIGHT_STATE,
            :preview,
            0.36,
            "Compiling package preview and uncertainty geometry…"
        )
        projection = preview_projection(design)
        boundaries = boundary_projection(design)
        wire_envelope = envelope_projection(
            design,
            DEFAULT_INPUT,
            :wire_diameter
        )
        insulation_envelope = envelope_projection(
            design,
            DEFAULT_INPUT,
            :insulation_thickness
        )
        yield()

        set_preflight!(
            PREFLIGHT_STATE,
            :uncertainty,
            0.62,
            "Compiling uncertain cable constants…"
        )
        constants = compute_constants(DEFAULT_INPUT)
        WARMUP_RESULT[] = (;
            design,
            projection,
            boundaries,
            wire_envelope,
            insulation_envelope,
            constants
        )
        yield()

        set_preflight!(
            PREFLIGHT_STATE,
            :dashboard_render,
            0.82,
            "Compiling the persistent WGLMakie dashboard…"
        )
        WGLMakie.activate!()
        warm_app = App() do session
            state = setup(session)
            return DOM.div(
                build_geometry_page(session, state).body...,
                build_page(session, state).body...
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
            "Cable-constants dashboard hot"
        )
    catch error
        set_preflight!(
            PREFLIGHT_STATE,
            :error,
            PREFLIGHT_STATE[].progress,
            "Cable-constants warmup failed: $(sprint(showerror, error))"
        )
        @error "ICHQP2026 cable-constants warmup failed" exception = (
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
            "Starting cable-constants warmup…"
        )
        WARMUP_TASK[] = @async warm_cable_constants!()
        errormonitor(WARMUP_TASK[])
        return WARMUP_TASK[]
    end
end

const PREFLIGHT_RESOURCE = preflight_resource(
    id = "cable-constants-dashboard",
    title = "Cable geometry and constants under uncertainty",
    activate = start_warmup!,
    state = PREFLIGHT_STATE
)

const DECK = deck_descriptor(
    id = "cable-constants",
    group = "ICHQP2026",
    title = "Cable constants",
    order = 30,
    render = true,
    setup = setup,
    resources = (PREFLIGHT_RESOURCE,),
    pages = (
        deck_page(
        "Geometry uncertainty envelopes";
        id = "geometry-uncertainty-envelopes",
        class = "lc-application-slide lc-fill-page",
        build = build_geometry_page
    ),
        deck_page(
        "Cable constants under uncertainty";
        id = "cable-constants-under-uncertainty",
        class = "lc-application-slide lc-fill-page",
        build = build_page
    ),
    )
)

end

ICHQP2026CableConstantsDeck.DECK
