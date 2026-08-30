module CableDesignPage

using Bonito
using LineCableModels
using Makie
using Markdown
using WGLMakie
using ..PageAuthoring

const CORE_MATERIAL, INSULATION_MATERIAL = let materials = MaterialsLibrary()
    Material(materials, :copper), Material(materials, :xlpe)
end

# Frozen visual samples from the current material-preview theme. This page
# intentionally does not call or extend the package's PlotBuilder machinery.
const PLOT_BACKGROUND = RGBf(0.9, 0.9, 0.9)
const COPPER_COLOR = RGBf(0.9093547, 0.8964516, 0.8706453)
const XLPE_COLOR = RGBf(0.06460433, 0.1270558, 0.1206059)

function design(core_radius_mm::Real, insulation_thickness_mm::Real)
    isfinite(core_radius_mm) && core_radius_mm > zero(core_radius_mm) ||
        throw(DomainError(core_radius_mm, "core radius must be positive and finite"))
    isfinite(insulation_thickness_mm) &&
    insulation_thickness_mm > zero(insulation_thickness_mm) || throw(DomainError(
        insulation_thickness_mm,
        "insulation thickness must be positive and finite"
    ))

    core_radius_m = float(core_radius_mm) * 1.0e-3
    insulation_thickness_m = float(insulation_thickness_mm) * 1.0e-3
    core = LineCableModels.Conductor.Solid(
        :core_conductor,
        CORE_MATERIAL;
        r = core_radius_m
    )
    insulation = LineCableModels.Insulator.Shell(
        :insulation,
        INSULATION_MATERIAL;
        t = insulation_thickness_m
    )
    root = Stack(Group(:core, core), insulation)
    return LineCableModels.build(CableDesign, "interactive-demo", root)
end

function geometry(cable::CableDesign)
    regions = cable.geometry.regions
    core = only(filter(region -> region.source.tag === :core_conductor, regions))
    insulation = only(filter(region -> region.source.tag === :insulation, regions))
    return (
        core_r_in_m = r_in(core.primitive),
        core_r_ex_m = r_ex(core.primitive),
        core_area_m2 = area(core.primitive),
        insulation_r_in_m = r_in(insulation.primitive),
        insulation_r_ex_m = r_ex(insulation.primitive)
    )
end

function cable_figure(design_state::Observable; session = nothing)
    observe(f) = isnothing(session) ? map(f, design_state) : map(f, session, design_state)
    insulation_circle = observe() do cable
        radius_mm = geometry(cable).insulation_r_ex_m * 1.0e3
        Circle(Point2f(0), Float32(radius_mm))
    end
    core_circle = observe() do cable
        radius_mm = geometry(cable).core_r_ex_m * 1.0e3
        Circle(Point2f(0), Float32(radius_mm))
    end

    figure = Figure(size = (760, 520), backgroundcolor = PLOT_BACKGROUND)
    axis = Axis(
        figure[1, 1];
        aspect = DataAspect(),
        backgroundcolor = PLOT_BACKGROUND
    )
    insulation_plot = poly!(
        axis,
        insulation_circle;
        color = XLPE_COLOR,
        strokecolor = :black,
        strokewidth = 0.5,
        label = "XLPE insulation"
    )
    core_plot = poly!(
        axis,
        core_circle;
        color = COPPER_COLOR,
        strokecolor = :black,
        strokewidth = 0.5,
        label = "Copper core"
    )
    hidedecorations!(axis)
    hidespines!(axis)
    limits!(axis, -50, 50, -50, 50)
    axislegend(
        axis;
        position = :rt,
        backgroundcolor = (:white, 0.9),
        framecolor = RGBf(0.55, 0.55, 0.55),
        framevisible = true,
        labelcolor = RGBf(0.15, 0.15, 0.15),
        padding = (8, 8, 6, 6)
    )
    return (;
        figure,
        axis,
        insulation_plot,
        core_plot,
        insulation_circle,
        core_circle
    )
end

function default_figure()
    state = Observable(design(12.5, 8.0))
    return cable_figure(state)
end

function build(session::Session)
    core_slider = Bonito.Slider(
        collect(5.0:0.1:25.0);
        value = 12.5,
        id = "core-radius-input",
        ariaLabel = "Core radius in millimetres"
    )
    insulation_slider = Bonito.Slider(
        collect(1.0:0.1:20.0);
        value = 8.0,
        id = "insulation-thickness-input",
        ariaLabel = "Insulation thickness in millimetres"
    )
    design_state = map(design, session, core_slider.value, insulation_slider.value)
    geometry_state = map(geometry, session, design_state)
    plot = cable_figure(design_state; session)

    core_value = map(
        values -> "$(round(values.core_r_ex_m * 1.0e3; digits = 1)) mm",
        session,
        geometry_state
    )
    insulation_value = map(
        values -> "$(round((values.insulation_r_ex_m - values.insulation_r_in_m) * 1.0e3; digits = 1)) mm",
        session,
        geometry_state
    )
    outer_diameter = map(
        values -> "$(round(2values.insulation_r_ex_m * 1.0e3; digits = 1)) mm",
        session,
        geometry_state
    )
    core_area = map(
        values -> "$(round(values.core_area_m2 * 1.0e6; digits = 2)) mm²",
        session,
        geometry_state
    )
    design_type = map(
        cable -> "$(nameof(typeof(cable))){$(eltype(cable)), …}",
        session,
        design_state
    )

    controls = webpart(
        control(
            "Core radius",
            core_slider;
            value = core_value,
            id = "core-radius-control",
            class = "lc-core-radius-control"
        ),
        control(
            "Insulation thickness",
            insulation_slider;
            value = insulation_value,
            id = "insulation-thickness-control",
            class = "lc-insulation-thickness-control"
        ),
        value_list(
            "Conductor radius" => core_value,
            "Insulation thickness" => insulation_value,
            "Outer cable diameter" => outer_diameter,
            "Conductor area" => core_area
        ),
        diagnostic(design_type; prefix = "Object: ");
        kind = :controls
    )
    plot_part = webpart(
        WGLMakie.WithConfig(plot.figure; resize_to = :parent);
        kind = :plot,
        id = "cable-plot"
    )
    content = webgrid(
        [:controls :plot];
        columns = "minmax(16rem, 3fr) minmax(0, 7fr)",
        compact_columns = "minmax(15rem, 2fr) minmax(0, 5fr)",
        rows = "minmax(0, 1fr)",
        height = "100%",
        stack_rows = "auto 26rem",
        class = "lc-interactive-grid",
        controls,
        plot = plot_part
    )
    return (;
        body = slide(
            "Core and insulation",
            content;
            lede = md"""
            Adjust the dimensions to reconstruct the domain object and update the cross-section in place.
            """
        ),
        design_state,
        geometry_state,
        plot,
        core_slider,
        insulation_slider
    )
end

const PAGE = page_descriptor(
    id = "core-and-insulation",
    group = "Cable design",
    title = "Core and insulation",
    order = 20,
    render = true,
    class = "lc-interactive-slide",
    build = build,
    export_figure = default_figure
)

end

CableDesignPage.PAGE
