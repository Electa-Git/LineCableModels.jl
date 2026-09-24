# Make the disposable script runnable directly, without requiring `--project=.`.
pushfirst!(LOAD_PATH, normpath(joinpath(@__DIR__, "..", "..", "..")))

using LineCableModels
using GLMakie

# For all you unbelievers: a six-terminal, genuinely surprised electric mouse.
const mm = 1.0e-3

function ellipse_polygon(cx, cy, rx, ry; vertices = 48, rotation = 0.0)
    cosine, sine = cos(rotation), sin(rotation)
    return [begin
                x, y = rx * cos(angle), ry * sin(angle)
                ((cx + cosine * x - sine * y) * mm,
                    (cy + sine * x + cosine * y) * mm)
            end
            for angle in range(0, 2pi; length = vertices + 1)[1:(end - 1)]]
end

polygon(points) = Polygon([(x * mm, y * mm) for (x, y) in points])

# The constitutive properties stay physical; the cartoon palette below is only
# presentation metadata applied to the native preview series.
yellow_rubber = Material(:insulator, 1.0e14, 3.2, 1.0, 20.0, 0.0)
copper = Material(:conductor, 1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)

# One concave arbitrary polygon supplies the insulating outer contour. The ear
# tips and indented cheeks deliberately make this much less cable-like.
body = Region(:pikachu,
    polygon([
        (-49, -36), (49, -36), (46, -17), (42, 5), (36, 29),
        (34, 38), (48, 52), (61, 73), (34, 62), (13, 51),
        (6, 50), (-8, 50), (-16, 52), (-38, 66), (-62, 80),
        (-53, 52), (-36, 31), (-33, 12), (-39, -12), (-45, -29)
    ]),
    yellow_rubber)

# Every facial feature is also an arbitrary polygon, wrapped in its own
# terminal so this ridiculous cross-section remains a valid cable design.
left_cheek = terminal(:left_cheek,
    Region(:left_cheek, Polygon(ellipse_polygon(-27.5, -5, 8.5, 8.2)), copper))
right_cheek = terminal(:right_cheek,
    Region(:right_cheek, Polygon(ellipse_polygon(27.5, -5, 8.5, 8.2)), copper))

left_eye = terminal(:left_eye,
    assembly(
        Region(:eye, Polygon(ellipse_polygon(-18, 18, 7.0, 7.4)), copper),
        Region(:eye_highlight, Polygon(ellipse_polygon(-20.2, 20.4, 2.5, 2.7)), copper)
    ))
right_eye = terminal(:right_eye,
    assembly(
        Region(:eye, Polygon(ellipse_polygon(18, 18, 7.0, 7.4)), copper),
        Region(:eye_highlight, Polygon(ellipse_polygon(15.8, 20.4, 2.5, 2.7)), copper)
    ))

nose = terminal(:nose,
    Region(:nose, Polygon(ellipse_polygon(0, 7.8, 2.5, 1.5)), copper))
mouth = terminal(:mouth,
    assembly(
        Region(:mouth_cavity, Polygon(ellipse_polygon(0, -13, 9.0, 9.5)), copper),
        Region(:mouth_inside, Polygon(ellipse_polygon(0.5, -13.8, 7.0, 7.2)), copper)
    ))

pikachu_assembly = assembly(
    body,
    left_cheek,
    right_cheek,
    left_eye,
    right_eye,
    nose,
    mouth
)

design = build(CableDesign, "surprised-pikachu", pikachu_assembly)

presentation_groups = Dict(
    :pikachu => :insulating_body,
    :left_cheek => :cheeks,
    :right_cheek => :cheeks,
    :left_eye => :eyes,
    :right_eye => :eyes,
    :eye => :eyes,
    :eye_highlight => :eye_highlights,
    :nose => :nose,
    :mouth_cavity => :mouth_cavity,
    :mouth_inside => :mouth_inside
)

# Group order follows first appearance in the assembly: body, cheeks, eyes,
# highlights, nose, mouth cavity, mouth interior. Colors are deliberately
# presentation-only overrides.
preview_plot = preview(
    design;
    backend = :gl,
    display_plot = true,
    controls = true,
    display_legend = false,
    display_colorbars = false,
    display_dielectric_pattern = false,
    title = "Unbelievers, lo and behold!",
    panel_titles = ("I told y'all that this could run a Pikachu-shaped cable model!",),
    size = (900, 760),
    legend_group = presentation_groups,
    series_attributes = (
        (; color = "#F6C945", strokecolor = "#2A2510", strokewidth = 5.0),
        (; color = "#E64A2E", strokecolor = "#C83A21", strokewidth = 2.0),
        (; color = "#211D08", strokecolor = "#211D08", strokewidth = 1.5),
        (; color = "#FFFDF4", strokecolor = "#FFFDF4", strokewidth = 0.5),
        (; color = "#5B0805", strokecolor = "#3B0503", strokewidth = 1.5),
        (; color = "#211D08", strokecolor = "#211D08", strokewidth = 2.0),
        (; color = "#EF7C5B", strokecolor = "#EF7C5B", strokewidth = 1.0)
    )
)

# axis = only(preview_plot.axes)
# hidedecorations!(axis)
# hidespines!(axis)

# output_path = joinpath(get(ENV, "LINECABLEMODELS_MANUAL_OUTPUT",
#     joinpath(tempdir(), "linecablemodels-manual")), "plotting", "surprised_pikachu_cable.png")
# CairoMakie.save(output_path, preview_plot.figure; px_per_unit = 2)
# println("Refactor morale conductor energized: $output_path")
