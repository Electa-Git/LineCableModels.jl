"""
    CableGeometry()

Construct a session-local coaxial cross-section with core-radius, insulation and
sheath thickness inputs in mm. Circles preserve radial proportions; their total
display radius is fitted to the frame. This illustration does not submit jobs,
modify the scientific views' inputs, or infer electrical/material properties.
"""
struct CableGeometry{F}
    "Existing typed numeric fields; instantiate once per owning session."
    fields::F
end

function CableGeometry()
    id = "geometry-$(uuid4())"
    return CableGeometry((core=number(id, :core, "Core radius", 10, 2, 30, 1, "mm"),
        insulation=number(id, :insulation, "Insulation thickness", 5, 1, 20, 1, "mm"),
        sheath=number(id, :sheath, "Sheath thickness", 1, 0.1, 4, 0.1, "mm")))
end

function Bonito.jsrender(session::Session, view::CableGeometry)
    core, insulation, sheath = (f.control.value for f in view.fields)
    radius = map(session, core, insulation, sheath) do c, i, s
        all(isfinite, (c,i,s)) && 2 <= c <= 30 && 1 <= i <= 20 && 0.1 <= s <= 4 || return (0.0,0.0)
        (140*c/(c+i+s), 140*(c+i)/(c+i+s))
    end
    diagram = SVG.svg(SVG.title("Coaxial cable radial construction"),
        SVG.circle(; cx=360, cy=175, r=140, class="lc-study-sheath"),
        SVG.circle(; cx=360, cy=175, r=last(radius[]), class="lc-study-insulation"),
        SVG.circle(; cx=360, cy=175, r=first(radius[]), class="lc-study-core");
        viewBox="0 0 720 350", role="img", class="lc-study-plot")
    Bonito.onload(session, diagram, js"""
        element => {
            const radius = $(radius);
            const update = value => {
                element.querySelector('.lc-study-core').setAttribute('r', String(value[0]));
                element.querySelector('.lc-study-insulation').setAttribute('r', String(value[1]));
            };
            radius.on(update);
            update(radius.value);
        }
    """)
    properties = PropertyGrid(PropertyItem("Conductor", "Core and sheath"),
        PropertyItem("Dielectric", "Insulation annulus"),
        PropertyItem("Electrical evaluation", "Separate worker operation"))
    inputs = ViewportFrame("Construction inputs", DOM.div(
        DOM.div(view.fields...; class="lc-form-fields lc-study-fields"), properties;
        class="lc-content-stack lc-panel-content"); sizing=:content)
    node = DOM.section(WorkspacePage("Cable construction",
        SplitPane(ViewportFrame("Core · insulation · sheath", diagram; sizing=:fill), inputs;
            ratio=.68, min_first="20rem", min_second="17rem", scroll=:parent);
        eyebrow="SCENE VIEWPORT", fill=true,
        description="Local geometry interaction only. Radii remain proportional. Scientific inputs are configured separately in the line-parameter view.");
        class="lc-study-geometry")
    return Bonito.jsrender(session, DOM.div(styles(), ComponentXRay.instrument(session, node, view); style="display: contents;"))
end

function ComponentXRay.inspection(view::CableGeometry)
    return ComponentXRay.ComponentInspection(view; name="CableGeometry",
        source=ComponentXRay.source_reference(@__MODULE__, @__FILE__, @__LINE__),
        css_scopes=[".lc-study-geometry",
            ".lc-study-plot", ".lc-study-sheath", ".lc-study-insulation", ".lc-study-core", ".lc-study-note"],
        notes=["Pure radial display. Typed fields own their bindings; no scientific execution."])
end
