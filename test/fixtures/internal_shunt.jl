# Small physical fixtures, intentionally unrelated to catalogue tags or IDs.
function internal_shunt_test_design(; radius=0.2e-3, epsilon=2.5, tapes=true,
        count=8, suffix="", core_radius=1e-3, wire_angle=0.0)
    copper = Material(:conductor,1.72e-8,1.0,1.0,20.0,0.0039)
    dielectric = Material(:insulator,Inf,epsilon)
    fill = Material(:insulator,Inf,1.0)
    inner = core_radius+1e-3
    wire_outer = inner+2radius
    outer = wire_outer+(tapes ? 0.1e-3 : 0.0)
    wires = Group(Symbol(:middle,suffix),Region(:round,Disk(radius),copper);
        pattern=Ring(count;r=inner+radius,φ0=wire_angle))
    contents = if tapes
        tape = Group(Symbol(:middle,suffix),Region(:strip,Rectangle(0.001,0.1e-3),copper);
            pattern=Ring(1;r=wire_outer+0.05e-3))
        Stack(wires,tape)
    else
        wires
    end
    return build(CableDesign,"internal-shunt"*suffix,
        Group(Symbol(:inner,suffix),Region(:metal,Disk(core_radius),copper)),
        Region(:inside,Shell(1e-3),dielectric),
        Enclosure(:host,contents;primitive=Annulus(inner,outer),fill),
        Region(:outside,Shell(0.5e-3),dielectric),
        Group(Symbol(:reference,suffix),Region(:shield,Shell(0.1e-3),copper)),
        Region(:jacket,Shell(0.5e-3),dielectric))
end

function internal_shunt_test_domain(design)
    engine = LineCableModels.Engine
    blueprint = engine.flatten(LineCableModelsCoaxial(),design)
    return only(engine.internal_shunt_domains([design],[blueprint])),blueprint
end
