# Current local-domain inputs. The open screen permits testing wire/tape charge
# scattering; the independent annular operator controls are in internal_shunt.jl.
function internal_shunt_test_design(;
        radius = 0.0003, epsilon = 3.0, tapes = true, count = 6,
        suffix = "", core_radius = 0.002, wire_angle = 0.17)
    metal=Material(kind = :conductor, rho = 2e-8, mu_r = 1.0, T0 = 20.0, alpha = 0.004)
    insulation_material=Material(kind = :insulator, rho = Inf, eps_r = epsilon)
    host_material=Material(kind = :insulator, rho = Inf, eps_r = 1.0)
    a=core_radius+0.002
    screen_outer=a+2radius
    t=0.00015
    b=screen_outer+(tapes ? t : zero(radius))
    open_parts=AbstractCablePart[terminal(Symbol(:middle, suffix),
        Group(:round_screen, Region(:wire, Disk(radius), metal);
            pattern = Ring(count; r = a+radius, φ0 = wire_angle)))]
    if tapes
        push!(open_parts,
            terminal(Symbol(:middle, suffix),
                Group(:strip_screen, Region(:strip, Rectangle(0.0012, t), metal);
                    pattern = Ring(1; r = screen_outer+t/2))))
    end
    origin=Stack(
        terminal(Symbol(:inner, suffix), Region(:axial, Disk(core_radius), metal)),
        Region(:inner_dielectric, Shell(0.002), insulation_material),
        Enclosure(:screen_domain, Stack(open_parts); primitive = Annulus(a, b), fill = host_material),
        Region(:outer_dielectric, Shell(0.0008), insulation_material),
        terminal(Symbol(:reference, suffix), Region(:closed_screen, Shell(0.0002), metal)),
        Region(:cover, Shell(0.0007), insulation_material))
    return build(CableDesign, "open-screen-control"*suffix, origin)
end
function internal_shunt_test_domain(design)
    owner=LineCableModels.Engine
    blueprint=owner.flatten(LineCableModelsCoaxial(), design)
    domains=owner.ShuntModel.internal_shunt_domains([design], [blueprint])
    length(domains)==1 || error("control requires exactly one resolved open-screen domain")
    return only(domains), blueprint
end
