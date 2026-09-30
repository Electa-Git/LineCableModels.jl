# Passive physical fixtures for the remaining conductor-mesh qualification.
# Include this file and call a constructor explicitly. No meshing, solving,
# plotting, acceptance policy or reference calculations run on inclusion.
# Dimensions and deliberate simplifications are recorded in fixtures.md.
module ConductorMeshFixtures

using LineCableModels

function screened_cable()
    materials = MaterialsLibrary(add_defaults=true)
    aluminum = Material(materials, :aluminum)
    copper = Material(materials, :copper)
    pe = Material(materials, :pe)
    wire_screen = terminal(:screen, wires(copper; shape=Disk(0.475e-3),
        n=49, r=29.025e-3, tag=:screen_wire))
    # The insulation, screen fill and bedding are all PE. Keep one material
    # region, without artificial tangent interfaces at the wire envelope.
    return build(CableDesign, "mesh-screen-49", Stack(
        Enclosure(:insulation_and_screen,
            Stack(terminal(:core,core(aluminum;r=19.05e-3)),wire_screen);
            primitive=Disk(29.9e-3),fill=pe),
        terminal(:foil, sheath(aluminum; t=0.15e-3)),
        jacket(pe; t=2.45e-3)))
end

function tubular_cable()
    materials = MaterialsLibrary(add_defaults=true)
    copper = Material(materials, :copper)
    lead = Material(materials, :lead)
    pe = Material(materials, :pe)
    return build(CableDesign, "mesh-lead-sheath", Stack(
        terminal(:core, core(copper; r=23.15e-3)),
        insulation(pe; t=30.1e-3),
        terminal(:sheath, sheath(lead; t=3.3e-3)),
        jacket(pe; t=3e-3)))
end

function sector_cable()
    materials = MaterialsLibrary(add_defaults=true)
    aluminum = Material(materials, :aluminum)
    copper = Material(materials, :copper)
    pvc = Material(kind=:insulator, rho=Inf, eps_r=8.0, mu_r=1.0)
    shape = Sector(span=deg2rad(119.0), r_base=1.10e-3,
        r_back=10.24e-3, fillet=1.02e-3)
    # Sleeves, bedding, neutral fill and jacket are the same PVC. One fill
    # avoids overlapping sleeves and artificial seams tangent to the wires.
    phases = assembly((at(terminal(name, core(aluminum, shape)),
        0.0, 0.0; φ=angle) for (name, angle) in
        zip((:a, :b, :c), (0.0, 2pi/3, 4pi/3)))...)
    neutral = terminal(:neutral, wires(copper; shape=Disk(0.79e-3),
        n=30, r=14.36e-3, tag=:neutral_wire))
    return build(CableDesign, "mesh-three-sector",
        Enclosure(:pvc, Stack(phases, neutral);
            primitive=Disk(17.25e-3), fill=pvc))
end

end
