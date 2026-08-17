function _fixture_sha(name::AbstractString)
    return bytes2hex(SHA.sha256(name))
end

"""Build the native EMT case used by Tutorial 3 without presentation or file I/O."""
function performance_case(::Val{:tutorial3})
    materials = MaterialsLibrary(add_defaults = true)
    copper = Material(materials, :copper)
    semicon1 = Material(materials, :semicon1)
    semicon2 = Material(materials, :semicon2)
    pe = Material(materials, :pe)
    polyacrylate = Material(materials, :polyacrylate)
    lead = Material(materials, :lead)
    pp = Material(materials, :pp)
    steel = Material(materials, :steel)

    core = (
        Conductor.Stranded(
            :core;
            layers = 7,
            wire_radius = 3.6649e-3 / 2,
            num_wires = 6,
            lay_ratio = 11.0,
            material = copper
        ),
        Insulator.Semicon(:core; thickness = 2e-3, material = semicon1),
        Insulator.Tubular(:core; thickness = 26e-3, material = pe),
        Insulator.Semicon(:core; thickness = 1.8e-3, material = semicon2),
        Insulator.Semicon(:core; thickness = 0.3e-3, material = polyacrylate)
    )
    sheath = (
        Conductor.Tubular(:sheath; thickness = 3.3e-3, material = lead),
        Insulator.Tubular(:sheath; thickness = 3e-3, material = pe),
        Insulator.Tubular(:sheath; thickness = 3e-3, material = pp)
    )
    armor = (
        Conductor.Wires(
            :armor;
            wire_radius = 5.827e-3 / 2,
            num_wires = 68,
            lay_ratio = 10.0,
            material = steel
        ),
        Insulator.Tubular(:armor; thickness = 10e-3, material = pp)
    )
    design = only(CableBuilder("tutorial3-525kv", core, sheath, armor))
    positions = (
        at(x = -0.5, y = -1.0, phases = (:core => 1, :sheath => 0, :armor => 0)),
        at(x = 0.5, y = -1.0, phases = (:core => 2, :sheath => 0, :armor => 0))
    )
    problem = only(SystemBuilder(
        "tutorial3-525kv-bipole",
        design,
        positions;
        length = 1000.0,
        temperature = 20.0,
        earth = Earth(rho = 100.0, eps_r = 10.0, mu_r = 1.0),
        frequencies = collect(10.0 .^ range(0, stop = 6, length = 61))
    ))
    identity = _fixture_sha("LineCableModels Tutorial 3 EMT performance case")
    provenance_value = Provenance(
        :linecablemodels,
        string(Base.pkgversion(LineCableModels)),
        "Tutorial 3",
        "525 kV 1600 mm² armored bipole",
        identity,
        identity,
        identity,
        nothing,
        Dict{String, String}(),
        identity
    )
    return Case(
        "tutorial3-emt",
        Coax(),
        Exact(),
        Reduced(),
        problem,
        Formulation(),
        Reference(),
        (),
        Assumption[],
        provenance_value
    )
end
