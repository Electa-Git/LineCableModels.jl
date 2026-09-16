@testitem "TextDisplay / physical declarations retain units and placement provenance" tags=[:unit] begin
    const DM = LineCableModels.DataModel
    const TD = LineCableModels.TextDisplay
    copper = Material(:conductor, 1.7241e-8)
    dielectric = Material(:insulator, Inf, 2.3)
    pose = Pose2(2e-3, -3e-3, pi / 6)
    sector = Sector(span=2pi / 3, r_base=0.6e-3, r_back=5e-3, fillet=0.2e-3)
    member = terminal(:core, solid(copper, sector), insulation(dielectric; t=0.2e-3))
    design = build(CableDesign, "display-sector", member)
    placed = at(member, pose)
    explicit = assembly(member, placed)
    repeated = Assembly(member; at=pose, pattern=Ring(3; r=8e-3),
        names=(:a, :b, :c), path=Helix(Pitch(0.25)), compact=FillFactor(0.9))
    wire = Region(:wire, Disk(0.5e-3), copper)
    patterned = Group(:core, wire; at=pose, pattern=Ring(6; r=1e-3),
        path=Helix(LayAngle(pi / 12)), compact=FillFactor(0.9))
    bounded = stranded(copper; shape=Disk(0.5e-3), boundary=Disk(1.5e-3), compact=true)
    enclosure = Enclosure(:duct, member; at=pose, primitive=Disk(10e-3),
        fill=Region(:bedding, Disk(10e-3), dielectric),
        wall=insulation(dielectric; t=1e-3))
    catalogue_record = DatasheetInfo(U0=76.0, U=132.0,
        conductor_cross_section=630.0, screen_cross_section=35.0,
        armor_cross_section=50.0, resistance=0.0283, capacitance=0.2,
        inductance=0.3, description="Reference construction")
    primitives = (
        Disk(1e-3), DM.resolve(pose, Disk(1e-3)), Rectangle(2e-3, 1e-3),
        Ellipse(3e-3, 2e-3), Annulus(1e-3, 2e-3), sector, Shell(0.2e-3),
        Ring(capacity()), Ring(6; r=1e-3, φ0=pi/6, span=pi, gap_frac=0.1),
        Polar(nr=2, nφ=3, r0=1e-3, dr=2e-3, φ0=pi/6, span=pi),
        Fill(r=3e-3, φ=pi/6, φ0=pi/12, span=pi),
        Lattice(nx=2, ny=3, dx=1e-3, dy=2e-3),
        FillFactor(0.9), LayRatio(12), Pitch(0.25), LayAngle(pi/12),
        Helix(LayRatio(12); dir=-1, φ0=pi/6),
        Helix(Pitch(0.25); dir=-1, φ0=pi/6),
        Helix(LayAngle(pi/12); dir=-1, φ0=pi/6),
        DM.EmptyBoundary(), DM.EnclosureBoundary(), catalogue_record,
    )
    physical_objects = (
        wire, Stack(wire), patterned, Group(:single, wire), bounded,
        placed, explicit, repeated, enclosure, design, design.geometry,
        design.geometry.regions..., (region.primitive for region in design.geometry.regions)...,
        DM.preview_shapes(first(design.geometry.regions))...,
        CablesLibrary(),
    )
    for object in (primitives..., physical_objects...)
        @test !isempty(TD.name(typeof(object)))
        @test !isempty(sprint(summary, object))
        compact = sprint(show, object)
        @test !occursin('\n', compact)
        @test sprint(show, MIME"text/plain"(), object; context=:compact=>true) == compact
        @test !endswith(sprint(show, MIME"text/plain"(), object), '\n')
    end

    # These are scientific descriptions, not dumps of implementation type
    # parameters. Check the quantities and identities, not cosmetic whitespace.
    @test occursin("mm", sprint(show, Disk(1e-3)))
    @test occursin("capacity()", sprint(show, Ring(capacity())))
    @test occursin("at=", sprint(show, DM.resolve(pose, Disk(1e-3))))
    @test occursin("−1", sprint(show, Helix(Pitch(0.25); dir=-1)))
    @test occursin("mm²", sprint(show, catalogue_record))
    @test occursin("Ω/km", sprint(show, catalogue_record))
    @test occursin("μF/km", sprint(show, catalogue_record))
    @test occursin("mH/km", sprint(show, catalogue_record))
    @test occursin("kV", sprint(show, catalogue_record))
    @test occursin("Reference construction", sprint(show, catalogue_record))
    @test occursin("2 explicit members", sprint(show, MIME"text/plain"(), explicit))
    @test occursin("250 mm", sprint(show, MIME"text/plain"(), repeated))
    @test occursin("FillFactor", sprint(show, MIME"text/plain"(), patterned;
        context=IOContext(IOBuffer(), :displaysize=>(40, 200))))
    @test occursin("fill · Region :bedding", sprint(show, MIME"text/plain"(), enclosure))
    @test occursin("wall", sprint(show, MIME"text/plain"(), enclosure))
    @test occursin("terminal=:core", sprint(show, first(design.geometry.regions)))
    @test !occursin("terminal=", sprint(show, last(design.geometry.regions)))
    @test occursin("inner", sprint(show, MIME"text/plain"(), last(design.geometry.regions).primitive))
    @test occursin("outer", sprint(show, MIME"text/plain"(), last(design.geometry.regions).primitive))
end
