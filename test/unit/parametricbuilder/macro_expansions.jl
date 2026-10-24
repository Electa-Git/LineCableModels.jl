# Each expected expression is the expansion recorded before the node checks of the
# construction macros became methods. Inputs are parsed under the file name "golden".
@testitem "ParametricBuilder / construction API / macro expansions are unchanged" tags=[
    :unit, :parametric
] begin
    PB = LineCableModels.ParametricBuilder
    DM = LineCableModels.DataModel
    expand(code) = macroexpand(LineCableModels,
        last(Meta.parseall(code; filename = "golden").args))
    build = GlobalRef(PB, :build)
    at = GlobalRef(PB, :at)

    # Line numbers, docstrings and `nothing` are not parts. A documented item is its expression.
    @test expand("""
        @cable identifier begin
            "documented"
            core(1)
            nothing
            @doc "explicit" sheath(2)
            "plain text"
        end""") == :($build($(GlobalRef(DM, :CableDesign)), identifier, core(1), sheath(2)))
    @test expand("""
        @system identifier line_length = 1000.0 begin
            @at design (0.0, -1.0)
            "documented"
            @trefoil design spacing = 0.1
        end""") == :($build($(GlobalRef(DM, :LineCableSystem)),
        ($at(design, 0.0, -1.0;), $(GlobalRef(PB, :trefoil))(design; spacing = 0.1));
        system_id = identifier, line_length = 1000.0))
    @test expand("""
        @earth vertical_layers = false begin
            layer(rho = 100.0)
            nothing
        end""") ==
          :($build($(GlobalRef(LineCableModels.Earth, :EarthModel)), (layer(rho = 100.0),);
        vertical_layers = false))
    @test expand("""
        @terminal core begin
            "documented"
            conductor
        end""") == :($(GlobalRef(PB, :terminal))(core, conductor))
    @test expand("""
        @assembly combine = :zip begin
            first_member
            nothing
            second_member
        end""") == :($(GlobalRef(PB, :assembly))(first_member, second_member; combine = :zip))
    @test expand("""
        @pipe shape = disk fill = air begin
            "documented"
            member
        end""") == :($(GlobalRef(PB, :pipe))(member; shape = disk, fill = air))
    @test expand("""
        @duct shape = disk fill = air wall = steel begin
            member
        end""") == :($(GlobalRef(PB, :duct))(member; shape = disk, fill = air, wall = steel))

    # A center given as a coordinate tuple becomes a placement. Any other center is unchanged.
    trefoil = GlobalRef(PB, :trefoil)
    @test expand("@trefoil design spacing = 0.1") == :($trefoil(design; spacing = 0.1))
    @test expand("@trefoil design spacing = 0.1 center = (0.0, -1.0) core = 1") ==
          :($trefoil(design; spacing = 0.1, center = $at(0.0, -1.0), connections = (; core = 1)))
    @test expand("@trefoil design spacing = 0.1 center = (0.0, -1.0, 30.0) φ0 = 10.0") ==
          :($trefoil(design; spacing = 0.1, center = $at(0.0, -1.0; φ = 30.0), φ0 = 10.0))
    @test expand("@hflat design spacing = 0.2 center = origin connections = (core = 1,)") ==
          :($(GlobalRef(PB, :hflat))(design; spacing = 0.2, center = origin,
        connections = (core = 1,)))
    @test expand("@vflat design spacing = 0.2 center = (x, y)") ==
          :($(GlobalRef(PB, :vflat))(design; spacing = 0.2, center = $at(x, y)))
end
