@testmodule CoatedArmour begin
    using LineCableModels

    # The cable has a core with its insulation and a lead sheath with its PE jacket. Then come
    # `n` steel armour wires of radius 1.5 mm, each with a coat of thickness `coat`, and a PP
    # jacket. The ring lies in a matrix fill that starts at its inner side and whose half-width
    # is `span`, the coated wire's radius by default. `coat=nothing` gives bare wires, and
    # `gap` adds matrix shells of that thickness around a bare ring, which then lies at the
    # pitch of the coated one. With `enclosed=false`, the ring takes its pitch from its context
    # and has no fill. All lengths are in m.
    function armoured(; coat = 0.2e-3, coat_material = Material(:insulator, 1.0e14, 2.3),
            n = 36, span = nothing, gap = nothing, enclosed = true, T = Float64)
        steel = Material(:conductor, T(1.38e-7), T(1.0), T(10.0), T(20.0), T(0.0045))
        lead = Material(:conductor, T(2.14e-7), T(1.0), T(1.0), T(20.0), T(0.004))
        copper = Material(:conductor, T(1.72e-8), T(1.0), T(1.0), T(20.0), T(0.00393))
        xlpe = Material(:insulator, T(1.0e14), T(2.4))
        pe = Material(:insulator, T(1.97e14), T(2.5))
        pp = Material(:insulator, T(1.0e15), T(2.8))
        matrix = Material(kind = :insulator, rho = T(Inf), eps_r = T(1.0), mu_r = T(1.0))
        wire = T(1.5e-3)
        radius = T(21.5e-3)
        parts = Any[Group(:core, Region(:core_conductor, Disk(T(10e-3)), copper)),
            Region(:insulation, Shell(T(8e-3)), xlpe),
            Group(:sheath, Region(:sheath_lead, Annulus(T(18e-3), T(19.5e-3)), lead)),
            Region(:sheath_pe, Shell(T(2e-3)), pe)]
        if gap !== nothing
            push!(parts, Region(:gap_inner, Shell(T(gap)), matrix))
            radius += T(gap)
        end
        member = coat === nothing ? Region(:armor_wire, Disk(wire), steel) :
                 Stack(Region(:armor_wire, Disk(wire), steel),
            Region(:armor_coat, Shell(coat), coat_material))
        path = Helix(LayRatio(T(10)))
        if enclosed
            half = coat === nothing ? wire : wire + coat
            armor = Group(:armor, member; pattern = Ring(n; r = radius + half), path)
            push!(parts, Enclosure(:armor_matrix, armor;
                primitive = Annulus(radius, radius + 2something(span, half)), fill = matrix))
        else
            push!(parts, Group(:armor, member; pattern = Ring(n; r = nothing), path))
        end
        gap === nothing || push!(parts, Region(:gap_outer, Shell(T(gap)), matrix))
        push!(parts, Region(:armor_jacket, Shell(T(3e-3)), pp))
        return build(CableDesign, "coated-armour", Stack(parts...))
    end

    # One cable 1 m deep in a homogeneous earth, each terminal a phase.
    function problem(design; frequencies = [50.0, 1.0e4])
        system = build(LineCableSystem, design, Pose2(0.0, -1.0);
            connections = Dict(:core => 1, :sheath => 2, :armor => 3))
        return LineParametersProblem(system; earth_props = homogeneous(rho = 100.0),
            frequencies)
    end
end
