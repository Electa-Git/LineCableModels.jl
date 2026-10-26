@testset "Engine matrices inside the cables" begin
    Engine = LineCableModels.Engine
    reductions = LineCableModels.FormulationOptions(CoaxialFormulation)
    copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
    # The wire cases take the copper of the materials library that their scripts read.
    library_copper = Material(:conductor, 1.7241e-8, 1, 0.999994, 20, 0.00393)
    insulation = Material(:insulator, 1.97e14, 2.3)
    polymer = Material(:insulator, 1.0e14, 2.3)
    semicon = Material(:semicon, 1000.0, 1000.0)
    steel = Material(:conductor, 1.7475e-7, 1.0, 36.6326, 20.0, 0.0045)
    air = material(kind = :insulator, rho = Inf)
    frequencies = 10.0 .^ range(-1, 7; length = 101)
    metal = Terminal(:core, Region(:core_metal, Circle(0.0425), library_copper))
    wires = [
        "two_bare_wires" => build(CableDesign, "bare_wire", metal),
        "two_insulated_wires" => build(CableDesign, "insulated_wire", metal,
            Region(:core_insulation, Offset(1.0e-3), insulation))]
    cases = [case => LineParametersProblem(
                 build(CableSystem, [design, design],
                     [Placement(0.0, -1.0), Placement(1.0, -1.0)];
                     connections = [(core = 1,), (core = 2,)], id = case, line_length = 1.0),
                 homogeneous(rho = 0.1); frequencies)
             for (case, design) in wires]
    for (case, strand) in (("coaxial_circular_strands", Circle(1.5e-3)),
        ("coaxial_rectangular_strands", Rectangle(3.0e-3, 3.0e-3)))
        courses = [Lay(Helix(LayRatio(12)),
                       Repeat(CircularPattern(n; r = 3.0e-3k),
                           Region(Symbol(:core_strands_, k), strand, copper)))
                   for (k, n) in enumerate((6, 12, 18))]
        armour = Lay(Helix(LayRatio(10)),
            Repeat(CircularPattern(40; r = 26.5e-3),
                Region(:armour_wires, Circle(2.0e-3), steel)))
        design = build(CableDesign, "coaxial",
            Terminal(:core,
                Enclosure(Region(:core_fill, Circle(10.5e-3), air),
                    Stack(Region(:core_centre, Circle(1.5e-3), copper), courses...))),
            Region(:core_semicon_inner, Offset(1.0e-3), semicon),
            Region(:core_insulation, Offset(8.0e-3), polymer),
            Region(:core_semicon_outer, Offset(1.0e-3), semicon),
            Terminal(:sheath, Region(:sheath_copper, Annulus(20.5e-3, 21.5e-3), copper)),
            Region(:sheath_bedding, Offset(3.0e-3), polymer),
            Terminal(:armour, Enclosure(Region(:armour_fill, Offset(4.0e-3), air), armour)),
            Region(:armour_jacket, Offset(4.0e-3), polymer))
        system = build(CableSystem, [design], [Placement(0.0, -1.0)];
            connections = [(core = 1, sheath = 2, armour = 3)], id = case,
            line_length = 1000.0)
        push!(cases,
            case => LineParametersProblem(system, homogeneous(rho = 100, eps_r = 10);
                frequencies))
    end
    internal = Engine.InternalImpedance.Formula(Val(:schelkunoff1934))
    term = Engine.InsulationImpedance.Formula(Val(:ametani1980))
    shunt = Engine.Shunt.Formula(Val(:series))
    law = Engine.Resistivity.Formula(Val(:linear))
    insulating = Engine.InsulationAdmittance.Formula(Val(:lossless))
    semiconducting = Engine.SemiconAdmittance.Formula(Val(:lossless))

    @testset "$case against its stage targets" for (case, problem) in cases
        input = Engine.CoaxialInput(problem, reductions)
        n = length(input.cable)
        c = length(input.radius)
        targets = Dict(matrix => zeros(ComplexF64, rows, rows, length(frequencies))
        for (matrix, rows) in
            (("Zin", n), ("Pin", n), ("Zg", c), ("Pg", c), ("Z", n), ("Y", n)))
        for line in eachline(joinpath(@__DIR__, "..", "references", case, "stages.csv"))
            fields = split(line, ',')
            startswith(line, "analytical") && haskey(targets, fields[3]) || continue
            k = searchsortedfirst(frequencies, parse(Float64, fields[2]))
            targets[fields[3]][parse(Int, fields[4]), parse(
                Int, fields[5]), k] = complex(
                parse(Float64, fields[6]), parse(Float64, fields[7]))
        end
        for line in eachline(joinpath(@__DIR__, "..", "references", case, "matrices.csv"))
            fields = split(line, ',')
            startswith(line, "analytical") || continue
            k = searchsortedfirst(frequencies, parse(Float64, fields[2]))
            targets[fields[3]][parse(Int, fields[4]), parse(
                Int, fields[5]), k] = complex(
                parse(Float64, fields[6]), parse(Float64, fields[7]))
        end
        rho = [Engine.resistivity(law, tube, problem.temperature)
               for tube in input.material]
        semicons = findall(layer -> layer isa Material{Float64, :semicon},
            input.layer_material)
        Zin = Array{ComplexF64}(undef, n, n, length(frequencies))
        Pin = similar(Zin)
        for (k, f) in enumerate(frequencies)
            s = im * (2π * f)
            kappa = ComplexF64[Engine.admittivity(
                                   i in semicons ? semiconducting : insulating, layer, f)
                               for (i, layer) in enumerate(input.layer_material)]
            Engine.impedance!(internal, term, view(Zin,:,:,k), input, rho, s)
            Engine.potential!(shunt, view(Pin,:,:,k), input, kappa, s)
        end
        for (computed, target) in ((Zin, targets["Zin"]), (Pin, targets["Pin"]))
            @test all(maximum(abs, computed[:, :, k] - target[:, :, k]) <=
                      1.0e-12 * maximum(abs(target[i, i, k]) for i in 1:n)
            for k in eachindex(frequencies))
        end
        result = compute(problem, CoaxialFormulation(); trace = true)
        for (matrix, computed, tolerance) in (("Z", result.Z, 1e-9), ("Y", result.Y, 1e-9),
            ("Zin", result.trace.Zin, 1e-12), ("Pin", result.trace.Pin, 1e-12),
            ("Zg", result.trace.Zg, 1e-9), ("Pg", result.trace.Pg, 1e-9))
            target = targets[matrix]
            errors = [maximum(abs, computed[:, :, k] - target[:, :, k])
                      for k in eachindex(frequencies)]
            scales = [maximum(abs(target[i, i, k]) for i in axes(target, 1))
                      for k in eachindex(frequencies)]
            @test all(errors .<= tolerance .* scales)
        end
    end
end
