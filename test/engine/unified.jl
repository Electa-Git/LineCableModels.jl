using LinearAlgebra: I, diag, lu!
using SpecialFunctions: besseli

# The targets come from the Unified formula of commit 989af96c of this repository: the state
# of the cables and the pair terms. It ran with Γ = 0 at 0.1 Hz, 1 kHz and
# 10 MHz. The earth-return terms read only the centers and the outer radii of the cables. A
# wire of the outer radius of a reference cable stands for that cable. The formulation fills
# the buffer `gamma` with its option Γ at each frequency. Here the test gives that buffer.
@testset "Engine Unified earth return" begin
    Engine = LineCableModels.Engine
    EarthReturn = Engine.EarthReturn
    formula = EarthReturn.Formula(Val(:unified))
    law = LineCableModels.Engine.EarthProperties.Formula(Val(:constant))
    reductions = LineCableModels.FormulationOptions(CoaxialFormulation)
    copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
    # Each case gives the outer radius of each cable in m, the placements and the earth. A state
    # row gives the frequency index, the cable, x, A and F. A term row gives the frequency index,
    # the receiver, the source, k in Ω/m and h in m/F.
    cases = [
        (name = "two buried bare wires", radii = [0.0425, 0.0425],
            placements = [Placement(0.0, -1.0), Placement(1.0, -1.0)],
            earth = homogeneous(rho = 0.1),
            state = [
                (1, 1, complex(8.444400025924338e-05, 8.444400025929034e-05),
                    complex(0.999915559565035, 3.5650935264235386e-09),
                    complex(0.05674501730546565, -1.0112762017883535e-10)),
                (1, 2, complex(8.444400025924338e-05, 8.444400025929034e-05),
                    complex(0.999915559565035, 3.5650935264235386e-09),
                    complex(0.05674501730546565, -1.0112762017883535e-10)),
                (2, 1, complex(0.00844440000243753, 0.00844440004941584),
                    complex(0.9915911534807212, 3.535413735005518e-05),
                    complex(0.05674501728143213, -1.0112762011926023e-06)),
                (2, 2, complex(0.00844440000243753, 0.00844440004941584),
                    complex(0.9915911534807212, 3.535413735005518e-05),
                    complex(0.05674501728143213, -1.0112762011926023e-06)),
                (3, 1, complex(0.8444165137639753, 0.8444634920747424),
                    complex(0.4161523516456092, 0.15270107075218897),
                    complex(0.05447840661298058, -0.009557481831077413)),
                (3, 2, complex(0.8444165137639753, 0.8444634920747424),
                    complex(0.4161523516456092, 0.15270107075218897),
                    complex(0.05447840661298058, -0.009557481831077413))],
            terms = [
                (1, 1, 1, complex(9.903332164789415e-08, 1.2122712849237515e-06),
                    complex(0.015708791953322424, 0.1444738754182352)),
                (1, 2, 1, complex(9.903215146015935e-08, 8.153603164141449e-07),
                    complex(0.015708540803114394, 0.11177291535348896)),
                (1, 1, 2, complex(9.903215146015935e-08, 8.153603164141449e-07),
                    complex(0.015708540803114394, 0.11177291535348896)),
                (1, 2, 2, complex(9.903332164789415e-08, 1.2122712849237515e-06),
                    complex(0.015708791953322424, 0.1444738754182352)),
                (2, 1, 1, complex(0.0011680804588461839, 0.006101515081950891),
                    complex(144.79464254274615, 533.9355326013214)),
                (2, 2, 1, complex(0.001112780868788763, 0.002112876932097273),
                    complex(137.72928290305848, 207.11125922198644)),
                (2, 1, 2, complex(0.001112780868788763, 0.002112876932097273),
                    complex(137.72928290305848, 207.11125922198644)),
                (2, 2, 2, complex(0.0011680804588461839, 0.006101515081950891),
                    complex(144.79464254274615, 533.9355326013214)),
                (3, 1, 1, complex(12.51689296620168, 5.02318787979271),
                    complex(996084.9289031818, 399677.1758765848)),
                (3, 2, 1, complex(1.455776974837841e-08, 8.193937638111909e-09),
                    complex(0.0011585067796976401, 0.00065198838864868)),
                (3, 1, 2, complex(1.455776974837841e-08, 8.193937638111909e-09),
                    complex(0.0011585067796976401, 0.00065198838864868)),
                (3, 2, 2, complex(12.51689296620168, 5.02318787979271),
                    complex(996084.9289031818, 399677.1758765848))]),
        (name = "two buried insulated wires", radii = [0.0435, 0.0435],
            placements = [Placement(0.0, -1.0), Placement(1.0, -1.0)],
            earth = homogeneous(rho = 0.1),
            state = [
                (1, 1, complex(8.643091791240205e-05, 8.643091791245012e-05),
                    complex(0.9999135728171318, 3.734828966945854e-09),
                    complex(0.05944678698755287, -1.1098831454672168e-10)),
                (1, 2, complex(8.643091791240205e-05, 8.643091791245012e-05),
                    complex(0.9999135728171318, 3.734828966945854e-09),
                    complex(0.05944678698755287, -1.1098831454672168e-10)),
                (2, 1, complex(0.008643091767200767, 0.00864309181528445),
                    complex(0.9913941520257152, 3.7030076383130116e-05),
                    complex(0.05944678695991977, -1.1098831447549736e-06)),
                (2, 2, complex(0.008643091767200767, 0.00864309181528445),
                    complex(0.9913941520257152, 3.7030076383130116e-05),
                    complex(0.05944678695991977, -1.1098831447549736e-06)),
                (3, 1, complex(0.8642851376172453, 0.8643332213000304),
                    complex(0.40666213121824585, 0.1567702794234474),
                    complex(0.05685507777132759, -0.01043369383465527)),
                (3, 2, complex(0.8642851376172453, 0.8643332213000304),
                    complex(0.40666213121824585, 0.1567702794234474),
                    complex(0.05685507777132759, -0.01043369383465527))],
            terms = [
                (1, 1, 1, complex(9.903351820293542e-08, 1.209350897524307e-06),
                    complex(0.015708823140074858, 0.14424157375524968)),
                (1, 2, 1, complex(9.903234809067121e-08, 8.153619364863723e-07),
                    complex(0.01570857199574813, 0.11177313743995486)),
                (1, 1, 2, complex(9.903234809067121e-08, 8.153619364863723e-07),
                    complex(0.01570857199574813, 0.11177313743995486)),
                (1, 2, 2, complex(9.903351820293542e-08, 1.209350897524307e-06),
                    complex(0.015708823140074858, 0.14424157375524968)),
                (2, 1, 1, complex(0.001168301159428558, 0.006073250378889802),
                    complex(144.82242471244095, 531.6960031463639)),
                (2, 2, 1, complex(0.0011129984037097581, 0.002113298674516299),
                    complex(137.7562996441238, 207.15264847001256)),
                (2, 1, 2, complex(0.0011129984037097581, 0.002113298674516299),
                    complex(137.7562996441238, 207.15264847001256)),
                (2, 2, 2, complex(0.001168301159428558, 0.006073250378889802),
                    complex(144.82242471244095, 531.6960031463639)),
                (3, 1, 1, complex(12.488532839756836, 4.715192317939723),
                    complex(993826.7382296916, 375167.79344574944)),
                (3, 2, 1, complex(1.4754996746279764e-08, 8.623953115574929e-09),
                    complex(0.0011742035091816137, 0.0006862070598237097)),
                (3, 1, 2, complex(1.4754996746279764e-08, 8.623953115574929e-09),
                    complex(0.0011742035091816137, 0.0006862070598237097)),
                (3, 2, 2, complex(12.488532839756836, 4.715192317939723),
                    complex(993826.7382296916, 375167.79344574944))]),
        (name = "a buried coaxial cable", radii = [0.0325],
            placements = [Placement(0.0, -1.0)],
            earth = homogeneous(rho = 100, eps_r = 10),
            state = [
                (1, 1, complex(2.04203521915319e-06, 2.0420352305135415e-06),
                    complex(0.9999979579668657, 2.0849496721851195e-12),
                    complex(3.31830724035422e-05, 1.8457114416005542e-13)),
                (2, 1, complex(0.00020419784238582406, 0.00020420920273885027),
                    complex(0.9997958230034147, 2.084528230099945e-08),
                    complex(3.318307240358068e-05, 1.8457114416016235e-09)),
                (3, 1, complex(0.015658665950048743, 0.026630032678149534),
                    complex(0.984349097084341, 0.00020524415729066272),
                    complex(3.318692131361649e-05, 1.845818432835047e-05))],
            terms = [
                (1, 1, 1, complex(9.870676853624838e-08, 1.6802324557995022e-06),
                    complex(15.707996048229546, 216.22224805981273)),
                (2, 1, 1, complex(0.000997426549997924, 0.011007121428432628),
                    complex(157156.91396419608, 1241433.1325170435)),
                (3, 1, 1, complex(14.818298565790123, 45.18913105671043),
                    complex(2435070348.5060434, 2088460127.983585))]),
        (name = "an overhead and a buried bare wire", radii = [0.0425, 0.0425],
            placements = [Placement(0.0, 10.0), Placement(1.0, -1.0)],
            earth = homogeneous(rho = 0.1),
            state = [
                (1, 1, complex(0.0, 8.907341340869984e-11),
                    complex(1.0, 0.0),
                    complex(0.0, 3.156867332565835e-14)),
                (1, 2, complex(8.444400025924338e-05, 8.444400025929034e-05),
                    complex(0.999915559565035, 3.5650935264235386e-09),
                    complex(0.05674501730546565, -1.0112762017883535e-10)),
                (2, 1, complex(0.0, 8.907341340869983e-07),
                    complex(0.9999999999998016, 0.0),
                    complex(0.0, 3.1568673325661486e-10)),
                (2, 2, complex(0.00844440000243753, 0.00844440004941584),
                    complex(0.9915911534807212, 3.535413735005518e-05),
                    complex(0.05674501728143213, -1.0112762011926023e-06)),
                (3, 1, complex(0.0, 0.008907341340869984),
                    complex(0.9999801649159178, 0.0),
                    complex(1.93304290827255e-22, 3.1568986414995888e-06)),
                (3, 2, complex(0.8444165137639753, 0.8444634920747424),
                    complex(0.4161523516456092, 0.15270107075218897),
                    complex(0.05447840661298058, -0.009557481831077413))],
            terms = [
                (1, 1, 1, complex(9.557821194533482e-08, 1.2157922597268011e-06),
                    complex(110618490449.433, -0.01993782365860378)),
                (1, 2, 1, complex(9.724616285403196e-08, 5.152690967698774e-07),
                    complex(0.015328831449057806, 0.07216474056647124)),
                (1, 1, 2, complex(9.725437689305671e-08, 5.153126096440578e-07),
                    complex(0.005396651426301309, -0.047179316525339444)),
                (1, 2, 2, complex(9.903332164789415e-08, 1.2122712849237515e-06),
                    complex(0.015708791953322424, 0.1444738754182352)),
                (2, 1, 1, complex(0.00024718208995873917, 0.008040575171935635),
                    complex(110618500094.26599, -5136.1552291958205)),
                (2, 2, 1, complex(0.0004014036917291605, 0.0003920139961180333),
                    complex(42.477367973575646, 34.14843650912342)),
                (2, 1, 2, complex(0.00040482175204522414, 0.000395323902011549),
                    complex(1742.8630998202075, -4120.284772104932)),
                (2, 2, 2, complex(0.0011680804588461839, 0.006101515081950891),
                    complex(144.79464254274615, 533.9355326013214)),
                (3, 1, 1, complex(27.089301992896484, 58.94174617873044),
                    complex(84348236879.01994, -38572301620.65046)),
                (3, 2, 1, complex(7.91314017905357e-11, -3.940728440021105e-10),
                    complex(6.231389788666503e-06, -3.145106422268027e-05)),
                (3, 1, 2, complex(-1.386464675648153e-10, -8.960505978313245e-10),
                    complex(-3.1904070725291853, 0.7285292482681873)),
                (3, 2, 2, complex(12.51689296620168, 5.02318787979271),
                    complex(996084.9289031818, 399677.1758765848))])
    ]
    @testset "the state and the pair terms of $(case.name)" for case in cases
        designs = [build(CableDesign, "wire",
                       Terminal(:core, Region(:core_metal, Circle(r), copper)))
                   for r in case.radii]
        system = build(CableSystem, designs, case.placements;
            connections = [(core = i,) for i in eachindex(designs)], id = "unified",
            line_length = 1.0)
        problem = LineParametersProblem(system, case.earth; frequencies = [
            0.1, 1.0e3, 1.0e7])
        @test Engine.validate(problem, formula) === problem
        input = Engine.CoaxialInput(problem, reductions)
        buffers = Engine.initialize_buffers(formula, problem, input,
            (gamma = zeros(ComplexF64, 3),))
        for (k, f) in enumerate(problem.frequencies)
            media = map(layer -> Engine.properties(law, layer, f), problem.earth.layers)
            Engine.prepare!(formula, buffers, input, media, im * (2 * π * f), k)
            for (index, p, x, A, F) in case.state
                index == k || continue
                @test isapprox(buffers.unified.x[p], x; rtol = 1e-12)
                @test isapprox(buffers.unified.A[p], A; rtol = 1e-12)
                @test isapprox(buffers.unified.F[p], F; rtol = 1e-12)
            end
            for (index, p, q, k_target, h_target) in case.terms
                index == k || continue
                @test isapprox(buffers.unified.K[p, q], k_target; rtol = 1e-12)
                @test isapprox(buffers.unified.H[p, q], h_target; rtol = 1e-12)
            end
        end
    end

    @testset "Air receiver at the Unified split" begin
        r = 0.0425
        wire = build(CableDesign, "wire",
            Terminal(:core, Region(:core_metal, Circle(r), copper)))
        for placements in ([Placement(0.0, 1.0), Placement(1.0, 1.0)],
            [Placement(0.0, 1.0), Placement(1.0, -1.0)])
            system = build(CableSystem, [wire, wire], placements;
                connections = [(core = 1,), (core = 2,)], id = "unified",
                line_length = 1.0)
            problem = LineParametersProblem(system, homogeneous(rho = 100, eps_r = 10);
                frequencies = [50.0, 1e4, 1e6])
            input = Engine.CoaxialInput(problem, reductions)
            buffers = [Engine.initialize_buffers(formula, problem, input,
                           (gamma = fill(im * scale / r, length(problem.frequencies)),))
                       for scale in (299.0, 299.9, 300.1, 301.0)]
            terms = [(H = zeros(ComplexF64, 2, 2, 3), K = zeros(ComplexF64, 2, 2, 3),
                         Pg = zeros(ComplexF64, 2, 2, 3)) for _ in buffers]
            Pg = zeros(ComplexF64, 2, 2)
            for (k, f) in enumerate(problem.frequencies)
                media = map(layer -> Engine.properties(law, layer, f), problem.earth.layers)
                for (buffer, term) in zip(buffers, terms)
                    Engine.prepare!(formula, buffer, input, media, im * (2 * π * f), k)
                    Engine.potential!(
                        formula, Pg, buffer, input, media, im * (2 * π * f), k)
                    term.H[:, :, k] = buffer.unified.H
                    term.K[:, :, k] = buffer.unified.K
                    term.Pg[:, :, k] = Pg
                end
                @test abs(real(buffers[1].unified.x[1])) < 300
                @test abs(real(buffers[2].unified.x[1])) < 300
                @test abs(real(buffers[3].unified.x[1])) > 300
                @test abs(real(buffers[4].unified.x[1])) > 300
            end
            for quantity in (:H, :K, :Pg)
                A₁ = getfield(terms[1], quantity)
                A₂ = getfield(terms[2], quantity)
                A₃ = getfield(terms[3], quantity)
                A₄ = getfield(terms[4], quantity)
                below = A₂ + (A₂ - A₁) * (0.1 / 0.9)
                above = A₃ - (A₄ - A₃) * (0.1 / 0.9)
                @test maximum(maximum(abs, above[:, :, k] - below[:, :, k]) /
                              maximum(abs, below[:, :, k])
                for k in eachindex(problem.frequencies)) <= 1e-7
            end
        end
    end

    # The earth-return matrices of the references at their 101 frequencies, against the stage
    # targets. The coaxial references share the cable, the earth and these matrices.
    frequencies = 10.0 .^ range(-1, 7; length = 101)
    zero_gamma = zeros(ComplexF64, length(frequencies))
    stages = [
        ("two_bare_wires", "stages.csv", 0.0425,
            [Placement(0.0, -1.0), Placement(1.0, -1.0)], homogeneous(rho = 0.1), zero_gamma),
        ("two_insulated_wires", "stages.csv", 0.0435,
            [Placement(0.0, -1.0), Placement(1.0, -1.0)], homogeneous(rho = 0.1), zero_gamma),
        ("coaxial_circular_strands", "stages.csv", 0.0325, [Placement(0.0, -1.0)],
            homogeneous(rho = 100, eps_r = 10), zero_gamma),
        ("overhead_bare_wires", "stages.csv", 0.0425,
            [Placement(0.0, 1.0), Placement(1.0, 1.0)], homogeneous(rho = 0.1), zero_gamma),
        ("mixed_bare_wires", "stages.csv", 0.0425,
            [Placement(0.0, 1.0), Placement(1.0, -1.0)], homogeneous(rho = 0.1), zero_gamma),
        ("two_insulated_wires", "stages-gamma.csv", 0.0435,
            [Placement(0.0, -1.0), Placement(1.0, -1.0)], homogeneous(rho = 0.1),
            [im * 1.01 * 2π * f / 299792458 for f in frequencies])]
    @testset "the earth-return matrices of $case, $file" for (
        case, file, radius, placements, earth, gamma) in stages

        designs = [build(CableDesign, "wire",
                       Terminal(:core, Region(:core_metal, Circle(radius), copper)))
                   for _ in placements]
        system = build(CableSystem, designs, placements;
            connections = [(core = i,) for i in eachindex(designs)], id = case,
            line_length = 1.0)
        problem = LineParametersProblem(system, earth; frequencies)
        input = Engine.CoaxialInput(problem, reductions)
        buffers = Engine.initialize_buffers(formula, problem, input, (gamma = gamma,))
        c = length(placements)
        targets = Dict(matrix => zeros(ComplexF64, c, c, length(frequencies))
        for matrix in ("Zg", "Pg"))
        for line in eachline(joinpath(@__DIR__, "..", "references", case, file))
            fields = split(line, ',')
            startswith(line, "analytical") && haskey(targets, fields[3]) || continue
            k = searchsortedfirst(frequencies, parse(Float64, fields[2]))
            targets[fields[3]][parse(Int, fields[4]), parse(
                Int, fields[5]), k] = complex(
                parse(Float64, fields[6]), parse(Float64, fields[7]))
        end
        Zg = zeros(ComplexF64, c, c)
        Pg = zeros(ComplexF64, c, c)
        for (k, f) in enumerate(frequencies)
            media = map(layer -> Engine.properties(law, layer, f), problem.earth.layers)
            Engine.prepare!(formula, buffers, input, media, im * (2 * π * f), k)
            Engine.impedance!(formula, Zg, buffers, input, media, im * (2 * π * f), k)
            Engine.potential!(formula, Pg, buffers, input, media, im * (2 * π * f), k)
            for (matrix, target) in
                ((Zg, targets["Zg"][:, :, k]), (Pg, targets["Pg"][:, :, k]))
                @test maximum(abs, matrix - target) <= 1e-9 * maximum(abs, diag(target))
            end
        end
    end

    @test isapprox(EarthReturn.i0m1(complex(0.6, 0.8)), besseli(0, complex(0.6, 0.8)) - 1;
        rtol = 1e-14)
    @test isapprox(EarthReturn.i0m1(1.0e-3), 1.0e-6 / 4 + 1.0e-12 / 64 + 1.0e-18 / 2304;
        rtol = 1e-15)

    # A cable that touches the interface leaves its circumference in two media, and a Γ
    # vector has one value per frequency.
    wire = build(CableDesign, "wire",
        Terminal(:core, Region(:core_metal, Circle(0.0425), copper)))
    touching = LineParametersProblem(
        build(CableSystem, [wire], [Placement(0.0, -0.0425)]; connections = [(core = 1,)],
            id = "touching", line_length = 1.0),
        homogeneous(rho = 0.1); frequencies = [1.0e3])
    @test_throws DomainError Engine.validate(touching, formula)
    buried = LineParametersProblem(
        build(CableSystem, [wire, wire], [Placement(0.0, -1.0), Placement(0.085, -1.0)];
            connections = [(core = 1,), (core = 2,)], id = "pair", line_length = 1.0),
        homogeneous(rho = 0.1); frequencies = [1.0e3, 1.0e4, 1.0e5])
    @test_throws DimensionMismatch Engine.validate(buried,
        CoaxialFormulation(; Γ = [0.1im, 0.2im]))
    # The formula computes cables in contact, with no guard between them. After the first call,
    # the preparation allocates only the pivot vector of its LU decomposition.
    contact = LineCableModels.Engine.EarthReturn.Formula(Val(:unified))
    input = LineCableModels.Engine.CoaxialInput(buried, reductions)
    buffers = LineCableModels.initialize_buffers(contact, buried, input,
        (gamma = zeros(ComplexF64, 3),))
    media = map(layer -> LineCableModels.properties(law, layer, 1.0e3), buried.earth.layers)
    s = im * (2 * π * 1.0e3)
    LineCableModels.prepare!(contact, buffers, input, media, s, 1)
    @test all(isfinite, buffers.unified.Ze) && all(isfinite, buffers.unified.Pe)
    block = Matrix{ComplexF64}(I, 2, 2)
    lu!(copy(block))
    pivots = @allocated lu!(block)
    # The arguments pass through an anonymous function, where their types are concrete.
    @test only(map(
        (f, b, c, m, z) -> @allocated(LineCableModels.prepare!(f, b, c, m, z, 1)),
        (contact,), (buffers,), (input,), (media,), (s,))) == pivots
end
