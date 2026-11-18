using LineCableModels.BaseGeometry: Contour, gmd
using LineCableModels.DataModel: ArcPolygon, gmr

@testset "Shape geometric mean distance" begin
    # Float64 bounds measured with gmd at τ = 1e-6 on Julia 1.12.7, x86_64 Linux.
    disk = Circle(3e-3)
    for (A, B) in ((disk, disk),
        (Annulus(1e-3, 3e-3), Annulus(1e-3, 3e-3)),
        (Annulus(2.7e-3, 3e-3), Annulus(2.7e-3, 3e-3)),
        (Annulus(4e-3, 5e-3), disk),
        (Annulus(4e-3, 5e-3), Rectangle(3e-3, 2e-3, Placement(0.5e-3, 0.3e-3))),
        (Annulus(4e-3, 5e-3), ArcPolygon(Contour(Rectangle(3e-3, 2e-3)).outer)),
        (Annulus(4e-3, 5e-3), Rectangle(3e-3, 2e-3, Placement(8e-3, 0))),
        (disk, Rectangle(3e-3, 2e-3, Placement(8e-3, 0))),
        (disk, Circle(2e-3, Placement(8e-3, 1e-3))))
        @test abs(log(gmd(A, B)) - log(gmd(Contour(A), Contour(B)))) <= 1e-9
        @test gmd(A, B) ≈ gmd(B, A) rtol = 1e-12
    end
    @test isempty(Test.detect_ambiguities(LineCableModels.BaseGeometry,
        LineCableModels.DataModel; recursive = true))

    @testset "Thin annulus own term" begin
        # On Julia 1.12.7, x86_64 Linux, the measured log g error was at most 0.583 ulp.
        # The 1 ulp bound leaves a margin of 0.417 ulp.
        # Compare log g and its inner-radius derivative with the closed form at 256 bits.
        b = 0.03
        for ε in (1e-12, 1e-8, 1e-4, 0.3, 0.5)
            a = b * (1 - ε)
            annulus = Annulus(a, b)
            expected = setprecision(BigFloat, 256) do
                let a = BigFloat(a), b = BigFloat(b)
                    log(b) - (a^4 * log(b / a) / (b^2 - a^2)^2 -
                     (3a^2 - b^2) / (4(b^2 - a^2)))
                end
            end
            @test abs(BigFloat(log(gmd(annulus, annulus))) - expected) <=
                  eps(abs(Float64(expected)))
        end
        # Julia 1.12.7 on x86_64 Linux with ForwardDiff 1.4.8 gives the errors below the 6 ulp bound.
        # The inner-radius partial errors were 1.05 ulps at ε = 1e-6 and 0.154 ulp at ε = 0.5.
        for ε in (1e-6, 0.5)
            a = b * (1 - ε)
            annulus = Annulus(ForwardDiff.Dual(a, 1.0), b)
            expected = setprecision(BigFloat, 256) do
                let b = BigFloat(b)
                    ForwardDiff.derivative(BigFloat(a)) do a
                        log(b) - (a^4 * log(b / a) / (b^2 - a^2)^2 -
                         (3a^2 - b^2) / (4(b^2 - a^2)))
                    end
                end
            end
            @test abs(BigFloat(only(ForwardDiff.partials(log(gmd(annulus, annulus))))) -
                      expected) <= 6eps(abs(Float64(expected)))
        end
        @testset "Two further terms preserve the inner-radius partial" begin
            # The largest improvement in the 400-width sweep occurs at its upper endpoint.
            # On Julia 1.12.7, x86_64 Linux, with ForwardDiff 1.4.8, the error was 0.154149 ulp.
            # Without the two further terms, it was 5.845851 ulps. The bound is 4 ulps.
            ε = 0.5
            a = b * (1 - ε)
            annulus = Annulus(ForwardDiff.Dual(a, 1.0), b)
            expected = setprecision(BigFloat, 256) do
                let b = BigFloat(b)
                    ForwardDiff.derivative(BigFloat(a)) do a
                        log(b) - (a^4 * log(b / a) / (b^2 - a^2)^2 -
                         (3a^2 - b^2) / (4(b^2 - a^2)))
                    end
                end
            end
            @test abs(BigFloat(only(ForwardDiff.partials(log(gmd(annulus, annulus))))) -
                      expected) <= 4eps(abs(Float64(expected)))
        end
    end

    copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
    polymer = Material(:insulator, 1e14, 2.5)
    bar = Region(:core_bar, Rectangle(10e-3, 5e-3), copper)
    plated = build(CableDesign, "plated",
        Terminal(:core, Stack(bar, Region(:core_plating, Offset(0.5e-3), copper))),
        Region(:core_insulation, Offset(1e-3), polymer))
    bare = build(CableDesign, "bare", Terminal(:core, bar),
        Region(:core_insulation, Offset(1.5e-3), polymer))
    @test only(coaxial_input(plated).material).mu_r ≈ 2.050 rtol = 1e-3
    @test only(coaxial_input(bare).material).mu_r ≈ 2.423 rtol = 1e-3
    shapes = [plated.regions[i].shape for i in only(plated.terminals).regions]
    shares = area.(shapes) ./ sum(area, shapes)
    weighted = sum(shares[i] * shares[j] * log(gmd(shapes[i], shapes[j]))
    for i in eachindex(shapes), j in eachindex(shapes))
    outline = Contour(last(shapes))
    outer = Contour(outline.outer, empty(outline.holes))
    @test abs(weighted - log(gmd(outer, outer))) <= 1e-9
end

@testset "The geometric mean radius of a conductor reads its shape" begin
    a, b = 1e-3, 3e-3
    # Maxwell's internal term of the annulus from a to b.
    internal = a^4 / (b^2 - a^2)^2 * log(b / a) - (3a^2 - b^2) / (4(b^2 - a^2))
    bar = Rectangle(10e-3, 5e-3)
    away = Placement(8e-3, 1e-3)
    for μ in (1.0, 36.6326)
        @test log(gmr(Circle(b), μ)) ≈ log(b) - μ / 4
        @test log(gmr(Annulus(a, b), μ)) ≈ log(b) - μ * internal
        @test log(gmr(bar, μ)) ≈ log(gmd(bar, bar)) - (μ - 1) / 4
        # A copy placed away from the origin keeps the radius.
        @test gmr(Circle(b, away), μ) == gmr(Circle(b), μ)
        @test gmr(Annulus(a, b, away), μ) == gmr(Annulus(a, b), μ)
    end
    @test gmr(Annulus(a, b), 1.0) ≈ gmd(Annulus(a, b), Annulus(a, b))
    # Each copy of a strand is a rigid image of the first copy of its course, with the same
    # radius. On Julia 1.12.7, x86_64 Linux, the bent rectangular strands differ by 1.3e-14
    # at most, relative.
    design = only(Dict(REFERENCES)["coaxial_rectangular_strands"].system.designs)
    for row in design.repeats
        copies = [design.regions[only(copy)] for copy in row.copies]
        g = [gmr(placed.shape, placed.source.material.mu_r) for placed in copies]
        @test g ≈ fill(first(g), length(g)) rtol = 1e-12
    end
end

@testset "The tube of a terminal reproduces its regions in parallel" begin
    # Maxwell's internal term of the tube from a to b, which is 1/4 for a solid one.
    internal(a, b) = (iszero(a) ? 0.0 : a^4 / (b^2 - a^2)^2 * log(b / a)) -
                     (3a^2 - b^2) / (4(b^2 - a^2))
    # The lay ratios that the reference designs declare, by the names of the laid regions.
    q = Dict(:core_strands_1 => 12, :core_strands_2 => 12, :core_strands_3 => 12,
        :armour_wires => 10, :screen => 5.5012442232492)
    for (case, problem) in REFERENCES, design in unique(problem.system.designs)

        input = coaxial_input(design)
        for (k, terminal) in enumerate(design.terminals)
            tube = input.material[k]
            axis = Placement(
                input.positions[findfirst(span -> k in span, input.assemblies)]...)
            regions = design.regions[terminal.regions]
            shapes = [LineCableModels.resolve(inv(axis), region.shape)
                      for region in regions]
            μ = [region.source.material.mu_r for region in regions]
            # A region laid at the ratio q is longer than the cable by √(1 + (π/q)²).
            ℓ = [sqrt(1 + (π / get(q, region.source.tag, Inf))^2) for region in regions]
            G = [area(shapes[i]) / (regions[i].source.material.rho * ℓ[i])
                 for i in eachindex(regions)]
            w = G ./ sum(G)
            # The own term of a circle or an annulus is Maxwell's, in closed form. A region
            # that is not round takes its geometric mean distance from itself, with 1/4 inside.
            own = [S isa Union{Circle, Annulus} ?
                   log(r_ex(S)) - μ_i * internal(r_in(S), r_ex(S)) :
                   log(gmd(S, S)) - (μ_i - 1) / 4 for (S, μ_i) in zip(shapes, μ)]
            log_g = sum(w .^ 2 .* own) +
                    2sum(
                (w[i] * w[j] * log(gmd(shapes[i], shapes[j]))
                for i in eachindex(shapes) for j in (i + 1):length(shapes));
                init = 0.0)
            r_i, r_e = input.r_in[k], input.r_ex[k]
            # Both sides evaluate `gmd` on the same centered shapes of the regions. Each copy of a strand has its own term here, where the tube takes one per course. The summation order and the own terms differ, with roundoff far below the default `≈` tolerance. A chain fault moves either quantity by 1e-3 or more.
            @test tube.rho / (π * (r_e^2 - r_i^2)) ≈ inv(sum(G))
            @test log(r_e) - tube.mu_r * internal(r_i, r_e) ≈ log_g
        end
    end
end
