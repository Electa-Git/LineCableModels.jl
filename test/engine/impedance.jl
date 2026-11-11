# The targets come from the schelkunoff1934 and ametani1980 formulas of commit 989af96c of
# this repository, with s = j 2π f. They are the equivalent tubes
# and the dielectric layers of the coaxial references.
using SpecialFunctions: besseli, besselk, besselkx

@testset "Engine impedances" begin
    Engine = LineCableModels.Engine
    schelkunoff = Engine.InternalImpedance.Formula(Val(:schelkunoff1934))
    ametani = Engine.InsulationImpedance.Formula(Val(:ametani1980))
    # The radii in m, the resistivity in Ω·m and the relative permeability of each tube, and the
    # radii and the relative permeability of each layer.
    tubes = Dict("core" => (0.0, 10.5e-3, 2.3580693965940444e-8, 1.0588802203492307),
        "sheath" => (20.5e-3, 21.5e-3, 1.7241e-8, 1.0),
        "armour" => (24.5e-3, 28.5e-3, 2.4270115498150639e-7, 5.8715848278702119))
    layers = Dict("insulation" => (10.5e-3, 20.5e-3, 1.416791518882097),
        "bedding" => (21.5e-3, 24.5e-3, 1.0),
        "jacket" => (28.5e-3, 32.5e-3, 1.1305521523424971))
    # Each row gives the tube, the frequency in Hz and the three impedances in Ω/m.
    surfaces = [
        ("core", 0.1, complex(0.0, 0.0),
            complex(6.808134793446447e-05, 3.3265701889145045e-08),
            complex(0.0, 0.0)),
        ("sheath", 0.1, complex(0.00013066620828468208, 2.0428433857038906e-09),
            complex(0.00013066620828439208, 1.947837934674311e-09),
            complex(0.00013066620827311955, -9.972179159207138e-10)),
        ("armour", 0.1, complex(0.00036440649618594623, 4.005977752554661e-08),
            complex(0.0003644064960718631, 3.444329799040416e-08),
            complex(0.0003644064947117064, -1.8540977975289263e-08)),
        ("core", 1000.0, complex(0.0, 0.0),
            complex(0.00016894965120533065, 0.000148661054855278),
            complex(0.0, 0.0)),
        ("sheath", 1000.0, complex(0.00013128849692723698, 2.0401301171856737e-05),
            complex(0.00013125955784695488, 1.9452508437422897e-05),
            complex(0.0001301346871555814, -9.946514884520603e-06)),
        ("armour", 1000.0, complex(0.0004392255487740415, 0.0003789784108287536),
            complex(0.00042874041454241986, 0.0003258431022972495),
            complex(0.00030444609800333043, -0.00016604709179759828)),
        ("core", 10000000.0, complex(0.0, 0.0),
            complex(0.015066199980529618, 0.01504915073770315),
            complex(0.0, 0.0)),
        ("sheath", 10000000.0, complex(0.006401857925416989, 0.006405120142563304),
            complex(0.006110177928862835, 0.006107207683649767),
            complex(-2.9188975931127455e-23, -1.6789477584050798e-24)),
        ("armour", 10000000.0, complex(0.048692353093421445, 0.04872449707673149),
            complex(0.041909772808230696, 0.04188597471717819),
            complex(-2.503604924833728e-55, 8.031441151984553e-56))
    ]
    # Each row gives the layer, the frequency in Hz and the magnetic term in Ω/m.
    terms = [
        ("insulation", 0.1, complex(0.0, 1.1911710960859268e-07)),
        ("bedding", 0.1, complex(0.0, 1.6414216219680304e-08)),
        ("jacket", 0.1, complex(0.0, 1.865882352450562e-08)),
        ("insulation", 1000.0, complex(0.0, 0.0011911710960859265)),
        ("bedding", 1000.0, complex(0.0, 0.00016414216219680301)),
        ("jacket", 1000.0, complex(0.0, 0.00018658823524505614)),
        ("insulation", 10000000.0, complex(0.0, 11.911710960859267)),
        ("bedding", 10000000.0, complex(0.0, 1.6414216219680304)),
        ("jacket", 10000000.0, complex(0.0, 1.8658823524505617))
    ]

    @testset "the surfaces of the tubes, and of the solid core" begin
        for (part, f, inner, outer, transfer) in surfaces
            surface = schelkunoff(tubes[part]..., im * (2π * f))
            @test isapprox(surface.inner, inner; rtol = 1.0e-13)
            @test isapprox(surface.outer, outer; rtol = 1.0e-13)
            @test isapprox(surface.transfer, transfer; rtol = 1.0e-13)
        end
    end

    @testset "the magnetic term across the layers" begin
        for (part, f, term) in terms
            @test isapprox(ametani(layers[part]..., im * (2π * f)), term;
                rtol = 1.0e-13)
        end
        @test iszero(ametani(10.5e-3, 10.5e-3, 1.0, im * (2π * 50.0)))
    end

    @testset "a tube has ordered radii, a conductive material and a frequency" begin
        s = im * (2π * 50.0)
        @test_throws DomainError schelkunoff(2.0e-3, 1.0e-3, 1.7241e-8, 1.0, s)
        @test_throws DomainError schelkunoff(0.0, 1.0e-3, 0.0, 1.0, s)
        @test_throws DomainError schelkunoff(0.0, 1.0e-3, 1.7241e-8, -1.0, s)
        @test_throws DomainError schelkunoff(0.0, 1.0e-3, 1.7241e-8, 1.0, 0.0im)
    end

    @testset "the impedances are inferred" begin
        s = im * (2π * 50.0)
        @test @inferred(schelkunoff(tubes["sheath"]..., s)).outer isa ComplexF64
        @test @inferred(ametani(layers["insulation"]..., s)) isa ComplexF64
    end

    # With SpecialFunctions 2.9.1 on Julia 1.12.7, the bounds below hold against the surfaces at
    # an order in Arb at 256 bits. Each wall gives its radii in m, its resistivity
    # in Ω·m and its relative permeability. The walls are a steel pipe of μ_r 10, the shield of
    # Table I of Kane 1995 and a steel pipe of μ_r 300. A lead sheath, an aluminum foil and a
    # thick copper tube follow. Each solid gives its radius, its resistivity and its
    # permeability.
    μ₀ = Engine.μ₀
    walls = ((62.0e-3, 65.0e-3, 7.0e-8, 10.0), (26.8e-3, 28.8e-3, 1 / 4.76e6, 1.0),
        (62.0e-3, 65.0e-3, 1.7e-7, 300.0), (30.0e-3, 31.0e-3, 2.14e-7, 1.0),
        (40.0e-3, 40.2e-3, 2.8e-8, 1.0), (20.0e-3, 60.0e-3, 1.7e-8, 1.0))
    solids = ((19.5e-3, 3.365e-8, 1.0), (7.0e-3, 1 / 3.33e7, 1.0), (10.0e-3, 1.7e-7, 300.0))

    @testset "the surfaces at order 0 are those of a net current" begin
        for part in ("core", "sheath"), f in (50.0, 1.0e6)

            s = im * (2π * f)
            @test schelkunoff(0, tubes[part]..., s) == schelkunoff(tubes[part]..., s)
        end
    end

    @testset "the surfaces of a wall that conducts nothing, at an order" begin
        # At 1e-6 Hz the field of the order n in the metal is the sum of the powers n and -n of
        # the radius. Measured: 5.6e-5 for a tube and 3.5e-7 for a solid.
        s = im * (2π * 1.0e-6)
        for (r_in, r_ex, ρ, μ_r) in walls, n in (1, 2, 3, 5, 10, 20, 40)

            surface = schelkunoff(n, r_in, r_ex, ρ, μ_r, s)
            ϱ = (r_in / r_ex)^n
            z = s * μ₀ * μ_r / (2π * n)
            @test isapprox(surface.inner, z * (1 + ϱ^2) / (1 - ϱ^2); rtol = 1.0e-4)
            @test isapprox(surface.outer, z * (1 + ϱ^2) / (1 - ϱ^2); rtol = 1.0e-4)
            @test isapprox(surface.transfer, z * 2ϱ / (1 - ϱ^2); rtol = 1.0e-4)
        end
        for (r_ex, ρ, μ_r) in solids, n in (1, 2, 3, 5, 10, 20, 40)

            @test isapprox(schelkunoff(n, 0.0, r_ex, ρ, μ_r, s).outer,
                s * μ₀ * μ_r / (2π * n); rtol = 1.0e-6)
        end
    end

    @testset "the inner surface of a thick tube is that of a wall without end" begin
        # A copper tube of 20 to 600 mm, at 1 kHz to 1 MHz. Measured: 5.2e-16.
        for f in (1.0e3, 1.0e5, 1.0e6), n in (1, 2, 3, 5, 10)

            s = im * (2π * f)
            w = 20.0e-3 * sqrt(s * μ₀ / 1.7e-8)
            @test isapprox(schelkunoff(n, 20.0e-3, 600.0e-3, 1.7e-8, 1.0, s).inner,
                s * μ₀ / 2π * besselkx(n, w) /
                (w * besselkx(n - 1, w) + n * besselkx(n, w));
                rtol = 1.0e-13)
        end
    end

    @testset "the terms of Kane, Ahmad and Auriol from the surfaces at an order" begin
        # The shield of their Table I, 26.8 to 28.8 mm, and its 7 mm core. The order n of their
        # (4) and (14), from their (5), (6) and Δ, is the wall's answer Z in parallel with the
        # medium z. The term n of their (10) is the core's outer surface in parallel with z.
        # Measured: 5.7e-15 and 2.7e-15.
        c₁, c₂, σ, μ_r = 26.8e-3, 28.8e-3, 4.76e6, 1.0
        b_i, b_j, θ = 12.0e-3, 13.0e-3, 2π / 3
        a, d, σ_k = 7.0e-3, 21.65e-3, 3.33e7
        for f in (1.0, 1.0e2, 1.0e3, 1.0e4, 1.0e5), n in (1, 2, 3, 5, 8)

            s = im * (2π * f)
            z = s * μ₀ / (2π * n)
            m = sqrt(s * μ₀ * μ_r * σ)
            u, v = c₁ * m, c₂ * m
            dI_u = besseli(n - 1, u) - n / u * besseli(n, u)
            dI_v = besseli(n - 1, v) - n / v * besseli(n, v)
            dK_u = -besselk(n - 1, u) - n / u * besselk(n, u)
            dK_v = -besselk(n - 1, v) - n / v * besselk(n, v)
            Δ = (μ_r * n / u * besseli(n, u) - dI_u) *
                (μ_r * n / v * besselk(n, v) + dK_v) -
                (μ_r * n / v * besseli(n, v) + dI_v) * (μ_r * n / u * besselk(n, u) - dK_u)
            P_n = Q_n = (b_i / c₁)^n / c₁
            A_n = m / (2π * Δ) * (P_n + Q_n) * (μ_r * n / v * besselk(n, v) + dK_v)
            B_n = -m / (2π * Δ) * (P_n + Q_n) * (μ_r * n / v * besseli(n, v) + dI_v)
            printed = (b_j / c₁)^n * cos(n * θ) *
                      (A_n * besseli(n, u) + B_n * besselk(n, u)) / σ
            wall = schelkunoff(n, c₁, c₂, 1 / σ, μ_r, s)
            Z = wall.inner - wall.transfer^2 / (wall.outer + z)
            @test isapprox(
                2 * (b_i * b_j / c₁^2)^n * cos(n * θ) * Z * z / (Z + z), printed;
                rtol = 1.0e-12)
            m_k = sqrt(s * μ₀ * σ_k)
            dI_k = besseli(n - 1, a * m_k) - n / (a * m_k) * besseli(n, a * m_k)
            printed = (s * μ₀ * σ_k * a^(n - 1) / π) * (a / d)^n * (1 / d)^n *
                      besseli(n, a * m_k) / (n / a * besseli(n, a * m_k) + m_k * dI_k) / σ_k
            core = schelkunoff(n, 0.0, a, 1 / σ_k, 1.0, s).outer
            @test isapprox(2 * (a / d)^(2n) * core * z / (core + z), printed; rtol = 1.0e-12)
        end
    end

    @testset "the surfaces at an order against their values in Arb" begin
        # Each row gives the radii, the resistivity, the relative permeability, the order, the
        # frequency in Hz and the three surfaces in Ω/m. Measured: 2.4e-14.
        values = [
            (62.0e-3, 65.0e-3, 7.0e-8, 10.0, 1, 50.0,
                complex(5.8821735037835256e-5, 1.0371099434188047e-5),
                complex(5.8808016218177365e-5, 9.904296433585952e-6),
                complex(5.81890343540421e-5, -4.672717825311791e-6)),
            (62.0e-3, 65.0e-3, 7.0e-8, 10.0, 5, 50.0,
                complex(5.809610929418076e-5, 1.6429650074515897e-5),
                complex(5.8131210964069594e-5, 1.596841785890669e-5),
                complex(5.749445229625907e-5, 1.4540664775631747e-6)),
            (62.0e-3, 65.0e-3, 7.0e-8, 10.0, 20, 1.0e6,
                complex(0.004265695218963713, 0.004267537014144765),
                complex(0.004071534562008365, 0.0040705573459651175),
                complex(3.065047929249126e-34, -1.310271962151012e-33)),
            (40.0e-3, 40.2e-3, 2.8e-8, 1.0, 1, 50.0,
                complex(0.0005545742617887872, 2.4565491139060058e-5),
                complex(0.0005545742847041231, 2.4564971158796086e-5),
                complex(0.0005545742658824611, 2.4408543236448823e-5)),
            (62.0e-3, 65.0e-3, 1.7e-7, 300.0, 3, 1.0e7,
                complex(0.11518080931048041, 0.11518433037961189),
                complex(0.10987133083380493, 0.10986813043983522),
                complex(0.0, 0.0)),
            (0.0, 19.5e-3, 3.365e-8, 1.0, 1, 50.0, complex(0.0, 0.0),
                complex(2.61923294994256e-5, 2.9381241181415082e-5), complex(0.0, 0.0)),
            (0.0, 19.5e-3, 3.365e-8, 1.0, 5, 1.0e6, complex(0.0, 0.0),
                complex(0.0029814388526867193, 0.002975209109192979), complex(0.0, 0.0)),
            (0.0, 10.0e-3, 1.7e-7, 300.0, 2, 1.0e7, complex(0.0, 0.0),
                complex(0.7142780412690753, 0.714142926211519), complex(0.0, 0.0))
        ]
        for (r_in, r_ex, ρ, μ_r, n, f, inner, outer, transfer) in values
            surface = schelkunoff(n, r_in, r_ex, ρ, μ_r, im * (2π * f))
            @test isapprox(surface.inner, inner; rtol = 1.0e-12)
            @test isapprox(surface.outer, outer; rtol = 1.0e-12)
            @test isapprox(surface.transfer, transfer; rtol = 1.0e-12)
        end
    end

    @testset "a tube at an order has a nonnegative order and the checks of a net current" begin
        s = im * (2π * 50.0)
        @test_throws DomainError schelkunoff(-1, 20.5e-3, 21.5e-3, 1.7241e-8, 1.0, s)
        @test_throws DomainError schelkunoff(2, 2.0e-3, 1.0e-3, 1.7241e-8, 1.0, s)
        @test_throws DomainError schelkunoff(2, 0.0, 1.0e-3, 0.0, 1.0, s)
        @test_throws DomainError schelkunoff(2, 0.0, 1.0e-3, 1.7241e-8, -1.0, s)
        @test_throws DomainError schelkunoff(2, 0.0, 1.0e-3, 1.7241e-8, 1.0, 0.0im)
        @test @inferred(schelkunoff(20, walls[1]..., s)).transfer isa ComplexF64
    end
end

@testset "Wedepohl and Wilcox internal approximation" begin
    wedepohl = LineCableModels.Engine.InternalImpedance.Formula(Val(:wedepohl1973))
    ρ, μ_r = 1.7241e-8, 1.0
    for (r2, r3) in ((20.0e-3, 21.0e-3), (10.0e-3, 21.0e-3))
        Δ = r3 - r2
        for f in (50.0, 1.0e4, 1.0e6)
            s = im * 2π * f
            m = sqrt(s * 4π * 1e-7 * μ_r / ρ)
            x = m * Δ
            inner = (ρ * m / (2π * r2)) * (cosh(x) / sinh(x)) -
                    ρ / (2π * r2 * (r2 + r3))
            transfer = (ρ * m / (π * (r2 + r3))) * (1 / sinh(x))
            outer = (ρ * m / (2π * r3)) * (cosh(x) / sinh(x)) +
                    ρ / (2π * r3 * (r2 + r3))
            surface = wedepohl(r2, r3, ρ, μ_r, s)
            @test isapprox(surface.inner, inner; rtol = 1.0e-12)
            @test isapprox(surface.outer, outer; rtol = 1.0e-12)
            @test isapprox(surface.transfer, transfer; rtol = 1.0e-12)
        end
        resistance = ρ / (π * Δ * (r2 + r3))
        surface = wedepohl(r2, r3, ρ, μ_r, im * 2π * 1.0e-6)
        @test isapprox(surface.inner, resistance; rtol = 1.0e-6)
        @test isapprox(surface.outer, resistance; rtol = 1.0e-6)
        @test isapprox(surface.transfer, resistance; rtol = 1.0e-6)
    end
    r_ex = 10.5e-3
    for f in (50.0, 1.0e4, 1.0e6)
        s = im * 2π * f
        m = sqrt(s * 4π * 1e-7 * μ_r / ρ)
        x = 0.777 * m * r_ex
        outer = (ρ * m / (2π * r_ex)) * (cosh(x) / sinh(x)) +
                0.356 * ρ / (π * r_ex^2)
        surface = wedepohl(0.0, r_ex, ρ, μ_r, s)
        @test isapprox(surface.outer, outer; rtol = 1.0e-12)
        @test iszero(surface.inner)
        @test iszero(surface.transfer)
    end

    copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
    polymer = Material(:insulator, 1.0e14, 2.3)
    thick = build(CableDesign, "thick_tube",
        Terminal(:core, Region(:core_metal, Circle(9.0e-3), copper)),
        Region(:core_insulation, Offset(1.0e-3), polymer),
        Terminal(:sheath, Region(:sheath_copper, Annulus(10.0e-3, 21.0e-3), copper)))
    system = build(CableSystem, [thick], [Placement(0.0, -1.0)];
        connections = [(core = 1, sheath = 2)], id = "thick_tube", line_length = 1.0)
    problem = LineParametersProblem(system, homogeneous(rho = 100);
        frequencies = [50.0, 1.0e4, 1.0e6])
    formulation = CoaxialFormulation(internal_impedance = :wedepohl1973)
    result = @test_logs (:warn,
        r"^the formula :wedepohl1973 holds for wall in 0\.0 \.\. 0\.125; the computation has wall in 0\.35") compute(
        problem, formulation)
    @test result isa LineParameters

    thin = build(CableDesign, "thin_tube",
        Terminal(:core, Region(:core_metal, Circle(19.0e-3), copper)),
        Region(:core_insulation, Offset(1.0e-3), polymer),
        Terminal(:sheath, Region(:sheath_copper, Annulus(20.0e-3, 21.0e-3), copper)))
    system = build(CableSystem, [thin], [Placement(0.0, -1.0)];
        connections = [(core = 1, sheath = 2)], id = "thin_tube", line_length = 1.0)
    problem = LineParametersProblem(system, homogeneous(rho = 100);
        frequencies = [50.0, 1.0e4, 1.0e6])
    @test_logs compute(problem, formulation)

    solid = build(CableDesign, "solid_conductor",
        Terminal(:core, Region(:core_metal, Circle(10.5e-3), copper)))
    system = build(CableSystem, [solid], [Placement(0.0, -1.0)];
        connections = [(core = 1,)], id = "solid_conductor", line_length = 1.0)
    problem = LineParametersProblem(system, homogeneous(rho = 100);
        frequencies = [50.0, 1.0e4, 1.0e6])
    result = @test_logs compute(problem, formulation)
    @test result isa LineParameters

    @test @inferred(wedepohl(20.0e-3, 21.0e-3, ρ, μ_r, im * 2π * 50.0)) isa NamedTuple
    @test @inferred(wedepohl(0.0, 10.5e-3, ρ, μ_r, im * 2π * 50.0)) isa NamedTuple
end
