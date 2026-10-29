# The targets come from the schelkunoff1934 and ametani1980 formulas of commit 989af96c of
# this repository, with s = j 2π f. They are the equivalent tubes
# and the dielectric layers of the coaxial references.
@testset "Engine impedances" begin
    Engine = LineCableModels.Engine
    impedance = Engine.impedance
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
            surface = impedance(schelkunoff, tubes[part]..., im * (2π * f))
            @test isapprox(surface.inner, inner; rtol = 1.0e-13)
            @test isapprox(surface.outer, outer; rtol = 1.0e-13)
            @test isapprox(surface.transfer, transfer; rtol = 1.0e-13)
        end
    end

    @testset "the magnetic term across the layers" begin
        for (part, f, term) in terms
            @test isapprox(impedance(ametani, layers[part]..., im * (2π * f)), term;
                rtol = 1.0e-13)
        end
        @test iszero(impedance(ametani, 10.5e-3, 10.5e-3, 1.0, im * (2π * 50.0)))
    end

    @testset "a tube has ordered radii, a conductive material and a frequency" begin
        s = im * (2π * 50.0)
        @test_throws DomainError impedance(schelkunoff, 2.0e-3, 1.0e-3, 1.7241e-8, 1.0, s)
        @test_throws DomainError impedance(schelkunoff, 0.0, 1.0e-3, 0.0, 1.0, s)
        @test_throws DomainError impedance(schelkunoff, 0.0, 1.0e-3, 1.7241e-8, -1.0, s)
        @test_throws DomainError impedance(schelkunoff, 0.0, 1.0e-3, 1.7241e-8, 1.0, 0.0im)
    end

    @testset "the impedances are inferred" begin
        s = im * (2π * 50.0)
        @test @inferred(impedance(schelkunoff, tubes["sheath"]..., s)).outer isa ComplexF64
        @test @inferred(impedance(ametani, layers["insulation"]..., s)) isa ComplexF64
        s32 = ComplexF32(s)
        @test real(@inferred(impedance(
            schelkunoff, 0.0f0, 1.0f-3, 1.7241f-8, 1.0f0, s32)).outer) > 0
        @test iszero(@inferred(impedance(ametani, 1.0f-3, 1.0f-3, 1.0f0, s32)))
    end
end
