@testset "Engine laws" begin
    Engine = LineCableModels.Engine
    select = LineCableModels.select
    copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)

    @testset "the linear resistivity law" begin
        law = Engine.Resistivity.Formula(Val(:linear))
        @test law(copper, copper.T0) == copper.rho
        @test law(copper, copper.T0 + 70) ≈
              copper.rho * (1 + 70 * copper.alpha) rtol = 1.0e-15
        @test_throws DomainError law(copper, copper.T0 + 150)
        @test_throws DomainError law(copper, copper.T0 - 150)
    end

    @testset "the constant properties of the earth" begin
        layer = last(homogeneous(rho = 100, eps_r = 10, mu_r = 2).layers)
        law = Engine.EarthProperties.Formula(Val(:constant))
        @test law(layer, 50.0) == layer
    end

    @testset "the deepest layer of a stack" begin
        model = Engine.EarthProperties.Formula(Val(:bottommost))
        upper = EarthLayer(100.0, 10.0, 1.0, 5.0)
        lower = EarthLayer(1000.0, 5.0, 1.0)
        @test model((upper, lower), 50.0) === lower
    end

    @testset "the lossless admittivity of insulators and semicons" begin
        f = 1.0e3
        polymer = Material(:insulator, 1.0e14, 2.3)
        semicon = Material(:semicon, 1000.0, 1000.0)
        insulation = Engine.InsulationAdmittance.Formula(Val(:lossless))
        @test insulation(polymer, f) ==
              im * 2π * f * 8.8541878128e-12 * polymer.eps_r
        semiconducting = Engine.SemiconAdmittance.Formula(Val(:lossless))
        @test semiconducting(semicon, f) ==
              im * 2π * f * 8.8541878128e-12 * semicon.eps_r
    end

    @testset "the lossy admittivity of insulators and semicons" begin
        insulating = Engine.InsulationAdmittance.Formula(Val(:lossy))
        semiconducting = Engine.SemiconAdmittance.Formula(Val(:lossy))
        polymer = Material(:insulator, 1.0e14, 2.3; tan_delta = 1.0e-3)
        semicon = Material(:semicon, 1000.0, 1000.0; tan_delta = 0.05)
        for f in (50.0, 1.0e3, 1.0e6),
            (law, material) in ((insulating, polymer), (semiconducting, semicon))

            ω = 2π * f
            ε = 8.8541878128e-12 * material.eps_r
            @test @inferred(law(material, f)) ≈
                  complex(1 / material.rho + ω * ε * material.tan_delta, ω * ε) rtol = 1.0e-15
        end
        polymer32 = Material(:insulator, 1.0f14, 2.3f0; tan_delta = 1.0f-3)
        semicon32 = Material(:semicon, 1000.0f0, 1000.0f0; tan_delta = 0.05f0)
        @test real(@inferred(insulating(polymer32, 50.0f0))) > 0
        @test real(@inferred(semiconducting(semicon32, 50.0f0))) > 0
    end

    @testset "the lossy admittivity without loss is the lossless one" begin
        insulator = Material(:insulator, Inf, 2.3)
        semicon = Material(:semicon, Inf, 1000.0)
        for f in (50.0, 1.0e3, 1.0e6)
            @test Engine.InsulationAdmittance.Formula(Val(:lossy))(insulator, f) ===
                  Engine.InsulationAdmittance.Formula(Val(:lossless))(insulator, f)
            @test Engine.SemiconAdmittance.Formula(Val(:lossy))(semicon, f) ===
                  Engine.SemiconAdmittance.Formula(Val(:lossless))(semicon, f)
        end
    end

    @testset "a lossy insulation gives the shunt admittance its conductance" begin
        polymer = Material(:insulator, 1.0e14, 2.3; tan_delta = 1.0e-3)
        design = build(CableDesign, "insulated_wire",
            Terminal(:core, Region(:core_metal, Circle(0.0425), copper)),
            Region(:core_insulation, Offset(1.0e-3), polymer))
        system = build(CableSystem, [design, design],
            [Placement(0.0, -1.0), Placement(1.0, -1.0)];
            connections = [(core = 1,), (core = 2,)], id = "two_insulated_wires",
            line_length = 1.0)
        problem = LineParametersProblem(system, homogeneous(rho = 0.1);
            frequencies = 10.0 .^ range(-1, 7; length = 101))
        lossless = compute(problem, CoaxialFormulation())
        lossy = compute(problem, CoaxialFormulation(insulation_admittance = :lossy))
        @test all(>(0), real.(lossy.Y[1, 1, :]))
        @test all(real.(lossy.Y[1, 1, :]) .> real.(lossless.Y[1, 1, :]))
    end

    @testset "a slot takes an identifier, an identifier with options or a formula" begin
        Formula = Engine.Resistivity.Formula
        law = Formula(Val(:linear))
        @test select(CoaxialFormulation, Val(:resistivity), :linear) === law
        @test select(CoaxialFormulation, Val(:resistivity), :linear => (;)) === law
        @test select(CoaxialFormulation, Val(:resistivity), law) === law
        @test_throws "the formulas are (:linear,)" select(CoaxialFormulation, Val(:resistivity), :quadratic)
        @test_throws "the formulas are (:martins2026, :ideal)" select(
            CoaxialFormulation, Val(:earth_impedance), :carson1926)
        @test_throws ArgumentError select(CoaxialFormulation, Val(:resistivity), :linear => (rate = 1.0,))
    end

    @testset "soil laws against 9676ae6c" begin
        for (id, rho, f, ρ, ε_r) in (
            (:alipio2014, 1.0, 10.0, 0.9999837683191205, 33106.910159656065),
            (:alipio2014, 1.0, 100.0, 0.999943721034889, 11487.224914768634),
            (:alipio2014, 1.0, 3000.0, 0.999646928503297, 2412.410010822532),
            (:alipio2014, 1.0, 100000.0, 0.9976592621167621, 490.3669935403034),
            (:alipio2014, 1.0, 1.0e7, 0.9725659696637746, 69.51236230283646),
            (:alipio2014, 100.0, 10.0, 99.95320846585999, 9556.676349359284),
            (:alipio2014, 100.0, 100.0, 99.83794360552436, 3321.4910159656065),
            (:alipio2014, 100.0, 3000.0, 98.99164245386598, 704.2858091632659),
            (:alipio2014, 100.0, 100000.0, 93.66224962294666, 149.96254794261966),
            (:alipio2014, 100.0, 1.0e7, 55.141250891350566, 28.586746470061136),
            (:alipio2014, 10000.0, 10.0, 9866.7870940653, 2764.714727869962),
            (:alipio2014, 10000.0, 100.0, 9552.80055520738, 966.4676349359287),
            (:alipio2014, 10000.0, 3000.0, 7729.3144403340375, 211.65740827943526),
            (:alipio2014, 10000.0, 100000.0, 3388.097479581717, 51.788833451812984),
            (:alipio2014, 10000.0, 1.0e7, 408.79233962015223, 16.783669935403033),
            (:cigre2019, 1.0, 10.0, 0.9999837036336047, 32952.000792990504),
            (:cigre2019, 1.0, 100.0, 0.999943496764369, 11433.51212886542),
            (:cigre2019, 1.0, 3000.0, 0.9996455219399029, 2401.174265122654),
            (:cigre2019, 1.0, 100000.0, 0.9976499557192541, 488.12787194590857),
            (:cigre2019, 1.0, 1.0e7, 0.9724596504314155, 69.24316067706397),
            (:cigre2019, 100.0, 10.0, 99.95302207962257, 9511.999999999998),
            (:cigre2019, 100.0, 100.0, 99.83729882839879, 3306.00007929905),
            (:cigre2019, 100.0, 3000.0, 98.98766462350514, 701.0453847073092),
            (:cigre2019, 100.0, 100000.0, 93.63859910453635, 149.31677822086309),
            (:cigre2019, 100.0, 1.0e7, 55.04285007542964, 28.50910787311906),
            (:cigre2019, 10000.0, 10.0, 9866.263313197007, 2751.829927970275),
            (:cigre2019, 10000.0, 100.0, 9551.098376313275, 961.9999999999997),
            (:cigre2019, 10000.0, 3000.0, 7722.326393949987, 210.7228596579871),
            (:cigre2019, 10000.0, 100000.0, 3379.193377819101, 51.60259142968185),
            (:cigre2019, 10000.0, 1.0e7, 407.2357653975235, 16.761278719459085),
            (:datsios2019, 1.0, 10.0, 1.0000001105085483, 1806.482324679456),
            (:datsios2019, 1.0, 100.0, 0.9999997997033178, 1806.482324679456),
            (:datsios2019, 1.0, 3000.0, 0.9999897849715096, 1806.482324679456),
            (:datsios2019, 1.0, 100000.0, 0.9996549249839193, 198.8093360351198),
            (:datsios2019, 1.0, 1.0e7, 0.966618999289488, 198.37628317197277),
            (:datsios2019, 100.0, 10.0, 100.00015082363583, 214.20702847212598),
            (:datsios2019, 100.0, 100.0, 99.99972663331963, 214.20702847212598),
            (:datsios2019, 100.0, 3000.0, 99.98606020463883, 214.20702847212598),
            (:datsios2019, 100.0, 100000.0, 99.53108398606008, 32.95824634211199),
            (:datsios2019, 100.0, 1.0e7, 67.96612590826454, 29.36262422627703),
            (:datsios2019, 10000.0, 10.0, 10000.208026503413, 25.400000000000002),
            (:datsios2019, 10000.0, 100.0, 9999.622974021395, 25.400000000000002),
            (:datsios2019, 10000.0, 3000.0, 9980.807854777806, 25.400000000000002),
            (:datsios2019, 10000.0, 100000.0, 9389.850984595923, 7.5438522561975),
            (:datsios2019, 10000.0, 1.0e7, 1333.2149976875323, 4.610191996639284),
            (:longmire1975, 1.0, 10.0, 0.9996662791406514, 3.5962513195532444e6),
            (:longmire1975, 1.0, 100.0, 0.99178601388505, 1.0921026825504529e6),
            (:longmire1975, 1.0, 3000.0, 0.9799562631280022, 34189.961250530556),
            (:longmire1975, 1.0, 100000.0, 0.9658964614184634, 1575.899859379897),
            (:longmire1975, 1.0, 1.0e7, 0.9195332260359066, 85.24043635092137),
            (:longmire1975, 100.0, 10.0, 97.08201016152198, 240404.23839363284),
            (:longmire1975, 100.0, 100.0, 95.4137763294489, 23741.091487743626),
            (:longmire1975, 100.0, 3000.0, 92.59407389337954, 1196.5570382796773),
            (:longmire1975, 100.0, 100000.0, 86.69792302757189, 139.76746827820944),
            (:longmire1975, 100.0, 1.0e7, 57.15065817158074, 21.28668616460591),
            (:longmire1975, 10000.0, 10.0, 8841.763945060515, 5528.738585709243),
            (:longmire1975, 10000.0, 100.0, 8451.961431970825, 881.8688386868981),
            (:longmire1975, 10000.0, 3000.0, 7311.458061999555, 114.6450286231031),
            (:longmire1975, 10000.0, 100000.0, 4854.649116685031, 26.25867072223481),
            (:longmire1975, 10000.0, 1.0e7, 692.1595886008422, 8.754367083659156),
            (:messier1985, 1.0, 10.0, 0.9999056627269164, 169596.22404647016),
            (:messier1985, 1.0, 100.0, 0.9997017401888505, 53636.50523297825),
            (:messier1985, 1.0, 3000.0, 0.9983685423445011, 9799.180680462012),
            (:messier1985, 1.0, 100000.0, 0.9906535627148881, 1703.8822404647015),
            (:messier1985, 1.0, 1.0e7, 0.9137876334438456, 177.58822404647015),
            (:messier1985, 100.0, 10.0, 99.90574275465987, 16966.822404647017),
            (:messier1985, 100.0, 100.0, 99.7025386756762, 5370.850523297825),
            (:messier1985, 100.0, 3000.0, 98.3921505885815, 987.1180680462011),
            (:messier1985, 100.0, 100000.0, 91.37876334438455, 177.58822404647015),
            (:messier1985, 100.0, 1.0e7, 51.4545762480152, 24.958822404647016),
            (:messier1985, 10000.0, 10.0, 9906.535627148882, 1703.8822404647017),
            (:messier1985, 10000.0, 100.0, 9710.29453119354, 544.2850523297825),
            (:messier1985, 10000.0, 3000.0, 8595.404598625762, 105.91180680462011),
            (:messier1985, 10000.0, 100000.0, 5145.457624801519, 24.958822404647016),
            (:messier1985, 10000.0, 1.0e7, 958.3485114897368, 9.695882240464702),
            (:portela1999, 1.0, 10.0, 0.9999980289908099, 7529.080638282218),
            (:portela1999, 1.0, 100.0, 0.9999896560833383, 3951.317687792627),
            (:portela1999, 1.0, 3000.0, 0.9998802805227083, 1524.5722032772887),
            (:portela1999, 1.0, 100000.0, 0.9985070658661395, 571.1391732787303),
            (:portela1999, 1.0, 1.0e7, 0.960448444985388, 157.30479046451094),
            (:portela1999, 100.0, 10.0, 99.98029375337818, 7529.080638282218),
            (:portela1999, 100.0, 100.0, 99.89666665166529, 3951.317687792627),
            (:portela1999, 100.0, 3000.0, 98.81682844639778, 1524.5722032772887),
            (:portela1999, 100.0, 100000.0, 86.99308668320748, 571.1391732787303),
            (:portela1999, 100.0, 1.0e7, 19.5387678797776, 157.30479046451094),
            (:portela1999, 10000.0, 10.0, 9806.708493380069, 7529.080638282218),
            (:portela1999, 10000.0, 100.0, 9062.56602610457, 3951.317687792627),
            (:portela1999, 10000.0, 3000.0, 4550.960964931424, 1524.5722032772887),
            (:portela1999, 10000.0, 100000.0, 626.8938782627822, 571.1391732787303),
            (:portela1999, 10000.0, 1.0e7, 24.22463020203948, 157.30479046451094),
            (:scott1967, 1.0, 10.0, 0.3483373150360116, 3.784425847170951e7),
            (:scott1967, 1.0, 100.0, 0.4943106869868354, 2.187761623949552e6),
            (:scott1967, 1.0, 3000.0, 0.7122956276975418, 57079.144223655516),
            (:scott1967, 1.0, 100000.0, 0.8590135215053943, 2691.534803926922),
            (:scott1967, 1.0, 1.0e7, 0.8203515443298177, 144.54397707459296),
            (:scott1967, 100.0, 10.0, 85.90135215053957, 230144.1817408509),
            (:scott1967, 100.0, 100.0, 98.62794856312104, 22490.54605835778),
            (:scott1967, 100.0, 3000.0, 103.93579039972226, 1274.2689454696138),
            (:scott1967, 100.0, 100000.0, 90.78205301781857, 133.6595516546442),
            (:scott1967, 100.0, 1.0e7, 56.75446054085471, 20.511621788255624),
            (:scott1967, 10000.0, 10.0, 10914.403364487567, 4988.844874600121),
            (:scott1967, 10000.0, 100.0, 10139.113857366789, 824.138115013003),
            (:scott1967, 10000.0, 3000.0, 7813.938119967247, 101.4015854388886),
            (:scott1967, 10000.0, 100000.0, 4943.106869868356, 23.65919697485761),
            (:scott1967, 10000.0, 1.0e7, 2023.0191786782718, 10.375284158180126),
            (:visacro1987, 1.0, 10.0, 1.1803206356517297, 591855.7311683979),
            (:visacro1987, 1.0, 100.0, 1.0, 149697.9515029397),
            (:visacro1987, 1.0, 3000.0, 0.7827936061004801, 19650.539049126954),
            (:visacro1987, 1.0, 100000.0, 0.6081350012787179, 2422.2326700296653),
            (:visacro1987, 1.0, 1.0e7, 0.43651583224016605, 154.95866186621237),
            (:visacro1987, 100.0, 10.0, 118.03206356517295, 50375.092592550325),
            (:visacro1987, 100.0, 100.0, 100.0, 12741.362076512656),
            (:visacro1987, 100.0, 3000.0, 78.279360610048, 1672.5321255893105),
            (:visacro1987, 100.0, 100000.0, 60.813500127871784, 206.16543628387532),
            (:visacro1987, 100.0, 1.0e7, 43.651583224016605, 13.189121146327338),
            (:visacro1987, 10000.0, 10.0, 11803.206356517298, 4287.615748348634),
            (:visacro1987, 10000.0, 100.0, 10000.0, 1084.4657921829137),
            (:visacro1987, 10000.0, 3000.0, 7827.9360610048, 142.3555712204536),
            (:visacro1987, 10000.0, 100000.0, 6081.3500127871785, 17.547524498379463),
            (:visacro1987, 10000.0, 1.0e7, 4365.15832240166, 1.1225762698098527),
            (:visacro2012, 1.0, 10.0, 1.0, 1205.8188262704462),
            (:visacro2012, 1.0, 100.0, 1.0, 1205.8188262704462),
            (:visacro2012, 1.0, 3000.0, 0.9997863774181489, 310.3047037316234),
            (:visacro2012, 1.0, 100000.0, 0.9978719900919923, 77.29999999999997),
            (:visacro2012, 1.0, 1.0e7, 0.9591614657926464, 13.34518826270446),
            (:visacro2012, 100.0, 10.0, 100.0, 1205.8188262704462),
            (:visacro2012, 100.0, 100.0, 100.0, 1205.8188262704462),
            (:visacro2012, 100.0, 3000.0, 99.3875481910401, 310.3047037316234),
            (:visacro2012, 100.0, 100000.0, 94.20601399456788, 77.29999999999997),
            (:visacro2012, 100.0, 1.0e7, 44.884439441982515, 13.34518826270446),
            (:visacro2012, 10000.0, 10.0, 10000.0, 1205.8188262704462),
            (:visacro2012, 10000.0, 100.0, 10000.0, 1205.8188262704462),
            (:visacro2012, 10000.0, 3000.0, 8490.971837004241, 310.3047037316234),
            (:visacro2012, 10000.0, 100000.0, 3605.1935793666357, 77.29999999999997),
            (:visacro2012, 10000.0, 1.0e7, 274.61753619048676, 13.34518826270446)
        )
            law = Engine.EarthProperties.Formula(Val(id))
            described = law(EarthLayer(rho, 10.0, 1.0), f)
            @test described.rho ≈ ρ rtol = 1.0e-12
            @test described.eps_r ≈ ε_r rtol = 1.0e-12
        end
    end

    @testset "soil laws against the survey fixture" begin
        for (id, f, ρ, ε_r) in (
            (:alipio2014, 100.0, 99.83794360552434, 3321.4910132744126),
            (:alipio2014, 3000.0, 98.99164245386598, 704.2858086003167),
            (:alipio2014, 1000000.0, 80.99638649461713, 59.83669931513075),
            (:cigre2019, 100.0, 99.83729882839879, 3306.0000792990504),
            (:cigre2019, 3000.0, 98.98766462350514, 701.0453847073092),
            (:cigre2019, 1000000.0, 80.9350916879637, 59.61278719459085),
            (:datsios2019, 100.0, 99.99972663331963, 214.20702847212598),
            (:datsios2019, 3000.0, 99.98606020463883, 214.20702847212598),
            (:datsios2019, 1000000.0, 95.49909557985406, 29.615091164516546),
            (:longmire1975, 100.0, 95.41377632589052, 23741.091487743626),
            (:longmire1975, 3000.0, 92.59407388780323, 1196.5570382796773),
            (:longmire1975, 1000000.0, 77.85741105911478, 44.00427501200052),
            (:messier1985, 100.0, 99.7025386755556, 5370.850521117358),
            (:messier1985, 3000.0, 98.39215058793827, 987.1180676481041),
            (:messier1985, 1000000.0, 77.02090551639128, 61.62850521117358),
            (:portela1999, 100.0, 99.89666665166529, 3951.3176845795156),
            (:portela1999, 3000.0, 98.81682844639778, 1524.5722020375456),
            (:portela1999, 1000000.0, 56.03262861300133, 299.7380987338258),
            (:scott1967, 100.0, 98.62794856312104, 22490.54605835778),
            (:scott1967, 3000.0, 103.93579039972226, 1274.2689454696138),
            (:scott1967, 1000000.0, 74.81695005111546, 44.87453899331316),
            (:visacro1987, 100.0, 100.0, 12741.362076512656),
            (:visacro1987, 3000.0, 78.279360610048, 1672.5321255893105),
            (:visacro1987, 1000000.0, 51.52286445817565, 52.145382492925116),
            (:visacro2012, 100.0, 100.0, 1205.8188262704462),
            (:visacro2012, 3000.0, 99.3875481910401, 310.3047037316234),
            (:visacro2012, 1000000.0, 78.43837052811071, 31.556144962065787))
            law = Engine.EarthProperties.Formula(Val(id))
            described = law(EarthLayer(100.0, 10.0, 1.0), f)
            @test described.rho ≈ ρ rtol = 2.0e-9
            @test described.eps_r ≈ ε_r rtol = 2.0e-9
        end
    end

    @testset "soil laws and the range of a scan" begin
        design = build(CableDesign, "bare_wire",
            Terminal(:core, Region(:core_metal, Circle(0.0425), copper)))
        system = build(CableSystem, [design, design],
            [Placement(0.0, -1.0), Placement(1.0, -1.0)];
            connections = [(core = 1,), (core = 2,)], id = "two_bare_wires",
            line_length = 1.0)
        outside = LineParametersProblem(system, homogeneous(rho = 0.1);
            frequencies = [1.0, 1.0e3])
        inside = LineParametersProblem(system, homogeneous(rho = 0.1);
            frequencies = [1.0e3, 1.0e5])
        @testset ":alipio2014" begin
            formulation = CoaxialFormulation(earth_properties = :alipio2014)
            result = @test_logs (:warn, r":alipio2014") compute(outside, formulation)
            @test all(isfinite, result.Z)
            @test_logs compute(inside, formulation)
        end
        @testset ":cigre2019" begin
            formulation = CoaxialFormulation(earth_properties = :cigre2019)
            result = @test_logs (:warn, r":cigre2019") compute(outside, formulation)
            @test all(isfinite, result.Z)
            @test_logs compute(inside, formulation)
        end
        @testset ":datsios2019" begin
            formulation = CoaxialFormulation(earth_properties = :datsios2019)
            result = @test_logs (:warn, r":datsios2019") compute(outside, formulation)
            @test all(isfinite, result.Z)
            @test_logs compute(inside, formulation)
        end
        @testset ":longmire1975" begin
            formulation = CoaxialFormulation(earth_properties = :longmire1975)
            result = @test_logs (:warn, r":longmire1975") compute(outside, formulation)
            @test all(isfinite, result.Z)
            @test_logs compute(inside, formulation)
        end
        @testset ":messier1985" begin
            formulation = CoaxialFormulation(earth_properties = :messier1985)
            result = @test_logs (:warn, r":messier1985") compute(outside, formulation)
            @test all(isfinite, result.Z)
            @test_logs compute(inside, formulation)
        end
        @testset ":portela1999" begin
            formulation = CoaxialFormulation(earth_properties = :portela1999)
            result = @test_logs (:warn, r":portela1999") compute(outside, formulation)
            @test all(isfinite, result.Z)
            @test_logs compute(inside, formulation)
        end
        @testset ":scott1967" begin
            formulation = CoaxialFormulation(earth_properties = :scott1967)
            result = @test_logs (:warn, r":scott1967") compute(outside, formulation)
            @test all(isfinite, result.Z)
            @test_logs compute(inside, formulation)
        end
        @testset ":visacro1987" begin
            formulation = CoaxialFormulation(earth_properties = :visacro1987)
            result = @test_logs (:warn, r":visacro1987") compute(outside, formulation)
            @test all(isfinite, result.Z)
            @test_logs compute(inside, formulation)
        end
        @testset ":visacro2012" begin
            formulation = CoaxialFormulation(earth_properties = :visacro2012)
            result = @test_logs (:warn, r":visacro2012") compute(outside, formulation)
            @test all(isfinite, result.Z)
            @test_logs compute(inside, formulation)
        end
        above = LineParametersProblem(system, homogeneous(rho = 0.1);
            frequencies = [1.0e3, 1.0e11])
        for id in
            (:alipio2014, :cigre2019, :datsios2019, :longmire1975, :messier1985,
            :portela1999, :scott1967, :visacro1987, :visacro2012)
            law = Engine.EarthProperties.Formula(Val(id))
            @test_logs (:warn, Regex(":$id")) LineCableModels.validate(above, law)
        end
    end

    @testset "soil laws refuse options outside their domains" begin
        @test_throws ArgumentError Engine.EarthProperties.Formula(Val(:alipio2014); exponent = 1.0)
        @test_throws ArgumentError Engine.EarthProperties.Formula(Val(:portela1999); exponent = 1.0)
        @test_throws ArgumentError Engine.EarthProperties.Formula(Val(:longmire1975); corner_scale = 0.0)
        @test_throws ArgumentError Engine.EarthProperties.Formula(Val(:messier1985); epsilon_high = -1.0)
        @test_throws ArgumentError Engine.EarthProperties.Formula(
            Val(:visacro1987); normalization_frequency = -1.0)
        @test_throws ArgumentError Engine.EarthProperties.Formula(Val(:visacro2012); frequency_boundary = 0.0)
    end

    @testset "soil laws through compute" begin
        design = build(CableDesign, "bare_wire",
            Terminal(:core, Region(:core_metal, Circle(0.0425), copper)))
        system = build(CableSystem, [design, design],
            [Placement(0.0, -1.0), Placement(1.0, -1.0)];
            connections = [(core = 1,), (core = 2,)], id = "two_bare_wires",
            line_length = 1.0)
        problem = LineParametersProblem(system, homogeneous(rho = 0.1);
            frequencies = [1.0e3, 1.0e4, 1.0e5])
        for id in (:alipio2014, :cigre2019, :datsios2019, :longmire1975, :messier1985,
            :portela1999, :scott1967, :visacro1987, :visacro2012)
            result = compute(problem, CoaxialFormulation(earth_properties = id))
            @test all(isfinite, result.Z)
            @test all(isfinite, result.Y)
        end
    end

    @testset "soil laws with Float32" begin
        layer = EarthLayer(100.0f0, 10.0f0, 1.0f0)
        for id in (:alipio2014, :cigre2019, :datsios2019, :longmire1975, :messier1985,
            :portela1999, :scott1967, :visacro1987, :visacro2012)
            law = Engine.EarthProperties.Formula(Val(id))
            @test @inferred(law(layer, 1000.0f0)) isa EarthLayer{Float32}
        end
    end
end
