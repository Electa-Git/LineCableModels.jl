function cases(; r_wire = 0.0425, t_ins = 1.0e-3, eps_ins = 2.3, rho_cu = 1.7241e-8,
        rho_earth_w = 0.1, x2 = 1.0, r_strand = 1.5e-3, w_strand = 3.0e-3,
        t_xlpe = 8.0e-3, mu_steel = 36.6326, rho_earth_c = 100.0, eps_earth_c = 10.0,
        temperature = 20.0, frequencies = last(first(REFERENCES)).frequencies[1:10:end],
        connections = [(core = 1,), (core = 2,)], split_wire = false)
    copper = Material(:conductor, rho_cu, 1, 1, 20, 0.00393)
    library_copper = Material(:conductor, rho_cu, 1, 0.999994, 20, 0.00393)
    insulation = Material(:insulator, 1.97e14, eps_ins)
    polymer = Material(:insulator, 1.0e14, 2.3)
    semicon = Material(:semicon, 1000.0, 1000.0)
    steel = Material(:conductor, 1.7475e-7, 1.0, mu_steel, 20.0, 0.0045)
    air = material(kind = :insulator, rho = Inf)
    metal = Terminal(:core, Region(:core_metal, Circle(r_wire), library_copper))
    metal2 = split_wire ?
             Terminal(:core,
        Region(:core_metal,
            Circle(measurement(nominal(r_wire), Measurements.uncertainty(r_wire))),
            library_copper)) : metal
    wires = [
        "two_bare_wires" => (build(CableDesign, "bare_wire", metal),
            build(CableDesign, "bare_wire", metal2)),
        "two_insulated_wires" => (
            build(CableDesign, "insulated_wire", metal,
                Region(:core_insulation, Offset(t_ins), insulation)),
            build(CableDesign, "insulated_wire", metal2,
                Region(:core_insulation, Offset(t_ins), insulation)))]
    result = [case => LineParametersProblem(
                  build(CableSystem, [a, b],
                      [Placement(0.0, -1.0), Placement(x2, -1.0)];
                      connections, id = case, line_length = 1.0),
                  homogeneous(rho = rho_earth_w); frequencies, temperature)
              for (case, (a, b)) in wires]
    for (case, strand) in (("coaxial_circular_strands", Circle(r_strand)),
        ("coaxial_rectangular_strands", Rectangle(w_strand, w_strand)))
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
            Region(:core_insulation, Offset(t_xlpe), polymer),
            Region(:core_semicon_outer, Offset(1.0e-3), semicon),
            Terminal(:sheath, Region(:sheath_copper, Annulus(20.5e-3, 21.5e-3), copper)),
            Region(:sheath_bedding, Offset(3.0e-3), polymer),
            Terminal(:armour, Enclosure(Region(:armour_fill, Offset(4.0e-3), air), armour)),
            Region(:armour_jacket, Offset(4.0e-3), polymer))
        system = build(CableSystem, [design], [Placement(0.0, -1.0)];
            connections = [(core = 1, sheath = 2, armour = 3)], id = case,
            line_length = 1000.0)
        push!(result,
            case => LineParametersProblem(system,
                homogeneous(rho = rho_earth_c, eps_r = eps_earth_c); frequencies, temperature))
    end
    return result
end

@testset "Measurement propagation" begin
    inputs = (
        r_wire = (0.0425, 1.0e-4, (1, 2), (1, -1), (1 / 2, -1 / 2)),
        t_ins = (1.0e-3, 2.0e-5, (2,), (1, -1), (1 / 2, -1 / 2)),
        eps_ins = (2.3, 0.05, (2,), (1, -1), (1 / 2, -1 / 2)),
        rho_cu = (1.7241e-8, 1.7e-10, (1, 2, 3, 4), (1, -1), (1 / 2, -1 / 2)),
        rho_earth_w = (0.1, 0.01, (1, 2), (1, -1), (1 / 2, -1 / 2)),
        x2 = (1.0, 1.0e-3, (1, 2), (1, -1), (1 / 2, -1 / 2)),
        r_strand = (1.5e-3, 1.0e-5, (3,), (0, -1, -2), (3 / 2, -2, 1 / 2)),
        w_strand = (3.0e-3, 2.0e-5, (4,), (0, -1, -2), (3 / 2, -2, 1 / 2)),
        mu_steel = (36.6326, 0.5, (3, 4), (1, -1), (1 / 2, -1 / 2)),
        rho_earth_c = (100.0, 10.0, (3, 4), (1, -1), (1 / 2, -1 / 2)),
        eps_earth_c = (10.0, 0.5, (3, 4), (1, -1), (1 / 2, -1 / 2)),
        temperature = (20.0, 2.0, (1, 3), (1, -1), (1 / 2, -1 / 2)))
    nominal_inputs = (t_xlpe = (8.0e-3, 1.0e-4, (3, 4)),)
    formulation = CoaxialFormulation()
    for (index, (case, reference)) in enumerate(REFERENCES)
        @testset "$case" begin
            frequencies = reference.frequencies[1:10:end]
            deterministic = compute(
                LineParametersProblem(reference.system, reference.earth;
                    frequencies, temperature = reference.temperature),
                formulation)
            sources = (;
                (name => measurement(v, σ)
            for (name, (v, σ, which)) in pairs(merge(inputs, nominal_inputs))
            if index in which)...)
            singles = (;
                (name => compute(
                     last(cases(; frequencies, (; name => source)...)[index]),
                     formulation) for (name, source) in pairs(sources))...)
            for (name, result) in pairs(singles)
                @testset "$name" begin
                    for (Q, Q0) in
                        ((result.Z, deterministic.Z), (result.Y, deterministic.Y))
                        @test maximum(maximum(abs, nominal.(Q[:, :, k]) - Q0[:, :, k]) /
                                      maximum(abs(Q0[i, i, k]) for i in axes(Q0, 1))
                        for k in eachindex(frequencies)) <= 1e-14
                    end
                end
            end
            for (name, (v, σ, which, offsets, weights)) in pairs(inputs)
                index in which || continue
                @testset "$name" begin
                    source = sources[name]
                    result = singles[name]
                    δ = 1e-5 * v
                    perturbed = [compute(
                                     last(cases(; frequencies,
                                         (; name => v + o * δ)...)[index]),
                                     formulation)
                                 for o in offsets]
                    differences = (
                        sum(w * p.Z for (p, w) in zip(perturbed, weights)) / δ,
                        sum(w * p.Y for (p, w) in zip(perturbed, weights)) / δ)
                    for (Q, Q0, D) in zip((result.Z, result.Y),
                        (deterministic.Z, deterministic.Y), differences)
                        scale = maximum(abs(Q0[i, i, k])
                        for i in axes(Q0, 1) for k in axes(Q0, 3))
                        @test maximum(
                            (abs(Measurements.derivative(part(Q[i]), source) - part(D[i])) /
                             abs(part(D[i]))
                            for i in eachindex(Q)
                            for part in (real, imag)
                            if abs(part(D[i])) * σ > 1e-9 * scale);
                            init = 0.0) <= 1e-5
                    end
                end
            end
            problem = last(cases(; frequencies, sources...)[index])
            together = @inferred compute(problem, CoaxialFormulation())
            for (Q, quantities) in ((together.Z, [single.Z for single in singles]),
                (together.Y, [single.Y for single in singles]))
                rss = [sqrt(sum(Measurements.uncertainty(part(single[i]))^2
                       for single in quantities))
                       for i in eachindex(Q), part in (real, imag)]
                @test maximum(
                    (abs(Measurements.uncertainty(part(Q[i])) - rss[i, j]) / rss[i, j]
                    for i in eachindex(Q)
                    for (j, part) in enumerate((real, imag))
                    if rss[i, j] > 0); init = 0.0) <= 1e-12
            end
            @test length(Measurements.uncertainty_components(real(together.Z[1, 1, 1]))) ==
                  count(single -> Measurements.uncertainty(real(single.Z[1, 1, 1])) > 0,
                singles)
        end
    end

    @testset "Repeated source" begin
        connections = [(core = 1,), (core = 1,)]
        one_source = compute(
            last(first(cases(; connections,
                r_wire = measurement(0.0425, 1e-4)))), formulation)
        two_sources = compute(
            last(first(cases(; connections, split_wire = true,
                r_wire = measurement(0.0425, 1e-4)))),
            formulation)
        @test maximum(abs(Measurements.uncertainty(part(one_source.Z[1, 1, k])) /
                          Measurements.uncertainty(part(two_sources.Z[1, 1, k])) -
                          sqrt(2)) /
                      sqrt(2) for k in axes(one_source.Z, 3) for part in (real, imag)) <=
              1e-12
        @test all(length(Measurements.uncertainty_components(part(one_source.Z[1, 1, k]))) ==
                  1 &&
                  length(Measurements.uncertainty_components(
                      part(two_sources.Z[1, 1, k]))) == 2
        for k in axes(one_source.Z, 3) for part in (real, imag))
    end

    @testset "Float32" begin
        wire = build(CableDesign, "float32",
            Terminal(:core,
                Region(:core_metal, Circle(measurement(1.5f-3, 1.0f-5)),
                    material(kind = :conductor, rho = 1.7241f-8))))
        @test wire isa CableDesign{Measurement{Float32}}
        system = build(CableSystem, [wire], [Placement(0.0f0, -1.0f0)];
            connections = [(core = 1,)], id = "float32", line_length = 1.0f0)
        problem = LineParametersProblem(system, homogeneous(rho = 100);
            frequencies = [50.0])
        @test problem isa LineParametersProblem{Measurement{Float64}}
        result = compute(problem, formulation)
        @test result isa LineParameters{Measurement{Float64}}
        @test all(isfinite, result.Z) && all(isfinite, result.Y)
    end

    @testset "JET" begin
        Engine = LineCableModels.Engine
        sources = (;
            (name => measurement(v, σ)
        for (name, (v, σ, which)) in pairs(inputs) if 2 in which)...)
        problem = last(cases(; sources...)[2])
        input = Engine.CoaxialInput(problem, formulation.options)
        buffers = Engine.initialize_buffers(formulation, problem, input)
        n = length(input.cable)
        c = length(input.radius)
        m = length(input.kept)
        nf = length(problem.frequencies)
        Z = Array{Complex{Measurement{Float64}}, 3}(undef, m, m, nf)
        Y = Array{Complex{Measurement{Float64}}, 3}(undef, m, m, nf)
        trace = (Zin = Array{Complex{Measurement{Float64}}, 3}(undef, n, n, nf),
            Pin = Array{Complex{Measurement{Float64}}, 3}(undef, n, n, nf),
            Zg = Array{Complex{Measurement{Float64}}, 3}(undef, c, c, nf),
            Pg = Array{Complex{Measurement{Float64}}, 3}(undef, c, c, nf))
        options = LineCableModels.ComputationOptions(Engine.LCMBackend; trace = true)
        JET.@test_opt target_modules=(LineCableModels,) Engine.solve!(
            Z, Y, trace, problem, input, formulation, buffers, options)
    end
end
