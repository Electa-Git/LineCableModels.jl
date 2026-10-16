@testitem "Engine / internal shunt / independent Green and integrated-face controls" tags=[:unit, :engine] begin
    using LinearAlgebra
    E = LineCableModels.Engine
    boundary = E.ShuntModel.Formula(:boundary)
    include(joinpath(pkgdir(LineCableModels), "test", "support", "internal_shunt.jl"))
    integration = (rtol = 1e-10, atol = 1e-12, maxevals = 100_000)
    # Independent Dirichlet-annulus separated solution, not the reflection form.
    h = annular_group(1.0, 3.0, 1.0)
    k = E.ShuntModel._shunt_kernel_coefficients(h, 256)
    @testset "Local Green function controls" begin
        z, s=1.4cis(0.3), 2.2cis(1.2)
        x, y, L=log(abs(z)), log(abs(s)), log(3.0)
        resolved=false
        for order in (32, 64, 128, 256, 512)
            exact=min(x, y)*(L-max(x, y))/L
            magnitude=abs(exact)
            for m in 1:order
                term=2sinh(m*min(x, y))*sinh(m*(L-max(x, y))) /
                     (m*sinh(m*L))*cos(m*(angle(z)-angle(s)))
                exact+=term;
                magnitude+=abs(term)
            end
            ratio=exp(-abs(y-x))
            tail=ratio^(order+1)/((order+1)*(1-ratio)*(1-exp(-2L)))
            coefficients=E.ShuntModel._shunt_kernel_coefficients(h, order)
            actual=reference_kernel(z, s, h, coefficients)
            budget=1e-10abs(exact)
            uncertainty=tail+(2order+20)*eps(Float64)*magnitude
            if uncertainty<=budget/4
                @test abs(actual-exact)+uncertainty<=budget
                resolved=true
            end
            @test reference_kernel(z, s, h, coefficients) ≈
                  reference_kernel(s, z, h, coefficients) atol=1e-12
            @test abs(reference_kernel(cis(0.3), s, h, coefficients))<1e-12
            @test abs(reference_kernel(3cis(0.3), s, h, coefficients))<1e-12
        end
        @test resolved
        h2 = annular_group(1.0, 3.0, 2.0)
        @test reference_kernel(1.4cis(0.2), 2.1cis(1.1), h2,
            E.ShuntModel._shunt_kernel_coefficients(h2, 256)) ≈
              reference_kernel(1.4cis(0.2), 2.1cis(1.1), h, k)/2
        targets, sources = [1.3cis(0.2), 1.4cis(2.1)], [2.1cis(1.1), 2.5cis(-0.4)]
        @test kernel_matrix(boundary, targets, sources, h, k) ≈
              [reference_kernel(z, s, h, k) for z in targets, s in sources]
        regular = [reference_kernel(z, s, h, k; regular = true) for z in targets, s in sources]
        @test kernel_matrix(boundary, targets, sources, h, k) ≈
              regular-[log(abs(z-s)) for z in targets, s in sources]
        @test isfinite(reference_kernel(2cis(0.3), 2cis(0.3), h, k; regular = true))
        layers = [(ri = 1.0, ro = 1.2, epsilon = 3.0), (ri = 1.2, ro = 2.0, epsilon = 3.0)]
        @test E.ShuntModel._shunt_load(layers, 5) ≈
              E.ShuntModel._shunt_load([(ri = 1.0, ro = 2.0, epsilon = 3.0)], 5)
        @test E.ShuntModel._shunt_load(layers, 1e-8) ≈ 3/log(2) rtol=1e-12

        # End-to-end charge extraction control: a complete metal annulus must
        # recover the two exact coaxial capacitances and perfect shielding.
        theta = 2pi .* (0:255) ./ 256
        targets = [1.8cis.(theta); 2.2cis.(theta)]
        sources = [1.98cis.(theta); 2.02cis.(theta)]
        kernel = kernel_matrix(boundary, targets, sources, h, k)
        q = kernel \
            hcat(-E.ShuntModel._shunt_core_voltage.(targets, Ref(h), Ref(k)),
                ones(length(targets)))
        screen = vec(sum(q; dims = 1))
        core = [1/k.Rtotal, 0.0] -
               vec(transpose(E.ShuntModel._shunt_core_voltage.(sources, Ref(h), Ref(k)))*q)
        capacitance = vcat(transpose(core), transpose(screen))
        left, right = 1/log(1.8), 1/log(3/2.2)
        @test capacitance ≈ [left -left; -left left+right] rtol=1e-9
    end
    @testset "Four-face tape controls" begin
        @test E.ShuntModel._shunt_junction_exponent(1.0, 1.0) ≈ 2/3
        @test E.ShuntModel._shunt_junction_exponent(1.0, 32.3) ≈
              2/pi*acos(sqrt(32.3/(2*33.3)))
        @test 0.5 < E.ShuntModel._shunt_junction_exponent(1.0, 32.3) < 2/3
        for (alpha, beta) in ((0.0, 0.0), (-1/3, -1/3), (-0.49, -1/3)),
            order in (16, 32, 64, 128)

            nodes, weights = E.ShuntModel._shunt_gauss_jacobi(order, alpha, beta)
            values = reduce(hcat, [E.ShuntModel._shunt_jacobi(t, alpha, beta, 8)
                                   for t in nodes])
            @test sum(weights) ≈ 1
            @test all(weights .> 0)
            @test dot(nodes, weights) ≈ (beta-alpha)/(alpha+beta+2) atol=1e-14
            @test values*weights ≈ [1.0; zeros(8)] atol=1e-13
        end
        nodes, weights = E.ShuntModel._shunt_gauss_jacobi(64, 0.0, 0.0)
        values = transpose(reduce(hcat, [E.ShuntModel._shunt_jacobi(t, 0.0, 0.0, 3)
                                         for t in nodes]))
        face = (; shape = (phi = 0.0, mid = 2.0, half = 0.1), alpha = 0.0, beta = 0.0, p = 3,
            nodes, weights, weighted = weights .* values, beta_norm = 1.0, terminal = 2,
            points = complex.(2 .+ 0.1 .* nodes),
            recurrence = E.ShuntModel._shunt_jacobi_recurrence(0.0, 0.0, 3))
        moments(z) = first(E.ShuntModel._shunt_log_moments(z, face; integration...,
            result = zeros(4), segments = nothing))
        @test moments(2.0+0im) ≈ [log(0.1)-1, 0, 1/3, 0] atol=2e-10
        endpoint = moments(2.1+0im)
        @test endpoint[1] ≈ log(0.2)-1 atol=2e-10
        for distance in (0.1, 1e-5)
            near = moments(2.0+im*distance)
            exact = log(hypot(0.1, distance))-1+distance/0.1*atan(0.1/distance)
            @test near[1] ≈ exact atol=2e-10
        end
        h = annular_group(1.0, 3.0, 1.0; right = [(ri = 3.0, ro = 4.0, epsilon = 4.0)])
        k = E.ShuntModel._shunt_kernel_coefficients(h, 512)
        for (z, s) in ((1.3cis(0.2), 2.1cis(1.1)), (3cis(0.2), 2.5cis(-0.4)))
            split = E.ShuntModel._shunt_kernel(z,
                s,
                h,
                k;
                split_images = true)-
            (log(abs(z-s))+k.ra0*log(abs(h.a^2/conj(z)-s))+
             k.rb0*log(abs(h.b^2/conj(z)-s)))/h.epsilon
            @test split ≈ E.ShuntModel._shunt_kernel(z, s, h, k) atol=2e-14
        end
        # Direct and coincident nearest image produce the interface coefficient.
        difference = reference_kernel(3.0+0im, 3cis(1e-7), h, k) -
                     reference_kernel(3.0+0im, 3cis(2e-7), h, k)
        @test difference ≈ (2/5)*log(2) rtol=1e-5
        # The charge measure contains no extra face-length Jacobian. Scaling all
        # radii and target-source positions leaves the potential columns unchanged.
        targets = [1.4cis(0.1), 2.5cis(0.3)]
        original = tape_columns(boundary, targets, [(face,)], h, k; integration).columns
        scaled_h = annular_group(10.0, 30.0, 1.0;
            right = [(ri = 30.0, ro = 40.0, epsilon = 4.0)])
        scaled_face = merge(face, (
            shape = (phi = 0.0, mid = 10face.shape.mid, half = 10face.shape.half),
            points = 10face.points))
        scaled = tape_columns(boundary, 10targets, [(scaled_face,)], scaled_h,
            E.ShuntModel._shunt_kernel_coefficients(scaled_h, 512); integration).columns
        @test scaled ≈ original atol=1e-12
    end
end
