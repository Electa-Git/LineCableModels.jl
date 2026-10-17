@testitem "Engine / internal shunt / independent Green and integrated-face controls" tags=[:unit, :engine] begin
    using LinearAlgebra
    E = LineCableModels.Engine
    boundary = E.ShuntModel.Formula(:boundary)
    include(joinpath(pkgdir(LineCableModels), "test", "support", "internal_shunt.jl"))
    integration = (rtol = 1e-10, atol = 1e-12, maxevals = 100_000)
    h = annular_group(1.0, 3.0, 1.0)
    k = E.ShuntModel._shunt_kernel_coefficients(h, 256)
    kernel(z, i, s, j, g, coefficients) =
        only(kernel_matrix(boundary, [z], [i], [s], [j], g, coefficients))
    @testset "Local Green function controls" begin
        # Independent Dirichlet-annulus separated solution, not the reflection form.
        z, s=1.4cis(0.3), 2.2cis(1.2)
        x, y, L=log(abs(z)), log(abs(s)), log(3.0)
        exact=min(x, y)*(L-max(x, y))/L+sum(1:200) do m
            2sinh(m*min(x, y))*sinh(m*(L-max(x, y)))/(m*sinh(m*L))*cos(m*(angle(z)-angle(s)))
        end
        @test kernel(z, 1, s, 1, h, k) ≈ exact rtol=1e-12
        @test kernel(s, 1, z, 1, h, k) ≈ exact rtol=1e-12
        @test reference_kernel(z, 1, s, 1, h; modes = 200) ≈ exact rtol=1e-12
        @test abs(kernel(cis(0.3), 1, s, 1, h, k))<1e-12
        @test abs(kernel(3cis(0.3), 1, s, 1, h, k))<1e-12
        h2 = annular_group(1.0, 3.0, 2.0)
        @test kernel(1.4cis(0.2), 1, 2.1cis(1.1), 1, h2,
            E.ShuntModel._shunt_kernel_coefficients(h2, 256)) ≈
              kernel(1.4cis(0.2), 1, 2.1cis(1.1), 1, h, k)/2
        targets, sources = [1.3cis(0.2), 1.4cis(2.1)], [2.1cis(1.1), 2.5cis(-0.4)]
        @test kernel_matrix(boundary, targets, [1, 1], sources, [1, 1], h, k) ≈
              [reference_kernel(z, 1, s, 1, h; modes = 400) for z in targets, s in sources]
        # Splitting a layer at an interface between equal permittivities leaves the kernel
        # unchanged.
        one = layered_group([(ri = 1.0, ro = 3.0, epsilon = 3.0)])
        divided = layered_group([(ri = 1.0, ro = 1.7, epsilon = 3.0),
            (ri = 1.7, ro = 3.0, epsilon = 3.0)])
        targets, sources = [1.2cis(0.2), 2.4cis(1.1), 1.69cis(-0.3)], [1.5cis(0.9), 2.9cis(2.0)]
        layers(points) = [abs(point) < 1.7 ? 1 : 2 for point in points]
        reference = kernel_matrix(boundary, targets, ones(Int, 3), sources, ones(Int, 2), one,
            E.ShuntModel._shunt_kernel_coefficients(one, 512))
        @test kernel_matrix(boundary, targets, layers(targets), sources, layers(sources), divided,
            E.ShuntModel._shunt_kernel_coefficients(divided, 512)) ≈ reference rtol=1e-12

        # End-to-end charge extraction control: a complete metal annulus must
        # recover the two exact coaxial capacitances and perfect shielding.
        theta = 2pi .* (0:255) ./ 256
        targets = [1.8cis.(theta); 2.2cis.(theta)]
        sources = [1.98cis.(theta); 2.02cis.(theta)]
        matrix = kernel_matrix(boundary, targets, ones(Int, 512), sources, ones(Int, 512), h, k)
        q = matrix \
            hcat(-E.ShuntModel._shunt_core_voltage.(targets, 1, Ref(k)),
                ones(length(targets)))
        screen = vec(sum(q; dims = 1))
        core = [1/k.Rtotal, 0.0] -
               vec(transpose(E.ShuntModel._shunt_core_voltage.(sources, 1, Ref(k)))*q)
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
            parts = E.ShuntModel._shunt_kernel(z, 1, s, 1, k; split_images = true) -
                    k.direct[1, 1]*log(abs(z-s)) -
                    k.inner[1]/k.epsilon[1]*log(abs(1/conj(z)-s)) -
                    k.outer[1]/k.epsilon[1]*log(abs(9/conj(z)-s))
            @test parts ≈ E.ShuntModel._shunt_kernel(z, 1, s, 1, k) atol=2e-14
        end
        # Direct and coincident nearest image produce the interface coefficient: a source on
        # the interface sees the mean permittivity, from either layer.
        for layer in (1, 2)
            difference = kernel(3.0+0im, 1, 3cis(1e-7), layer, h, k) -
                         kernel(3.0+0im, 1, 3cis(2e-7), layer, h, k)
            @test difference ≈ (2/5)*log(2) rtol=1e-5
        end
        # The charge measure contains no extra face-length Jacobian. Scaling all
        # radii and target-source positions leaves the potential columns unchanged.
        targets = [1.4cis(0.1), 2.5cis(0.3)]
        face = merge(face, (; layer = 1))
        original = tape_columns(boundary, targets, [1, 1], [(face,)], h, k; integration).columns
        scaled_h = annular_group(10.0, 30.0, 1.0;
            right = [(ri = 30.0, ro = 40.0, epsilon = 4.0)])
        scaled_face = merge(face, (
            shape = (phi = 0.0, mid = 10face.shape.mid, half = 10face.shape.half),
            points = 10face.points))
        scaled = tape_columns(boundary, 10targets, [1, 1], [(scaled_face,)], scaled_h,
            E.ShuntModel._shunt_kernel_coefficients(scaled_h, 512); integration).columns
        @test scaled ≈ original atol=1e-12
    end
end

@testitem "Engine / internal shunt / layered Green function controls" tags=[:unit, :engine] begin
    using LinearAlgebra
    E = LineCableModels.Engine
    boundary = E.ShuntModel.Formula(:boundary)
    include(joinpath(pkgdir(LineCableModels), "test", "support", "internal_shunt.jl"))
    g = layered_group([(ri = 1.0, ro = 1.5, epsilon = 2.0), (ri = 1.5, ro = 2.4, epsilon = 5.0),
        (ri = 2.4, ro = 3.2, epsilon = 1.0)])
    k = E.ShuntModel._shunt_kernel_coefficients(g, 256)
    first_points = [1.2cis(0.3), 1.9cis(-0.8), 2.8cis(2.1)]
    second_points = [1.4cis(1.0), 2.1cis(0.1), 2.6cis(-2.5)]
    layers = [1, 2, 3]
    forward = kernel_matrix(boundary, first_points, layers, second_points, layers, g, k)
    # Every pair of layers against the independent transfer-matrix reference.
    reference = [reference_kernel(z, i, s, j, g; modes = 600)
                 for (z, i) in zip(first_points, layers), (s, j) in zip(second_points, layers)]
    @test forward ≈ reference rtol=1e-10
    # Reciprocity across layers, to roundoff.
    backward = kernel_matrix(boundary, second_points, layers, first_points, layers, g, k)
    @test maximum(abs, forward-transpose(backward)) <= 1e-13maximum(abs, forward)
    # A line charge outside a dielectric cylinder, with the anchor and the reference far.
    # Inside, it gives the transmitted logarithm with 2/(epsilon_1+epsilon_2). Outside, it
    # gives the direct logarithm and the image with (epsilon_1-epsilon_2)/(epsilon_1+epsilon_2).
    # Classical image solution of the line charge and the dielectric cylinder, cited from
    # Wang (2010), College Physics 29(12), https://dxwl.bnu.edu.cn/EN/Y2010/V29/I12/12.
    e1, e2 = 4.0, 1.5
    cylinder = layered_group([(ri = 1e-6, ro = 1.0, epsilon = e1), (ri = 1.0, ro = 1e6, epsilon = e2)])
    coefficients = E.ShuntModel._shunt_kernel_coefficients(cylinder, 512)
    # The zero mode, from the log-radius sum R(r) of the two layers.
    R(r) = r < 1 ? log(r/1e-6)/e1 : log(1/1e-6)/e1+log(r)/e2
    total = R(1e6)
    zero_mode(z, s) = min(R(abs(z)), R(abs(s)))*(total-max(R(abs(z)), R(abs(s))))/total
    source = 2.0cis(0.4)
    for (z, layer, expected) in (
        (0.5cis(-0.7), 1, 2/(e1+e2)*log(2.0/abs(0.5cis(-0.7)-source))),
        (3.0cis(1.9), 2, (log(3.0/abs(3.0cis(1.9)-source)) +
                          (e1-e2)/(e1+e2)*log(abs(1-1/(conj(3.0cis(1.9))*source))))/e2))
        actual = only(kernel_matrix(boundary, [z], [layer], [source], [2], cylinder, coefficients)) -
                 zero_mode(z, source)
        @test actual ≈ expected rtol=1e-10
    end
end
