@testitem "Engine / complete manuscript earth matrices reproduce the accepted snapshot" tags=[:unit] begin
    using TOML, LinearAlgebra
    const E=LineCableModels.Engine
    fixture=TOML.parsefile(joinpath(
        pkgdir(LineCableModels), "test/fixtures/reference/unified_earth_return.toml"))
    for row in fixture["cases"]
        row["model"]=="proposed_full_current" || continue
        positions=row["positions_m"]
        n=length(positions)
        geometry=E.EarthReturnGeometry(first.(positions), last.(positions), fill(row["radius_m"], n))
        workspace=E.EarthReturnWorkspace(geometry)
        s=complex(0.0, 2pi*row["frequency"])
        sigma=[0.0, inv(row["earth_resistivity_ohm_m"])]
        epsilon=fill(fixture["analytic_epsilon"], 2)
        mu=fill(fixture["analytic_mu"], 2)
        state=(jω = s, Γ = zero(s), sigma, epsilon, mu,
            gamma_medium_squared = s .* mu .* (sigma .+ s .* epsilon))
        for method in (:quad, :trapz, :cim)
            rtol=method===:quad ? 1e-9 : 1e-6
            controls=E.computation_options(E.SpectralIntegral, (method, options = (; rtol)))
            E.unified_earth!(workspace, state, controls)
            @testset "$n wires / $(row["frequency"]) Hz / $method" begin
                for (key, floor) in ((:Ze, 1e-14), (:Pe, 1e-14), (:Ye, 1e-10))
                    expected=reshape(complex.(row[string(key) * "_real"], row[string(key) * "_imag"]), n, n)
                    actual=getproperty(workspace, key)
                    for i in eachindex(actual), component in (identity, real, imag)

                        @test isapprox(component(actual[i]), component(expected[i]);
                            rtol = method===:quad ? 1e-7 : 1e-5, atol = floor)
                    end
                end
                @test workspace.Ye*workspace.Pe≈s*I rtol=1e-11 atol=1e-12*abs(s)
                @test workspace.Pe*workspace.L≈workspace.H rtol=1e-11
                @test workspace.Ze*workspace.L≈workspace.K rtol=1e-11
                # Independent frozen source kernels also validate the current
                # map before the complete matrix elimination.
                u=E.unified_earth_state(state, geometry)
                for key in (:K, :H, :L)
                    expected=reshape(complex.(row[string(key) * "_real"], row[string(key) * "_imag"]), n, n)
                    expected .*= reshape(exp.(u.scaling), 1, n)
                    for i in eachindex(expected)
                        @test isapprox(getproperty(workspace, key)[i], expected[i];
                            rtol = method===:quad ? 1e-7 : 1e-5, atol = 1e-14)
                    end
                end
            end
        end
    end
end

@testitem "Engine / manuscript geometry and reference identities" tags=[:unit] begin
    using LinearAlgebra
    const E=LineCableModels.Engine
    @test_throws DomainError E.EarthReturnGeometry([0.0], [-0.01], [0.02])
    @test_throws DomainError E.EarthReturnGeometry([0.0, 0.01], [-1.0, -1.0], [0.02, 0.02])
    geometry=E.EarthReturnGeometry([0.0, 1.0, 2.0], [1.2, -0.9, -1.4], [0.01, 0.025, 0.04])
    s=complex(0.0, 2pi*1e4)
    Γ=complex(1e-4, 2e-4)
    sigma=[0.0, 0.1]
    epsilon=8.8541878128e-12 .* [1.0, 8.0]
    mu=4pi*1e-7 .* [1.0, 3.0]
    state=(jω = s, Γ, sigma, epsilon, mu,
        gamma_medium_squared = s .* mu .* (sigma .+ s .* epsilon))
    integration=E.computation_options(E.SpectralIntegral, (
        method = :quad, options = (rtol = 1e-9,)))
    deep=E.unified_earth!(E.EarthReturnWorkspace(geometry), state, integration)
    for reference in (:interface, 3.0, :scalar)
        shifted=E.unified_earth!(E.EarthReturnWorkspace(geometry), state, integration; reference)
        @test shifted.K≈deep.K rtol=1e-12
        @test shifted.L≈deep.L rtol=1e-12
        @test shifted.Ze-deep.Ze≈Γ^2/s*(shifted.Pe-deep.Pe) rtol=1e-8 atol=1e-13
        @test shifted.Ye*shifted.Pe≈s*I rtol=1e-10 atol=1e-10*abs(s)
    end
    @test_throws DomainError E.unified_earth!(
        E.EarthReturnWorkspace(geometry), state, integration; reference = 1.0)
    # Identical media erase the interface. The unbounded-medium direct field
    # provides an independent closed-form kernel for both ordered mixed cases.
    sigma=fill(0.1, 2)
    epsilon=fill(8.8541878128e-12, 2)
    mu=fill(4pi*1e-7, 2)
    state=(jω = s, Γ = zero(s), sigma, epsilon, mu,
        gamma_medium_squared = s .* mu .* (sigma .+ s .* epsilon))
    free=E.unified_earth!(E.EarthReturnWorkspace(geometry), state, integration)
    u=E.unified_earth_state(state, geometry)
    k=u.k[1]
    for p in 1:3, q in 1:3

        D=p==q ? geometry.radius[p] :
          hypot(geometry.horizontal[p]-geometry.horizontal[q], geometry.height[p]-geometry.height[q])
        trace=E.special_besselk(0, k*D)*exp(u.scaling[q])
        p==q || (trace*=u.A[p]*exp(u.scaling[p]))
        @test free.K[p, q]≈s*mu[1]/(2pi)*trace rtol=1e-8
        @test free.H[p, q]≈s/(2pi*u.sh[1])*trace rtol=1e-8
    end
    # Large electrical sizes require the spectral decay and analytic weight
    # in one exponential. Equal media give an independent one-wire K0 field.
    single=E.EarthReturnGeometry([0.0], [-1.0], [0.0425])
    s=2pi*1e12im
    sigma=fill(10.0, 2)
    state=(jω = s, Γ = zero(s), sigma, epsilon, mu,
        gamma_medium_squared = s .* mu .* (sigma .+ s .* epsilon))
    u=E.unified_earth_state(state, single)
    direct=E.special_besselkx(0, u.x[1])*exp(u.scaling[1]-u.x[1])
    K=s*mu[1]/(2pi)*direct
    H=s/(2pi*u.sh[1])*direct
    L=inv(u.A[1])-u.F[1]*K
    for method in (:quad, :trapz)
        integration=E.computation_options(E.SpectralIntegral, (
            method, options = (rtol = 1e-8,)))
        actual=E.unified_earth!(E.EarthReturnWorkspace(single), state, integration)
        @test actual.Ze[1, 1]≈K/L rtol=1e-8
        @test actual.Pe[1, 1]≈H/L rtol=1e-8
        @test actual.Ye[1, 1]≈s*L/H rtol=1e-8
    end
end

@testitem "Engine / public earth defaults preserve uncertainty and share matching preparations" tags=[:unit] begin
    using Measurements, LinearAlgebra
    const E=LineCableModels.Engine
    radius=measurement(0.0425, 0.0001)
    depth=measurement(1.0, 0.001)
    spacing=measurement(1.0, 0.001)
    rho=measurement(0.1, 0.001)
    metal=Material(kind = :conductor, rho = 1.7241e-8)
    design=build(CableDesign, "uncertain earth circles",
        Stack(Group(:phase, Region(:core, Disk(radius), metal))))
    system=build(
        LineCableSystem, fill(design, 2), [(zero(spacing), -depth), (spacing, -depth)];
        connections = [(phase = i,) for i in 1:2])
    problem=LineParametersProblem(
        system; earth_props = homogeneous(rho = rho, eps_r = 1.0, mu_r = 1.0),
        frequencies = [1e6])
    responses=map((:quad, :trapz)) do method
        definition=formula(:default; options = (integration = (
            method, options = (rtol = 1e-8,)),))
        selected=Formulation(earth_impedance = definition, earth_admittance = definition,
            options = (ideal_transposition = false,))
        execution=computation_options(LineCableModelsCoaxial, (trace = true,))
        T=eltype(problem)
        blueprints=E.CableBlueprint{T}[E.flatten(LineCableModelsCoaxial(), d, T)
                                       for d in problem.system.designs]
        workspace=E.LineParametersWorkspace(problem, selected, execution, blueprints)
        z=only(workspace.buffers.earth_numerical.earth_impedance.systems).response
        p=only(workspace.buffers.earth_numerical.earth_admittance.systems).response
        @test z===p
        result=@inferred E._solve!(workspace, selected)
        @test z.assemblies[]==length(problem.frequencies)
        @test eltype(result.Z.values)===Complex{Measurement{Float64}}
        trace=workspace.capture
        Ze=trace.Zg[:, :, 1]
        Pe=trace.Pg[:, :, 1]
        Ye=(2pi*1e6*im)*(Pe\Matrix{eltype(Pe)}(I, 2, 2))
        @test maximum(E.spectral_magnitude.(Ye*Pe-(2pi*1e6*im)*I))<1e-6
        for matrix in (Ze, Pe, Ye), parameter in (rho, radius, depth, spacing)

            @test any(matrix) do value
                !iszero(Measurements.derivative(real(value),
                    parameter)) ||
                    !iszero(Measurements.derivative(imag(value), parameter))
            end
        end
        (Ze, Pe, Ye)
    end
    for (quad, trapz) in zip(responses...)
        @test maximum(E.spectral_magnitude.(quad-trapz) ./ E.spectral_magnitude.(quad))<1e-6
    end
    @test !E.earth_state_equal(rho, measurement(0.1, 0.001))
    @test E.earth_state_equal(rho, 2rho-rho)
    # Check the propagated physical derivatives against independent central
    # perturbations of the material/geometry inputs to the complete closure.
    function perturbed(values)
        rho, r, h, d=values
        geometry=E.EarthReturnGeometry([0.0, d], [-h, -h], [r, r])
        s=2pi*1e6im
        sigma=[0.0, inv(rho)]
        epsilon=fill(8.8541878128e-12, 2)
        mu=fill(4pi*1e-7, 2)
        state=(jω = s, Γ = zero(s), sigma, epsilon, mu,
            gamma_medium_squared = s .* mu .* (sigma .+ s .* epsilon))
        integration=E.computation_options(E.SpectralIntegral, (
            method = :quad, options = (rtol = 1e-10,)))
        w=E.unified_earth!(E.EarthReturnWorkspace(geometry), state, integration)
        return w.Ze, w.Pe, w.Ye
    end
    parameters=(rho, radius, depth, spacing)
    for (index, parameter) in enumerate(parameters)
        values=Measurements.value.(collect(parameters))
        step=1e-4*values[index]
        plus=copy(values)
        minus=copy(values)
        plus[index]+=step
        minus[index]-=step
        for (actual, upper, lower) in zip(first(responses), perturbed(plus), perturbed(minus))
            for i in eachindex(actual)
                derivative=complex(Measurements.derivative(real(actual[i]), parameter),
                    Measurements.derivative(imag(actual[i]), parameter))
                finite_difference=(upper[i]-lower[i])/(2step)
                @test derivative≈finite_difference rtol=1e-4 atol=1e-10
            end
        end
    end
end

@testitem "Engine / complete earth scalar types and air-root limit" tags=[:unit] begin
    using LinearAlgebra
    const E=LineCableModels.Engine
    setprecision(BigFloat, 128) do
        for T in (Float32, Float64, BigFloat)
            geometry=E.EarthReturnGeometry(T[0, 1], T[1, -1], T[0.02, 0.03])
            s=complex(zero(T), T(2pi*50))
            sigma=T[0, 0.01]
            epsilon=T(8.8541878128e-12) .* T[1, 10]
            mu=fill(T(4pi*1e-7), 2)
            state=(jω = s, Γ = zero(s), sigma, epsilon, mu,
                gamma_medium_squared = s .* mu .* (sigma .+ s .* epsilon))
            methods=T===BigFloat ? (:quad, :trapz) : (:quad, :trapz, :cim)
            reference=nothing
            for method in methods
                tolerance=T===Float32 ? 2e-4 : 1e-6
                integration=E.computation_options(E.SpectralIntegral, (method,
                    options = (rtol = tolerance,)))
                workspace=E.unified_earth!(E.EarthReturnWorkspace(geometry), state, integration)
                @test eltype(workspace.Ze)===Complex{T}
                for matrix in (workspace.Ze, workspace.Pe, workspace.Ye)
                    @test all(isfinite, matrix)
                end
                if reference===nothing
                    reference=(copy(workspace.Ze), copy(workspace.Pe), copy(workspace.Ye))
                else
                    for (actual, expected, floor) in zip(
                        (workspace.Ze, workspace.Pe, workspace.Ye),
                        reference, (
                            1e-14, 1e-14, 1e-10))
                        for i in eachindex(actual), component in (real, imag)

                            @test isapprox(component(actual[i]),
                                component(expected[i]); rtol = 10tolerance,
                                atol = floor+64eps(T)*abs(expected[i]))
                        end
                    end
                end
            end
        end
        for z in (big"0.001"+big"0.002"*im, big"1.0"+big"2.0"*im, big"10.0"*im, big"300.0"+big"300.0"*im),
            order in (0, 1)

            for f in (E.special_besselix, E.special_besseljx, E.special_besselkx)
                @test f(order, z)≈f(order, ComplexF64(z)) rtol=1e-12
            end
        end
    end
    geometry=E.EarthReturnGeometry([0.0, 1.0], [1.0, -1.0], [0.02, 0.03])
    s=2pi*1000im
    sigma=[0.0, 0.01]
    epsilon=8.8541878128e-12 .* [1.0, 10.0]
    mu=fill(4pi*1e-7, 2)
    gamma2=s .* mu .* (sigma .+ s .* epsilon)
    Γ=sqrt(gamma2[1])
    gamma2[1]=Γ^2
    state=(jω = s, Γ, sigma, epsilon, mu, gamma_medium_squared = gamma2)
    quad=E.computation_options(E.SpectralIntegral, (
        method = :quad, options = (rtol = 1e-9,)))
    reference=E.unified_earth!(E.EarthReturnWorkspace(geometry), state, quad)
    for method in (:trapz, :cim)
        integration=E.computation_options(E.SpectralIntegral, (
            method, options = (rtol = 1e-6,)))
        actual=E.unified_earth!(E.EarthReturnWorkspace(geometry), state, integration)
        for (value, expected) in zip(
            (actual.Ze, actual.Pe, actual.Ye), (
                reference.Ze, reference.Pe, reference.Ye))
            @test value≈expected rtol=1e-5 atol=1e-10
        end
    end
end
