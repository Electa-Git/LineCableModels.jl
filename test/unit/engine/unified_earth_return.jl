# Closure identities establish current algebra and conventions; independent
# cylindrical controls below provide physical validation for equal media.
@testitem "Engine / unified retains five-layout FEM baseline agreement" tags=[:unit] setup=[TestFixtures] begin
    root=joinpath(pkgdir(LineCableModels), "test/fixtures/reference/three_bare_wires",
        "capture-20260917T102438-Nt9nic")
    selected=Formulation(
        earth_impedance = formula(:unified), earth_admittance = formula(:unified),
        options = (
            reduce_bundle = false, kron_reduction = false, ideal_transposition = false))
    for layout in (:all_air, :all_earth, :air_1, :air_2, :air_3)
        directory=joinpath(root, "three_bare_wires_$layout")
        # The original problem JSON is historical input with problem-owned Γ.
        # Rebuild the same physical fixture; the retained numerical CSVs stay fixed.
        problem=TestFixtures.three_bare_wires_problem(
            heights=TestFixtures.three_bare_wires_layouts[layout])
        @test problem.frequencies==10.0 .^ (-1:7)
        result=compute(problem, selected; options = (verbosity = (default = 0,),))
        for (name, selector) in (("Z", Z), ("Y", Y))
            rows=map(
                backend->split.(readlines(joinpath(directory, backend, "$name.csv"))[2:end], ','),
                ("unified", "fem"))
            @test all(a[1:3]==b[1:3] for (a, b) in zip(rows...))
            old,
            fem=map(rs->[complex(parse(Float64, r[4]), parse(Float64, r[5])) for r in rs], rows)
            values=observe(result, selector)
            actual=[values[i, j, k] for k in 1:9 for i in 1:3 for j in 1:3]
            @test length(actual)==length(old)==length(fem)==81
            # Preserve the observed worst discrepancy; this does not claim
            # exact FEM agreement or impose a solver-side acceptance policy.
            roundoff=64eps(Float64)*max(maximum(abs, old), maximum(abs, fem))
            @test maximum(abs.(actual-fem))<=maximum(abs.(old-fem))+roundoff
            # With Γ=0 the reference correction does not change Z.
            name=="Z" && (@test actual≈old rtol=1e-11 atol=1e-13)
        end
    end
end

@testitem "Engine / full earth / current closure across layouts" tags=[:unit] setup=[UnifiedFormulaFixtures] begin
    using LinearAlgebra
    const E=LineCableModels.Engine
    layouts=([(0.0, -1.0)], [(0.0, 1.0), (0.4, 1.5)],
        [(0.0, -1.0), (0.4, -1.5)], [(0.0, 1.0), (0.4, -1.0), (0.9, -1.5)])
    for positions in layouts, f in (50.0, 10000.0)

        n=length(positions)
        geometry=(horizontal = first.(positions), height = last.(positions),
            radius = [0.005, 0.007, 0.009][1:n])
        for gamma in (n==3 ? (0im, 1e-4*(1+im)) : (0im,))
            s=2pi*im*f;
            sigma=[0.0, 0.01];
            epsilon=8.8541878128e-12 .* [1, 10];
            mu=fill(4pi*1e-7, 2)
            state=(jω = s, Γ = gamma, sigma, epsilon, mu)
            @testset "n=$n / $f Hz / gamma=$gamma / $method" for method in (:quad,)
                limits=(; maxevals = 10^6)
                controls=E.formulation_options(
                    E.SpectralIntegral, (method, options = merge((rtol = 1e-9,), limits)))
                w=E.EarthImpedance._unified_current!(
                    UnifiedFormulaFixtures.buffers(geometry),
                    geometry,
                    state,
                    controls)
                # Normalized backward errors use independent operand norms.
                for (lhs, rhs, scale) in ((w.Pe*w.L, w.H, norm(w.Pe)*norm(w.L)+norm(w.H)),
                    (w.Ze*w.L, w.K+gamma^2*w.H/s,
                        norm(w.Ze)*norm(w.L)+norm(w.K)+abs(gamma^2/s)*norm(w.H)),
                    ((s*inv(w.Pe))*w.H, s*w.L,
                        norm((s*inv(w.Pe)))*norm(w.H)+abs(s)*norm(w.L)),
                    ((s*inv(w.Pe))*w.Pe, s*I,
                        norm((s*inv(w.Pe)))*norm(w.Pe)+abs(s)*sqrt(n)))
                    @test norm(lhs-rhs)/scale <= 1e-10
                end
            end
        end
    end
end

@testitem "Engine / manuscript geometry and reference identities" tags=[:unit] setup=[UnifiedFormulaFixtures] begin
    using LinearAlgebra
    const E=LineCableModels.Engine
    crossing=[E.EarthPair(1, 1, (-0.01, -0.01), 0.0, (2, 2); radius = 0.02)]
    overlapping=[E.EarthPair(1, 1, (-1.0, -1.0), 0.0, (2, 2); radius = 0.02),
        E.EarthPair(2, 2, (-1.0, -1.0), 0.0, (2, 2); radius = 0.02),
        E.EarthPair(1, 2, (-1.0, -1.0), 0.01, (2, 2))]
    @test_throws DomainError E.EarthImpedance._unified_geometry(crossing)
    @test_throws DomainError E.EarthImpedance._unified_geometry(overlapping)
    geometry=(horizontal = [0.0, 1.0, 2.0], height = [1.2, -0.9, -1.4],
        radius = [
            0.01, 0.025, 0.04])
    s=complex(0.0, 2pi*1e4)
    Γ=complex(1e-4, 2e-4)
    sigma=[0.0, 0.1]
    epsilon=8.8541878128e-12 .* [1.0, 8.0]
    mu=4pi*1e-7 .* [1.0, 3.0]
    state=(jω = s, Γ, sigma, epsilon, mu)
    integration=E.formulation_options(E.SpectralIntegral, (
        method = :quad, options = (rtol = 1e-9,)))
    @test_throws ArgumentError E.EarthImpedance.Formula(:unified; parameters = (reference = :deep,))
    @test_throws ArgumentError E.EarthAdmittance.Formula(:unified; parameters = (reference = :interface,))
    # Identical media erase the interface. The unbounded-medium direct field
    # provides a closed-form consistency check using production special functions.
    sigma=fill(0.1, 2)
    epsilon=fill(8.8541878128e-12, 2)
    mu=fill(4pi*1e-7, 2)
    state=(jω = s, Γ = zero(s), sigma, epsilon, mu)
    free=E.EarthImpedance._unified_current!(
        UnifiedFormulaFixtures.buffers(geometry),
        geometry,
        state,
        integration)
    u=E.EarthImpedance._unified_state!(
        UnifiedFormulaFixtures.buffers(geometry).unified.current, state, geometry)
    k=u.k[1]
    for p in 1:3, q in 1:3

        D=p==q ? geometry.radius[p] :
          hypot(geometry.horizontal[p]-geometry.horizontal[q], geometry.height[p]-geometry.height[q])
        trace=E.special_besselk(0, k*D)*exp(u.scaling[q])
        p==q||(trace*=u.A[p]*exp(u.scaling[p]))
        @test free.K[p, q]≈s*mu[1]/(2pi)*trace rtol=1e-8
        geometry.height[p]<0&&(@test free.H[p, q]≈s/(2pi*u.sh[1])*trace rtol=1e-8)
    end
    # Large electrical sizes require the spectral decay and analytic weight
    # in one exponential. This is a consistency check with the production K0 path.
    single=(horizontal = [0.0], height = [-1.0], radius = [0.0425])
    s=2pi*1e12im
    sigma=fill(10.0, 2)
    state=(jω = s, Γ = zero(s), sigma, epsilon, mu)
    u=E.EarthImpedance._unified_state!(UnifiedFormulaFixtures.buffers(single).unified.current, state, single)
    direct=E.special_besselkx(0, u.x[1])*exp(u.scaling[1]-u.x[1])
    K=s*mu[1]/(2pi)*direct
    H=s/(2pi*u.sh[1])*direct
    L=inv(u.A[1])-u.F[1]*K
    for method in (:quad,)
        local integration=E.formulation_options(E.SpectralIntegral, (
            method, options = (rtol = 1e-8,)))
        actual=E.EarthImpedance._unified_current!(
            UnifiedFormulaFixtures.buffers(single),
            single,
            state,
            integration)
        @test actual.Ze[1, 1]≈K/L rtol=1e-8
        @test actual.Pe[1, 1]≈H/L rtol=1e-8
        @test (s * inv(actual.Pe))[1, 1]≈s*L/H rtol=1e-8
    end
end

@testitem "Engine / raw integral trace and warnings do not reject physical computation" tags=[:unit] setup=[TestFixtures] begin
    using Logging
    const E=LineCableModels.Engine
    problem=TestFixtures.three_bare_wires_problem(frequencies = [50.0])
    controls=(integration = (method = :quad, options = (rtol = 1e-14, maxevals = 15)),)
    definition=formula(:unified; options = controls)
    selected=Formulation(earth_impedance = definition, earth_admittance = definition,
        options = (
            reduce_bundle = false, kron_reduction = false, ideal_transposition = false))
    # Public compute installs its normal console logger; observe its real
    # stderr output rather than replacing the outer task logger.
    result,
    warning=mktemp() do path, io
        computed=redirect_stderr(io) do
            compute(problem, selected; options = (trace = true, verbosity = (default = 0,)))
        end
        flush(io);
        seekstart(io)
        (computed, read(io, String))
    end
    @test occursin("QuadGK returned an estimated error", warning)
    @test occursin("estimated_error", warning)&&occursin("context", warning)
    records=details(result).data.trace.integrals
    @test length(records)==18 # One shared assembly, not one per owner or entry request.
    @test all(r->isfinite(r.value)&&isfinite(r.estimated_error)&&r.estimated_error>=0, records)
    @test Set(r.context.formula for r in records)==Set([:unified])
    @test Set(r.context.term for r in records)==Set([:Z, :air_voltage])
    @test all(r->r.context.frequency==50.0, records)
    @test Set((r.context.receiver, r.context.source) for r in records) ==
          Set((p, q) for p in 1:3 for q in 1:3)
    @test any(r->r.estimated_error>1e-14*abs(r.value), records)
    @test all(isfinite, result.Z.values)&&all(isfinite, result.Y.values)
    execution=computation_options(LineCableModelsCoaxial,
        ComputationOptions(trace = true, verbosity = (default = 0,)))
    blueprints=[E.flatten(LineCableModelsCoaxial(), d, Float64)
                for d in problem.system.designs]
    w=E.LineParametersWorkspace(problem, Formulation(), execution, blueprints)
    E._solve!(w, Formulation());
    E._solve!(w, Formulation())
    @test length(w.capture.integrals)==18
    quiet=compute(problem; options = (trace = false, verbosity = (default = 0,)))
    @test !haskey(details(quiet).data, :trace)
end

@testitem "Engine / public earth defaults preserve uncertainty and share matching preparations" tags=[:unit] setup=[UnifiedFormulaFixtures] begin
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
    responses=map((:quad,)) do method
        definition=formula(:default; options = (integration = (
            method, options = (rtol = 1e-8,)),))
        selected=Formulation(earth_impedance = definition, earth_admittance = definition,
            options = (ideal_transposition = false,))
        execution=computation_options(LineCableModelsCoaxial, ComputationOptions((trace = true,)))
        T=eltype(problem)
        blueprints=E.CableBlueprint{T}[E.flatten(LineCableModelsCoaxial(), d, T)
                                       for d in problem.system.designs]
        workspace=E.LineParametersWorkspace(problem, selected, execution, blueprints)
        @test only(workspace.invariants.earth_bindings.earth_impedance.cases).partner == 1
        @test (@inferred E._solve!(workspace, selected)) === workspace
        @test all(isfinite, workspace.buffers.unified.Ze)
        @test !isempty(workspace.capture.integrals)
        @test eltype(workspace.buffers.Zout)===Complex{Measurement{Float64}}
        trace=workspace.capture
        Ze=trace.Zg[:, :, 1]
        Pe=trace.Pg[:, :, 1]
        Ye=(2pi*1e6*im)*(Pe\Matrix{eltype(Pe)}(I, 2, 2))
        @test maximum(E.numerical_magnitude.(Ye*Pe-(2pi*1e6*im)*I))<1e-6
        for matrix in (Ze, Pe, Ye), parameter in (rho, radius, depth, spacing)

            @test any(matrix) do value
                !iszero(Measurements.derivative(real(value),
                    parameter)) ||
                    !iszero(Measurements.derivative(imag(value), parameter))
            end
        end
        (Ze, Pe, Ye)
    end
    @test !E.same_physical_state(rho, measurement(0.1, 0.001))
    @test E.same_physical_state(rho, 2rho-rho)
    # Check the propagated physical derivatives against independent central
    # perturbations of the material/geometry inputs to the complete closure.
    function perturbed(values)
        rho, r, h, d=values
        geometry=(horizontal = [0.0, d], height = [-h, -h], radius = [r, r])
        s=2pi*1e6im
        sigma=[0.0, inv(rho)]
        epsilon=fill(8.8541878128e-12, 2)
        mu=fill(4pi*1e-7, 2)
        state=(jω = s, Γ = zero(s), sigma, epsilon, mu)
        integration=E.formulation_options(E.SpectralIntegral, (
            method = :quad, options = (rtol = 1e-10,)))
        w=E.EarthImpedance._unified_current!(
            UnifiedFormulaFixtures.buffers(geometry),
            geometry,
            state,
            integration)
        return w.Ze, w.Pe, (s*inv(w.Pe))
    end
    parameters=(rho, radius, depth, spacing)
    for (index, parameter) in enumerate(parameters)
        values=Measurements.value.(collect(parameters))
        step=1e-4*values[index]
        plus=copy(values)
        minus=copy(values)
        plus[index]+=step
        minus[index]-=step
        for (actual, upper,
            lower) in
        zip(first(responses), perturbed(plus), perturbed(minus))
            for i in eachindex(actual)
                derivative=complex(Measurements.derivative(real(actual[i]), parameter),
                    Measurements.derivative(imag(actual[i]), parameter))
                finite_difference=(upper[i]-lower[i])/(2step)
                @test derivative≈finite_difference rtol=1e-4 atol=1e-10
            end
        end
    end
end

@testitem "Engine / complete earth scalar types and air-root limit" tags=[:unit] setup=[UnifiedFormulaFixtures] begin
    using LinearAlgebra
    const E=LineCableModels.Engine
    setprecision(BigFloat, 128) do
        for T in (Float32, Float64, BigFloat)
            geometry=(horizontal = T[0, 1], height = T[1, -1], radius = T[0.02, 0.03])
            s=complex(zero(T), T(2pi*50))
            sigma=T[0, 0.01]
            epsilon=T(8.8541878128e-12) .* T[1, 10]
            mu=fill(T(4pi*1e-7), 2)
            state=(jω = s, Γ = zero(s), sigma, epsilon, mu)
            for method in (:quad,)
                tolerance=T===Float32 ? 2e-4 : 1e-6
                integration=E.formulation_options(E.SpectralIntegral, (method,
                    options = (rtol = tolerance,)))
                workspace=E.EarthImpedance._unified_current!(
                    UnifiedFormulaFixtures.buffers(geometry),
                    geometry,
                    state,
                    integration)
                @test eltype(workspace.Ze)===Complex{T}
                for matrix in (workspace.Ze, workspace.Pe, (s*inv(workspace.Pe)))
                    @test all(isfinite, matrix)
                end
            end
        end
        for z in (big"0.001"+big"0.002"*im, big"1.0"+big"2.0"*im,
                big"10.0"*im, big"300.0"+big"300.0"*im),
            order in (0, 1)

            for f in (E.special_besselix, E.special_besseljx, E.special_besselkx)
                @test f(order, z)≈f(order, ComplexF64(z)) rtol=1e-12
            end
        end
    end
    geometry=(horizontal = [0.0, 1.0], height = [1.0, -1.0], radius = [0.02, 0.03])
    s=2pi*1000im
    sigma=[0.0, 0.01]
    epsilon=8.8541878128e-12 .* [1.0, 10.0]
    mu=fill(4pi*1e-7, 2)
    gamma2=s .* mu .* (sigma .+ s .* epsilon)
    Γ=sqrt(gamma2[1])
    gamma2[1]=Γ^2
    state=(jω = s, Γ, sigma, epsilon, mu)
    quad=E.formulation_options(E.SpectralIntegral, (
        method = :quad, options = (rtol = 1e-9,)))
    reference=E.EarthImpedance._unified_current!(
        UnifiedFormulaFixtures.buffers(geometry),
        geometry,
        state,
        quad)
    for method in (:quad,)
        integration=E.formulation_options(E.SpectralIntegral, (
            method, options = (rtol = 1e-6,)))
        actual=E.EarthImpedance._unified_current!(
            UnifiedFormulaFixtures.buffers(geometry),
            geometry,
            state,
            integration)
        for (value,
            expected) in zip(
            (actual.Ze, actual.Pe, (s*inv(actual.Pe))), (
                reference.Ze, reference.Pe, (s*inv(reference.Pe))))
            @test value≈expected rtol=1e-5 atol=1e-10
        end
    end
end

@testitem "Engine / full earth / independent equal-medium cylindrical control" tags=[:unit] setup=[UnifiedFormulaFixtures] begin
    using LinearAlgebra, QuadGK
    include(joinpath(pkgdir(LineCableModels), "test/support/radial_control.jl"))
    const E=LineCableModels.Engine
    # I0/I1 use an independently coded Frobenius recurrence; K0 uses its
    # decaying integral, not the engine's selected spectral integrator.
    function reference(positions, radii, f)
        s=2big(pi)*im*BigFloat(f);
        mu=4big(pi)*big"1e-7"
        kappa=big".01"+s*3big"8.8541878128e-12";
        k=sqrt(s*mu*kappa)
        n=length(radii);
        A=zeros(Complex{BigFloat}, n);
        F=similar(A)
        eA=zeros(BigFloat, n);
        eF=similar(eA)
        for p in 1:n
            x=k*radii[p];
            v=RadialControl.fundamental(x)
            A[p]=v.u;
            F[p]=2big(pi)*kappa*radii[p]^2*v.du/(x*v.u)
            eA[p]=v.bound
            eA[p]<abs(A[p])/2||error("inconclusive cylindrical reference denominator")
            eF[p]=abs(2big(pi)*kappa*radii[p]^2/x)*v.bound*
                  (1+abs(v.du/A[p]))/(abs(A[p])-eA[p])
        end
        T=zeros(Complex{BigFloat}, n, n);
        eT=zeros(BigFloat, n, n)
        for p in 1:n, q in 1:n

            distance=p==q ? radii[p] : hypot((positions[p] .- positions[q])...)
            z=k*distance;
            limit=log(160/real(z))
            value,
            error=quadgk(
                t->exp(-z*cosh(t)), big"0", limit; rtol = big"1e-12", maxevals = 10^6)
            tail=2exp(-real(z)*exp(limit)/2)/(real(z)*exp(limit))
            T[p, q]=(p==q ? one(k) : A[p])*value
            eT[p, q]=p==q ? error+tail :
                     abs(A[p])*(error+tail)+eA[p]*(abs(value)+error+tail)
        end
        K=s*mu*T/(2big(pi));
        H=s*T/(2big(pi)*kappa)
        L=diagm(inv.(A))-diagm(F)*K
        Ze=K/L;
        Pe=H/L;
        Ye=s*inv(Pe)
        # Normwise perturbation bounds for right solves X*L=B. For
        # eta=||L^-1||*deltaL<1, deltaX <=
        # (deltaB+||X||*deltaL)*||L^-1||/(1-eta).
        # Absolute input bounds retain units and apply to either component.
        norminf(x)=opnorm(x, Inf)
        roundoff(x)=128eps(BigFloat)*norminf(x)
        eK=abs(s*mu/(2big(pi)))*norminf(eT)+roundoff(K)
        eH=abs(s/(2big(pi)*kappa))*norminf(eT)+roundoff(H)
        eL=maximum(eA ./ (abs.(A) .* (abs.(A) .- eA)))+
           maximum(abs, F)*eK+maximum(eF)*(norminf(K)+eK)+roundoff(L)
        inverseL=norminf(inv(L));
        eta=inverseL*eL
        eta<1||error("inconclusive cylindrical closure conditioning")
        eZe=(eK+norminf(Ze)*eL)*inverseL/(1-eta)+roundoff(Ze)
        ePe=(eH+norminf(Pe)*eL)*inverseL/(1-eta)+roundoff(Pe)
        inverseP=norminf(inv(Pe));
        etaP=inverseP*ePe
        etaP<1||error("inconclusive cylindrical admittance conditioning")
        eYe=abs(s)*inverseP^2*ePe/(1-etaP)+roundoff(Ye)
        bounds=(K = eK, H = eH, L = eL, Ze = eZe, Pe = ePe, Ye = eYe)
        return (; K, H, L, Ze, Pe, Ye, k, bounds)
    end
    layouts=([(0.0, -1.0)], [(0.0, 1.0), (0.4, 1.5)], [(0.0, -1.0), (0.4, -1.5)],
        [(0.0, 1.0), (0.4, -1.0), (0.9, -1.5)])
    for positions in layouts, f in (50.0, 10000.0)

        n=length(positions);
        radii=[0.005, 0.007, 0.009][1:n]
        refs=[setprecision(BigFloat, bits) do
                  reference(map(p->BigFloat.(p), positions), BigFloat.(radii), f)
              end
              for bits in (128, 256, 512)]
        expected=last(refs)
        geometry=(horizontal = first.(positions), height = last.(positions), radius = radii)
        s=2pi*im*f;
        sigma=fill(0.01, 2);
        epsilon=fill(3*8.8541878128e-12, 2);
        mu=fill(4pi*1e-7, 2)
        state=(jω = s, Γ = zero(s), sigma, epsilon, mu)
        controls=E.formulation_options(
            E.SpectralIntegral, (method = :quad, options = (rtol = 1e-10, maxevals = 10^6)))
        w=E.EarthImpedance._unified_current!(
            UnifiedFormulaFixtures.buffers(geometry),
            geometry,
            state,
            controls)
        column_scale=reshape(exp.(abs.(real.(expected.k .* radii))), 1, :)
        # Equal-medium impedance/current maps erase the interface. Air
        # voltages deliberately retain an interface reference, so only buried
        # layouts have the unbounded-medium P/Y used by this independent control.
        for quantity in
        (all(p->last(p)<0, positions) ? (:K, :H, :L, :Ze, :Pe, :Ye) : (:K, :L, :Ze))
            target=getproperty(expected, quantity)
            actual=quantity in (:K, :H, :L) ? getproperty(w, quantity) ./ column_scale :
                   quantity===:Ye ? s*inv(w.Pe) : getproperty(w, quantity)
            for index in eachindex(target), component in (real, imag)

                qstar=component(target[index]);
                budget=1e-5*abs(qstar)
                u=getproperty(expected.bounds,
                    quantity)+
                abs(component(target[index]-getproperty(first(refs), quantity)[index]))
                @test u<=budget/4
                @test abs(component(actual[index])-qstar)+u<=budget
            end
        end
    end
end

@testitem "Engine / combined voltage paths equal the manuscript endpoint terms" tags=[:unit] setup=[UnifiedFormulaFixtures] begin
    const E=LineCableModels.Engine
    geometry=(horizontal = [0.0, 1.0], height = [1.2, -0.9], radius = [0.01, 0.025])
    s=complex(0.0, 2pi*1e4)
    sigma=[0.0, 0.1]
    epsilon=8.8541878128e-12 .* [1.0, 8.0]
    mu=4pi*1e-7 .* [1.0, 3.0]
    u=E.EarthImpedance._unified_state!(
        UnifiedFormulaFixtures.buffers(geometry).unified.current,
        (jω = s, Γ = 1e-4+2e-4im, sigma, epsilon, mu),
        geometry)
    for P in (1,), Q in (1, 2)

        hp=abs(geometry.height[P])
        hq=abs(geometry.height[Q])
        radius=geometry.radius[P]
        padding=hq/2
        sq=u.scaling[Q]
        g=(; hp, hq, radius, padding, logscale = sq,
            i0minus = E.EarthImpedance.bessel_i0m1(u.x[P]))
        combined=E.EarthImpedance.AirVoltageSpectrum{Q, typeof(u), typeof(g)}(u, g)
        voltage=E.EarthImpedance.earth_spectrum(Val(:voltage), Val(P), Val(Q), u, hp, hq,
            u.scaling[P]+sq)
        for λ in (0.001, 0.1, 3.0, 40.0) .* exp(0.03im)
            j0=E.SpecialFunctions.besselj(0, radius*λ)
            original=u.A[P]*voltage(λ)*exp(-(hp+hq)*λ)
            a0=sqrt(λ^2+u.k2[1]);
            ag=sqrt(λ^2+u.k2[2])
            aq=Q==1 ? a0 : ag
            original-=ag*j0*exp(sq-hq*aq)/(a0*(u.sh[2]*a0+u.sh[1]*ag))
            actual=combined(λ)*exp(-(hq-padding)*λ)
            @test actual≈original rtol=2e-10 atol=1e-20
        end
        # At a_receiver=0, I0(κ_receiver*r)=J0(r*λ). The removable
        # quotient tends to -h_receiver*J0, including the spectral padding.
        roots=ntuple(m->m==P ? complex(-1.0) : u.k2[m], 2)
        branch=merge(u, (; k2 = roots))
        kernel=E.EarthImpedance.AirVoltageSpectrum{Q, typeof(branch), typeof(g)}(branch, g)
        a=map(k2->E.EarthImpedance.outgoing_root(1+k2), roots)
        other=a[3 - P]
        j0=E.SpecialFunctions.besselj(0, radius)
        expected=-hp*other*j0*exp(sq-hq*a[Q])/(u.sh[2]*a[1]+u.sh[1]*a[2])
        @test kernel(1.0)*exp(-(hq-padding))≈expected rtol=1e-12
    end
end

@testitem "Engine / prescribed longitudinal constants keep spectral contours on their branch" tags=[:unit] setup=[UnifiedFormulaFixtures] begin
    const E=LineCableModels.Engine
    geometry=(horizontal = [0.0, 1.0], height = [1.0, -1.0], radius = [0.02, 0.03])
    s=complex(0.0, 2pi*1e4)
    sigma=[0.0, 0.1]
    epsilon=8.8541878128e-12 .* [1.0, 8.0]
    mu=4pi*1e-7 .* [1.0, 3.0]
    controls=E.formulation_options(E.SpectralIntegral, (
        method = :quad, options = (rtol = 1e-9,))).options
    for Γ in (0.0+0.0im, 1e-4+2e-4im, 3e-4+1e-4im)
        state=E.EarthImpedance._unified_state!(
            UnifiedFormulaFixtures.buffers(geometry).unified.current,
            (jω = s, Γ, sigma, epsilon, mu),
            geometry)
        angle=E.EarthImpedance.earth_contour_angle(state, pi/6)
        @test 0<angle<=pi/6
        iszero(Γ)||@test angle<pi/6
        for P in (1, 2), Q in (1, 2), kind in (:Z, :phi, :voltage)
            # At Γ=0 the separate air-voltage term has a branch-point
            # singularity; production combines its cancelling endpoints.
            # Here the added zero-Γ check concerns the scalar coordinate.
            iszero(Γ)&&kind!==:phi&&continue
            kernel=E.EarthImpedance.earth_spectrum(
                Val(kind), Val(P), Val(Q), state, 1.0, 1.0, 0.0)
            integral=E.SpectralIntegral(lambda->kernel(lambda)*exp(-2lambda)*cos(lambda))
            points=E.EarthImpedance._unified_points!(
                UnifiedFormulaFixtures.buffers(geometry).unified, state, 2.0, 1.0, 0.0, 0.0)
            push!(points, 1.0)
            real_axis, _=E.integrate(Val(:quad), integral, controls; points)
            rotated=E.EarthImpedance.earth_spectral_term(
                Val(kind), Val(P), Val(Q), state, 1.0, 1.0, 1.0, 0.0, 0.0,
                Val(:quad), controls, UnifiedFormulaFixtures.buffers(geometry))
            @test rotated≈real_axis rtol=1e-7
        end
    end
end
