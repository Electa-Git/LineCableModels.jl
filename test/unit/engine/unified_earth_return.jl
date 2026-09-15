# Closure identities establish current algebra and conventions; unlike-media
# cross-integrator agreement is consistency evidence, not physical validation.
@testitem "Engine / full earth / current closure across layouts" tags=[:unit] begin
    using LinearAlgebra
    const E=LineCableModels.Engine
    layouts=([(0.0,-1.0)],[(0.0,1.0),(0.4,1.5)],
        [(0.0,-1.0),(0.4,-1.5)],[(0.0,1.0),(0.4,-1.0),(0.9,-1.5)])
    for positions in layouts, f in (50.0,10000.0)
        n=length(positions)
        geometry=E.EarthReturnGeometry(first.(positions),last.(positions),[.005,.007,.009][1:n])
        for gamma in (n==3 ? (0im,1e-4*(1+im)) : (0im,))
            s=2pi*im*f; sigma=[0.0,.01]; epsilon=8.8541878128e-12.*[1,10]; mu=fill(4pi*1e-7,2)
            state=(jω=s,Γ=gamma,sigma,epsilon,mu,gamma_medium_squared=s.*mu.*(sigma.+s.*epsilon))
            @testset "n=$n / $f Hz / gamma=$gamma / $method" for method in (:quad,:trapz,:cim)
                # Closure algebra does not require a fixed construction sample
                # cap. Adaptive construction and explicit exhaustion have
                # separate contracts under the spectral owner.
                limits=method===:trapz ? (;max_refinements=14) : (;maxevals=10^6)
                controls=E.computation_options(E.SpectralIntegral,(method,options=merge((rtol=1e-9,),limits)))
                w=E.unified_earth!(E.EarthReturnWorkspace(geometry),state,controls)
                # Normalized backward errors use independent operand norms.
                for (lhs,rhs,scale) in ((w.Pe*w.L,w.H,norm(w.Pe)*norm(w.L)+norm(w.H)),
                    (w.Ze*w.L,w.K+gamma^2*w.H/s,norm(w.Ze)*norm(w.L)+norm(w.K)+abs(gamma^2/s)*norm(w.H)),
                    (w.Ye*w.H,s*w.L,norm(w.Ye)*norm(w.H)+abs(s)*norm(w.L)),
                    (w.Ye*w.Pe,s*I,norm(w.Ye)*norm(w.Pe)+abs(s)*sqrt(n)))
                    @test norm(lhs-rhs)/scale <= 1e-10
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
    # provides a closed-form consistency check using production special functions.
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
    # in one exponential. This is a consistency check with the production K0 path.
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

@testitem "Engine / full earth / independent equal-medium cylindrical control" tags=[:unit] begin
    using LinearAlgebra, QuadGK
    include(joinpath(pkgdir(LineCableModels),"test/support/radial_control.jl"))
    const E=LineCableModels.Engine
    # I0/I1 use an independently coded Frobenius recurrence; K0 uses its
    # decaying integral, not the engine's selected spectral integrator.
    function reference(positions,radii,f)
        s=2big(pi)*im*BigFloat(f); mu=4big(pi)*big"1e-7"
        kappa=big".01"+s*3big"8.8541878128e-12"; k=sqrt(s*mu*kappa)
        n=length(radii); A=zeros(Complex{BigFloat},n); F=similar(A)
        eA=zeros(BigFloat,n);eF=similar(eA)
        for p in 1:n
            x=k*radii[p]; v=RadialControl.fundamental(x)
            A[p]=v.u; F[p]=2big(pi)*kappa*radii[p]^2*v.du/(x*v.u)
            eA[p]=v.bound
            eA[p]<abs(A[p])/2 || error("inconclusive cylindrical reference denominator")
            eF[p]=abs(2big(pi)*kappa*radii[p]^2/x)*v.bound*
                (1+abs(v.du/A[p]))/(abs(A[p])-eA[p])
        end
        T=zeros(Complex{BigFloat},n,n); eT=zeros(BigFloat,n,n)
        for p in 1:n,q in 1:n
            distance=p==q ? radii[p] : hypot((positions[p].-positions[q])...)
            z=k*distance; limit=log(160/real(z))
            value,error=quadgk(t->exp(-z*cosh(t)),big"0",limit;rtol=big"1e-12",maxevals=10^6)
            tail=2exp(-real(z)*exp(limit)/2)/(real(z)*exp(limit))
            T[p,q]=(p==q ? one(k) : A[p])*value
            eT[p,q]=p==q ? error+tail :
                abs(A[p])*(error+tail)+eA[p]*(abs(value)+error+tail)
        end
        K=s*mu*T/(2big(pi)); H=s*T/(2big(pi)*kappa)
        L=diagm(inv.(A))-diagm(F)*K
        Ze=K/L; Pe=H/L; Ye=s*inv(Pe)
        # Normwise perturbation bounds for right solves X*L=B. For
        # eta=||L^-1||*deltaL<1, deltaX <=
        # (deltaB+||X||*deltaL)*||L^-1||/(1-eta).
        # Absolute input bounds retain units and apply to either component.
        norminf(x)=opnorm(x,Inf)
        roundoff(x)=128eps(BigFloat)*norminf(x)
        eK=abs(s*mu/(2big(pi)))*norminf(eT)+roundoff(K)
        eH=abs(s/(2big(pi)*kappa))*norminf(eT)+roundoff(H)
        eL=maximum(eA./(abs.(A).*(abs.(A).-eA)))+
            maximum(abs,F)*eK+maximum(eF)*(norminf(K)+eK)+roundoff(L)
        inverseL=norminf(inv(L));eta=inverseL*eL
        eta<1 || error("inconclusive cylindrical closure conditioning")
        eZe=(eK+norminf(Ze)*eL)*inverseL/(1-eta)+roundoff(Ze)
        ePe=(eH+norminf(Pe)*eL)*inverseL/(1-eta)+roundoff(Pe)
        inverseP=norminf(inv(Pe));etaP=inverseP*ePe
        etaP<1 || error("inconclusive cylindrical admittance conditioning")
        eYe=abs(s)*inverseP^2*ePe/(1-etaP)+roundoff(Ye)
        bounds=(K=eK,H=eH,L=eL,Ze=eZe,Pe=ePe,Ye=eYe)
        return (;K,H,L,Ze,Pe,Ye,k,bounds)
    end
    layouts=([(0.0,-1.0)],[(0.0,1.0),(.4,1.5)],[(0.0,-1.0),(.4,-1.5)],
        [(0.0,1.0),(.4,-1.0),(.9,-1.5)])
    for positions in layouts,f in (50.0,10000.0)
        n=length(positions);radii=[.005,.007,.009][1:n]
        refs=[setprecision(BigFloat,bits) do
            reference(map(p->BigFloat.(p),positions),BigFloat.(radii),f)
        end for bits in (128,256,512)]
        expected=last(refs)
        geometry=E.EarthReturnGeometry(first.(positions),last.(positions),radii)
        s=2pi*im*f; sigma=fill(.01,2);epsilon=fill(3*8.8541878128e-12,2);mu=fill(4pi*1e-7,2)
        state=(jω=s,Γ=zero(s),sigma,epsilon,mu,gamma_medium_squared=s.*mu.*(sigma.+s.*epsilon))
        controls=E.computation_options(E.SpectralIntegral,(method=:quad,options=(rtol=1e-10,maxevals=10^6)))
        w=E.unified_earth!(E.EarthReturnWorkspace(geometry),state,controls)
        column_scale=reshape(exp.(abs.(real.(expected.k.*radii))),1,:)
        for quantity in (:K,:H,:L,:Ze,:Pe,:Ye)
            target=getproperty(expected,quantity)
            actual=quantity in (:K,:H,:L) ? getproperty(w,quantity)./column_scale : getproperty(w,quantity)
            for index in eachindex(target),component in (real,imag)
                qstar=component(target[index]);budget=1e-5*abs(qstar)
                u=getproperty(expected.bounds,quantity)+
                    abs(component(target[index]-getproperty(first(refs),quantity)[index]))
                @test u<=budget/4
                @test abs(component(actual[index])-qstar)+u<=budget
            end
        end
    end
end
