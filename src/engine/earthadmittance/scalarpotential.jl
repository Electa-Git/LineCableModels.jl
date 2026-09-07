# Shared numerical form of the scalar potentials printed by Kikuchi (5.1)–
# (5.4) and Mariscotti (14)–(15). This is not the mixed cable potential kernel.
function _scalar_potential_functor(identifier,formula,rho,epsilon,mu,s,Γ,segments)
    _check(rho,epsilon,mu)
    length(rho)==2 || throw(ArgumentError("this scalar potential requires exactly two half-spaces"))
    all(x->x>0 && !isnan(x),rho) && all(x->isfinite(x) && x>0,epsilon) &&
        all(x->isfinite(x) && x>0,mu) && isfinite(s) &&
        iszero(real(s)) && !iszero(s) || throw(DomainError((rho,epsilon,mu,s),
            "scalar potential requires positive material values and nonzero real frequency"))
    longitudinal=Γ===nothing ? formula.routes.Γ(s,mu[1],epsilon[1]) : (Γ=Γ,squared=Γ^2)
    basic=_homogeneous_functor(identifier,formula,rho,epsilon,mu,s,longitudinal.Γ,segments)
    state=basic.state
    bulk=map(1:2) do i
        s*state.mu[i]*(state.sigma[i]+s*state.epsilon[i])
    end
    squared=longitudinal.squared
    if Γ===nothing && isinf(rho[1]) && formula.routes.Γ==routes(identifier).Γ
        squared=-bulk[1]
    end
    all(isfinite,(longitudinal.Γ,squared)) ||
        throw(DomainError(longitudinal,"longitudinal input must be finite"))
    values=merge(state,(gamma_medium_squared=Tuple(bulk),gamma_squared=squared))
    return Functor{formula_id(formula),typeof(formula.routes),typeof(values)}(formula.routes,values)
end

function _scalar_potential_integral(integrand,state,height,source)
    T=typeof(state.tolerance); scale=max(height,one(T))
    bulk=state.gamma_medium_squared
    transverse=bulk .+ state.gamma_squared
    points=T[0,1]
    for a in transverse
        root=abs(spectral_root(a,state.jω))*scale
        isfinite(root) && root>0 && push!(points,root)
    end
    for medium in 1:2
        other=3-medium
        transition=abs(bulk[medium]/bulk[other]*
            spectral_root(transverse[other],state.jω))*scale
        isfinite(transition) && transition>0 && push!(points,transition)
    end
    push!(points,T(Inf)); sort!(unique!(points))
    guarded=t->begin
        lambda=t/scale
        lambda>sqrt(floatmax(T))/2 && return zero(state.jω)
        integrand(lambda)/scale
    end
    result=quadgk(guarded,points...;rtol=state.tolerance)[1]
    return _complex_result(state.jω,result)
end

function _scalar_potential_coefficient(state,source,d,D,height,lateral)
    all(isfinite,(d,D,height,lateral)) && d>0 && D>0 && height>0 &&
        lateral>=0 || throw(DomainError((d,D,height,lateral),"invalid scalar-potential geometry"))
    bulk=state.gamma_medium_squared
    transverse=bulk .+ state.gamma_squared
    all(iszero,transverse) && throw(DomainError(transverse,
        "the homogeneous zero-transverse-constant potential has no finite remote reference"))
    chi=spectral_root(transverse[source],state.jω)
    direct=iszero(chi) ? log(D/d) : bessel_difference(chi,d,D)
    integral=_scalar_potential_integral(state,height,source) do lambda
        q1=spectral_root(lambda^2+transverse[1],state.jω)
        q2=spectral_root(lambda^2+transverse[2],state.jω)
        qsource=source==1 ? q1 : q2
        bulk[source]*exp(-height*qsource)*cos(lateral*lambda)/(bulk[1]*q2+bulk[2]*q1)
    end
    kappa=state.sigma[source]+state.jω*state.epsilon[source]
    return _complex_result(state.jω,state.jω/(2*(one(d)*π)*kappa)*(direct+2integral))
end
