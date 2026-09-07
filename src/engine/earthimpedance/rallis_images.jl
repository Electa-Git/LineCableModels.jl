# Rallis2013, Appendix A and (4.5), (4.11), (4.20).
# The finite fit is intentionally computed in ComplexF64: its empirical error
# dominates the coefficient roundoff. Final arithmetic retains the caller type.
function _rallis_gpof(function_value,span,origin,samples,terms)
    t=range(0.0,1.0;length=samples)
    values=ComplexF64[function_value(span*u+origin) for u in t]
    all(isfinite,values) || throw(DomainError(values,"nonfinite DCIM samples"))
    emptyfit=(residues=ComplexF64[],poles=ComplexF64[],origins=ComplexF64[])
    all(iszero,values) && return emptyfit
    pencil=samples÷2
    y1=[values[i+j-1] for i in 1:samples-pencil,j in 1:pencil]
    y2=[values[i+j] for i in 1:samples-pencil,j in 1:pencil]
    decomposition=svd(y1)
    cutoff=samples*eps(Float64)*first(decomposition.S)
    rank=min(terms,count(>(cutoff),decomposition.S))
    rank>0 || throw(ArgumentError("DCIM samples have no resolved pencil rank"))
    u=decomposition.U[:,1:rank]
    v=decomposition.V[:,1:rank]
    z=eigvals(Diagonal(inv.(decomposition.S[1:rank]))*(u'*y2*v))
    poles=log.(z)*(samples-1)/span
    # Growing exponentials cannot be integrated on the infinite physical path.
    filter!(p->isfinite(p)&&real(p)<0,poles)
    isempty(poles) && throw(ArgumentError("DCIM fit has no decaying image poles"))
    vandermonde=[exp(p*span*t_i) for t_i in t,p in poles]
    residues=vandermonde\values
    return (;residues,poles,origins=fill(ComplexF64(origin),length(poles)))
end

function _rallis_spectrum(fit,variable)
    return sum((c*exp(p*(variable-b)) for (c,p,b) in
        zip(fit.residues,fit.poles,fit.origins));init=0.0+0.0im)
end

function _rallis_fit(kind,kappa,depth,levels,samples,terms)
    if kind===:underground
        function_value=q->inv(sqrt(q^2-kappa^2)+q)
        origin=kappa
        endpoint=sqrt(100+kappa^2)
        splitpoint=sqrt(0.22^2+kappa^2)
    else
        function_value=λ->exp(-depth*sqrt(λ^2+kappa^2))/
            (λ+sqrt(λ^2+kappa^2))
        origin=0.0
        endpoint=levels==1 ? 28.0 : 100.0
        splitpoint=5.0
    end
    levels==1 && return _rallis_gpof(function_value,endpoint-origin,origin,samples,terms)
    far=_rallis_gpof(function_value,endpoint-splitpoint,splitpoint,samples,terms)
    remainder=q->function_value(q)-_rallis_spectrum(far,q)
    near=_rallis_gpof(remainder,splitpoint-origin,origin,samples,terms)
    return (residues=vcat(far.residues,near.residues),
        poles=vcat(far.poles,near.poles),origins=vcat(far.origins,near.origins))
end

function _rallis_functor(formula,base)
    state=base.state
    s=state.jω
    isfinite(s) && iszero(real(s)) && !iszero(s) ||
        throw(DomainError(s,"DCIM requires a finite nonzero real frequency"))
    scale=abs(state.gamma[2])
    isfinite(scale) && scale>0 ||
        throw(DomainError(scale,"DCIM requires a finite positive earth conductivity"))
    kappa=ComplexF64(state.gamma[2]/scale)
    a=formula.assumptions
    overhead=_rallis_fit(:overhead,kappa,0.0,a.dcim_levels,a.dcim_samples,a.dcim_terms)
    underground=_rallis_fit(:underground,kappa,0.0,a.dcim_levels,a.dcim_samples,a.dcim_terms)
    extended=merge(state,(dcim_scale=scale,dcim_kappa=kappa,
        dcim_overhead=overhead,dcim_underground=underground))
    return Functor{:Pollaczek1926,typeof(formula.routes),typeof(extended)}(
        formula.routes,extended)
end

function _rallis_correction(fit,kappa,height,lateral,::Val{:underground})
    T=typeof(height)
    total=complex(zero(T),zero(T))
    for (residue,pole,origin) in zip(fit.residues,fit.poles,fit.origins)
        c=Complex{T}(residue);p=Complex{T}(pole);b=Complex{T}(origin)
        H=height-p
        D=sqrt(H^2+lateral^2)
        argument=kappa*D
        # Keep the shifted residue and scaled K together; expanding the
        # source coefficient c_n first can overflow for deep complex images.
        total+=2c*kappa*H/D*special_besselkx(1,argument)*exp(-p*b-argument)
    end
    isfinite(total) || throw(DomainError(total,"nonfinite DCIM image sum"))
    return total
end

function _rallis_correction(fit,kappa,height,lateral,::Union{Val{:overhead},Val{:mixed}})
    T=typeof(height)
    total=complex(zero(T),zero(T))
    for (residue,pole,origin) in zip(fit.residues,fit.poles,fit.origins)
        c=Complex{T}(residue);p=Complex{T}(pole);b=Complex{T}(origin)
        H=height-p
        total+=2c*exp(-p*b)*H/(H^2+lateral^2)
    end
    isfinite(total) || throw(DomainError(total,"nonfinite DCIM image sum"))
    return total
end

function earth_impedance(::Val{:Pollaczek1926},::Val{:dcim},placement,functor,pair)
    _require(pair,placement)
    state=functor.state;geometry=_geometry(pair)
    scale=state.dcim_scale
    T=typeof(scale);kappa=state.gamma[2]/scale
    lateral=geometry.y_ij*scale
    if placement===Val(:overhead)
        self=pair.row==pair.column
        self && (lateral=zero(lateral))
        height=geometry.H*scale
        fit=state.dcim_overhead
        direct=self ? log(geometry.H/geometry.y_ij) :
            log(geometry.D_ij/geometry.d_ij)
    elseif placement===Val(:underground)
        height=geometry.H*scale
        fit=state.dcim_underground
        direct=bessel_difference(state.gamma[2],geometry.d_ij,geometry.D_ij)
    else
        air=pair.layers[1]==1 ? 1 : 2
        earth=3-air
        height=abs(pair.heights[air])*scale
        depth=Float64(abs(pair.heights[earth])*scale)
        a=state.formula.assumptions
        fit=_rallis_fit(:mixed,state.dcim_kappa,depth,a.dcim_levels,
            a.dcim_samples,a.dcim_terms)
        direct=zero(state.jω)
    end
    correction=_rallis_correction(fit,kappa,height,lateral,placement)
    return _complex_result(state.jω,
        state.jω*state.mu[1]/(2*(one(T)*π))*(direct+correction))
end
