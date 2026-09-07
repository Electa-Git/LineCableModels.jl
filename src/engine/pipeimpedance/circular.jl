function _validate_pipe(radius::T,outer_radius::T,rho::T,mur::T,s::Complex{T}) where {T}
    all(isfinite,(radius,outer_radius,rho,mur,s)) &&
        zero(T)<radius<outer_radius && rho>zero(T) && mur>zero(T) ||
        throw(DomainError((radius,outer_radius,rho,mur,s),
            "finite-pipe radii and material values must be positive and finite"))
    !iszero(s) && iszero(real(s)) ||
        throw(DomainError(s,"finite-pipe evaluation requires nonzero real frequency"))
    return nothing
end

function _geometry(pair::Pair{T}, radius::T) where {T}
    z1=complex(pair.positions[1]...)
    z2=complex(pair.positions[2]...)
    for (z,r) in zip((z1,z2),pair.radii)
        isfinite(z) && isfinite(r) && r>zero(T) && abs(z)+r<radius ||
            throw(DomainError((z,r), "each cable boundary must fit strictly inside the pipe"))
    end
    w=z1*conj(z2)/radius^2
    if pair.row==pair.column
        z1==z2 && pair.radii[1]==pair.radii[2] ||
            throw(ArgumentError("pipe self geometry must repeat the same cable"))
        geometric=log(radius/pair.radii[1])+log1p(-abs2(z1)/radius^2)
    else
        distance=abs(z1-z2)
        distance >= sum(pair.radii) ||
            throw(DomainError(distance, "enclosed cable boundaries must not overlap"))
        geometric=log(radius/distance)+log(abs(one(T)-w))
    end
    return (;w,geometric)
end

# At high order, direct Float64 K values can overflow even though their
# boundary ratios remain finite. Upward recurrence retains BigFloat range.
function _k_pair(n, x::Complex{BigFloat})
    previous=special_besselkx(0,x)
    current=special_besselkx(1,x)
    for k in 1:(n-1)
        previous,current=current,previous+2k/x*current
    end
    return current,previous
end

_k_pair(n,x) = (special_besselkx(n,x),special_besselkx(n-1,x))

function _finite_coefficient(n,x1,x2,mur)
    i1=special_besselix(n,x1); inext1=special_besselix(n+1,x1)
    i2=special_besselix(n,x2); inext2=special_besselix(n+1,x2)
    k1,kprev1=_k_pair(n,x1); k2,kprev2=_k_pair(n,x2)
    all(isfinite,(i1,inext1,i2,inext2,k1,kprev1,k2,kprev2)) || return nothing
    any(iszero,(i1,i2,k1,k2)) && return nothing
    # Recurrences eliminate derivatives without subtracting leading I terms.
    ai1=n*(mur-1)/x1*i1-inext1
    ai2=n*(mur+1)/x2*i2+inext2
    ak1=n*(mur+1)/x1*k1+kprev1
    ak2=n*(mur-1)/x2*k2-kprev2
    decay=exp(x1-x2+abs(real(x1))-abs(real(x2)))
    numerator=decay*ak2*i1-ai2*k1
    denominator=decay*ai1*ak2-ai2*ak1
    value=2mur/x1*numerator/denominator
    return isfinite(value) ? value : nothing
end

function _finite_harmonic(n,x1::Complex{T},x2::Complex{T},mur::T) where {T <: Real}
    value=try
        _finite_coefficient(n,x1,x2,mur)
    catch error
        error isa AmosException || rethrow()
        nothing
    end
    value === nothing || return Complex{T}(value)
    T <: Union{Float32,Float64} || throw(ErrorException(
        "finite-pipe Bessel ratio could not be evaluated at order $n"
    ))
    # This evaluates the same boundary ratio, not a large-order approximation.
    return setprecision(BigFloat,max(precision(BigFloat),256)) do
        wide=_finite_coefficient(n,Complex{BigFloat}(x1),Complex{BigFloat}(x2),BigFloat(mur))
        wide === nothing && throw(ErrorException(
            "finite-pipe Bessel ratio failed at order $n with wider precision"
        ))
        Complex{T}(wide)
    end
end

function _infinite_pipe_harmonic(n,x::Complex{T},mur::T) where {T}
    # Upward recurrence for K_(n-1)/K_n avoids large-order K overflow.
    ratio=special_besselkx(0,x)/special_besselkx(1,x)
    for k in 1:(n-1)
        ratio=inv(ratio+2k/x)
    end
    return Complex{T}(2mur/(n*(1+mur)+x*ratio))
end

function _static_pipe_harmonic(n,log_ratio,mur)
    contrast=(mur-1)/(mur+1)
    loss=-expm1(-2n*log_ratio)
    numerator=2/(mur+1)+contrast*loss
    denominator=4mur/(mur+1)^2+contrast^2*loss
    return 2mur/(n*(mur+1))*numerator/denominator
end

_harmonic_sum(state,w) = _harmonic_sum(
    n->_finite_harmonic(n,state.x1,state.x2,state.mur),state,w
)

function _harmonic_sum(coefficient_at,state,w)
    T=typeof(state.radius)
    iszero(w) && return zero(state.s)
    q=abs(w); power=one(w); result=zero(state.s); small=0
    for n in 1:state.max_terms
        power*=w
        coefficient=coefficient_at(n)
        result+=coefficient*real(power)
        # Test the amplitude before the angular cosine; a zero harmonic does
        # not establish convergence. The threshold is numerical, not a
        # source-supplied global error bound.
        envelope=abs(coefficient)*abs(power)*q/(1-q)
        if envelope <= state.rtol*max(one(T),abs(result))
            small+=1
            small>=3 && return result
        else
            small=0
        end
    end
    throw(ErrorException("finite-pipe harmonic sum reached max_terms=$(state.max_terms)"))
end
