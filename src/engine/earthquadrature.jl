"""
Evaluate a semi-infinite integral with the double-exponential node sum
used by [Pires2026](@cite), equations (12)–(14).

The positive interval map is λ=(1+x)/(1-x), dλ/dx=2/(1-x)².
Composing it analytically with x=tanh(π*sinh(t)/2) gives
λ=exp(π*sinh(t)); no rounded endpoint subtraction is performed.

The implementation halves the mesh, sums both tails, and checks two
successive refinements. These stopping controls are numerical choices,
not a source-wide error bound or the paper's particular flowchart.
The returned error estimate is the change between successive sums.
Failure to establish convergence is reported rather than replaced
silently by a different quadrature.
"""
function double_exponential(integrand,::Type{T};rtol=T(1e-8),atol=zero(T),
        maxlevel::Int=12,breakpoints=()) where {T <: AbstractFloat}
    rtol>0 && isfinite(rtol) && atol>=0 && isfinite(atol) && maxlevel>=3 ||
        throw(ArgumentError("double-exponential quadrature requires positive tolerance and at least three levels"))
    if !isempty(breakpoints)
        points=sort!(unique!(T[breakpoints...]))
        all(x->isfinite(x) && x>0,points) ||
            throw(ArgumentError("quadrature breakpoints must be finite and positive"))
        lower=zero(T); values=Complex{T}[]; errors=T[]; level=0
        for upper in points
            width=upper-lower
            mapped=t->begin
                v=inv(one(T)+t)
                x=t<=1 ? lower+width*t*v : upper-width*v
                (x==lower || x==upper) && return zero(Complex{T})
                integrand(x)*width*v*v
            end
            part=double_exponential(mapped,T;rtol,atol,maxlevel)
            push!(values,part.value); push!(errors,part.error)
            level=max(level,part.level); lower=upper
        end
        tail=t->begin
            x=lower+t
            x==lower && return zero(Complex{T})
            integrand(x)
        end
        part=double_exponential(tail,T;rtol,atol,maxlevel)
        push!(values,part.value); push!(errors,part.error)
        return (value=sum(values),error=sum(errors),level=max(level,part.level))
    end
    # BigFloat's exponent range is vastly larger than its significand precision.
    # Visiting its entire range can make oscillatory argument reduction unbounded.
    # Retain precision-scaled endpoints and require the omitted tails to be small.
    limit=min(log(floatmax(T))-T(4),max(T(64),-T(4)*log(eps(T))))
    extent=asinh(limit/T(π))
    previous=nothing; previous_error=T(Inf); successes=0
    for level in 1:maxlevel
        h=ldexp(one(T),-level); count=ceil(Int,extent/h)
        central=integrand(one(T))*T(π)*h
        total=central; compensation=zero(central)
        magnitude=abs(central); left_tail=zero(T); right_tail=zero(T)
        for k in 1:count
            t=T(k)*h; u=T(π)*sinh(t)
            u>=limit && break
            jacobian=T(π)*cosh(t)*h
            left=exp(-u); right=exp(u)
            additions=(integrand(left)*left*jacobian,integrand(right)*right*jacobian)
            left_tail=abs(additions[1]); right_tail=abs(additions[2])
            for term in additions
                isfinite(term) || throw(DomainError(term,"nonfinite double-exponential quadrature term"))
                adjusted=term-compensation
                next=total+adjusted
                compensation=(next-total)-adjusted
                total=next; magnitude+=abs(term)
            end
        end
        isfinite(total) || throw(DomainError(total,"nonfinite double-exponential integral"))
        target=max(atol,rtol*abs(total),8eps(T)*magnitude)
        tails=left_tail+right_tail
        if previous!==nothing
            previous_error=abs(total-previous)
            converged=previous_error<=target && tails<=target
            successes=converged ? successes+1 : 0
            successes>=2 && return (value=total,error=max(previous_error,tails),level)
        end
        previous=total
    end
    throw(ErrorException("double-exponential quadrature did not converge within $(maxlevel) levels (last change $(previous_error))"))
end
