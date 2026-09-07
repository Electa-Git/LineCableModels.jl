"""
Physical round-core data for a pipe proximity calculation.
Indices are coaxial assembly labels; radii are metal, not coating, radii.
"""
struct Cores{T <: Real}
    indices::Vector{Int}
    positions::Vector{Tuple{T,T}}
    radii::Vector{T}
    rho::Vector{T}
    mur::Vector{T}
    function Cores(indices,positions::Vector{Tuple{T,T}},radii::Vector{T},
            rho::Vector{T},mur::Vector{T}) where {T <: Real}
        n=length(indices)
        n>0 && all(x->x isa Integer && x>0,indices) && length(unique(indices))==n &&
            all(length(x)==n for x in (positions,radii,rho,mur)) ||
            throw(ArgumentError("core arrays must have equal nonzero lengths and unique indices"))
        all(x->isfinite(x) && x>0,Iterators.flatten((radii,rho,mur))) &&
            all(p->all(isfinite,p),positions) ||
            throw(DomainError((radii,rho,mur),"core geometry and material values must be finite and positive"))
        for i in 1:n,j in i+1:n
            hypot((positions[i].-positions[j])...)>=radii[i]+radii[j] ||
                throw(DomainError((i,j),"metal core sections must not overlap"))
        end
        new{T}(collect(Int,indices),copy(positions),copy(radii),copy(rho),copy(mur))
    end
end

function with_cores(functor::Functor{ID},cores::Cores) where {ID}
    selection=functor.state.proximity
    selection===:none && throw(ArgumentError("select a proximity formula before supplying core data"))
    all(==(first(cores.radii)),cores.radii) &&
        all(==(first(cores.rho)),cores.rho) &&
        all(==(first(cores.mur)),cores.mur) ||
        throw(ArgumentError("reciprocal pipe assembly requires equal core radii and materials; the source's unequal-core terms are directional"))
    for (position,radius) in zip(cores.positions,cores.radii)
        hypot(position...)+radius<functor.state.radius ||
            throw(DomainError(position,"each metal core must fit inside the pipe"))
    end
    if selection===:Hoidalen2013
        length(cores.indices)==3 && all(isone,cores.mur) ||
            throw(ArgumentError("Hoidalen proximity matrix requires three equal nonmagnetic cores"))
        d=[hypot((cores.positions[i].-cores.positions[j])...) for (i,j) in ((1,2),(1,3),(2,3))]
        all(x->isapprox(x,d[1];rtol=64eps(typeof(x))),d) ||
            throw(ArgumentError("Hoidalen's published mode factors require symmetrical three-core spacing"))
    end
    state=merge(functor.state,(;cores))
    return Functor{ID,typeof(functor.routes),typeof(state)}(functor.routes,state)
end

# I_n(z)/(z I_(n-1)(z)), with its difference from 1/(2n).
# The continued fraction retains the low-frequency difference directly.
function _core_i_ratio(n,z::Complex{T}) where {T}
    if abs(z)<one(T)/2
        q=zero(z)
        for order in n+max(24,precision(T)÷4+8):-1:n+1
            q=inv(2order+z*z*q)
        end
        correction=z*z*q
        denominator=2n+correction
        return (ratio=inv(denominator),difference=-correction/(2n*denominator))
    end
    numerator=special_besselix(n,z); denominator=z*special_besselix(n-1,z)
    if iszero(denominator) || !isfinite(numerator/denominator)
        q=zero(z)
        for order in n+ceil(Int,abs(z))+max(32,precision(T)÷4):-1:n
            q=inv(2order+z*z*q)
        end
    else
        q=Complex{T}(numerator/denominator)
    end
    return (ratio=q,difference=q-inv(T(2n)))
end

"""
Evaluate the directional solid-core proximity increment of
[Kane1995](@cite), (10),(17), or [Hoidalen2013](@cite), (36).
This is not the isolated skin term or a general unequal-core matrix.
"""
function core_proximity(::Val{ID},radius::T,distance::T,rho::T,mur::T,
        s::Complex{T};rtol::Real=1e-10,max_terms::Int=1024) where {ID,T <: Real}
    ID in (:Kane1995,:Hoidalen2013) || throw(ArgumentError("unknown core proximity formula"))
    all(isfinite,(radius,distance,rho,mur,s)) &&
        0<radius<distance && rho>0 && mur>0 && iszero(real(s)) ||
        throw(DomainError((radius,distance,rho,mur,s),"invalid solid-core proximity parameters"))
    0<rtol<1 && max_terms>0 || throw(ArgumentError("invalid proximity sum controls"))
    ID===:Hoidalen2013 && !isone(mur) &&
        throw(DomainError(mur,"Hoidalen (36) requires nonmagnetic cores"))
    iszero(s) && return zero(s)
    piT=one(T)*π; mu0=T(4)*piT/T(10)^7
    z=sqrt(s*mu0*mur/rho)*radius
    g=(radius/distance)^2; power=one(T); result=zero(s); small=0
    q1=_core_i_ratio(1,z).ratio; q2=_core_i_ratio(2,z).ratio
    k1=g*z*z*q1*q2; A=k1/(1-k1)
    for n in 1:max_terms
        power*=g
        ratio=_core_i_ratio(n,z)
        coefficient=if ID===:Kane1995
            2mur*ratio.ratio/(1+n*(mur-1)*ratio.ratio)
        else
            # Subtract the entire static logarithm analytically, term by term.
            excess=2n*real(A)+n*n*abs2(A)
            ratio.difference*(1+excess)+excess/(2n)
        end
        term=power*coefficient
        result+=term
        if abs(term)*g/(1-g)<=max(T(rtol),eps(T))*abs(result)
            small+=1
            small>=3 && return Complex{T}(s*mu0/(2piT)*result)
        else
            small=0
        end
    end
    throw(ErrorException("core proximity sum reached max_terms=$max_terms"))
end

function _core_increment(functor,pair)
    state=functor.state
    state.proximity===:none && return zero(state.s)
    hasproperty(state,:cores) || throw(ArgumentError("this pipe formula requires physical core data through with_cores"))
    cores=state.cores
    i=findfirst(==(pair.row),cores.indices); j=findfirst(==(pair.column),cores.indices)
    (i===nothing || j===nothing) && throw(ArgumentError("pipe pair indices are absent from core data"))
    cores.positions[i]==pair.positions[1] && cores.positions[j]==pair.positions[2] &&
        cores.radii[i]<=pair.radii[1] && cores.radii[j]<=pair.radii[2] ||
        throw(ArgumentError("pipe pair and physical core geometry disagree"))
    function increment(receiver,scatterer)
        d=hypot((cores.positions[receiver].-cores.positions[scatterer])...)
        return core_proximity(Val(state.proximity),cores.radii[scatterer],d,
            cores.rho[scatterer],cores.mur[scatterer],state.s;
            rtol=state.rtol,max_terms=state.max_terms)
    end
    if i==j
        return sum((increment(i,k) for k in eachindex(cores.indices) if k!=i);init=zero(state.s))
    end
    # Source Z_(j,i) uses the receiver core's material; reciprocal assembly
    # has already rejected unequal cores instead of averaging their values.
    return increment(j,i)
end
