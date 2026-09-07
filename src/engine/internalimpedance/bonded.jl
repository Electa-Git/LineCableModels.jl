# Algebraic elimination of a bonded cylindrical interface. Positive transfer
# impedances use the engine's opposing-surface current convention.
function _bonded_combine(first,second)
    denominator=first.outer+second.inner
    return (
        inner=first.inner-first.mutual^2/denominator,
        outer=second.outer-second.mutual^2/denominator,
        mutual=first.mutual*second.mutual/denominator)
end

function _bonded_work(layers,s)
    first_layer=first(layers)
    result=surface_impedances(Val(:Schelkunoff1934),first_layer.r_in,
        first_layer.r_ex,first_layer.rho,first_layer.mu_r,s)
    for layer in Iterators.drop(layers,1)
        next=surface_impedances(Val(:Schelkunoff1934),layer.r_in,
            layer.r_ex,layer.rho,layer.mu_r,s)
        result=_bonded_combine(result,next)
    end
    return result
end

function _bonded_surfaces(layers,s::Complex{T}) where {T <: AbstractFloat}
    isempty(layers) && throw(ArgumentError("bonded impedance needs physical radial layers"))
    all(isfinite,(real(s),imag(s))) && iszero(real(s)) ||
        throw(DomainError(s,"bonded impedance requires a finite real frequency"))
    for (n,layer) in pairs(layers)
        all(isfinite,(layer.r_in,layer.r_ex,layer.rho,layer.mu_r)) &&
            0<=layer.r_in<layer.r_ex && layer.rho>0 && layer.mu_r>0 ||
            throw(DomainError(layer,"invalid bonded-layer geometry or material"))
        n==1 || layers[n-1].r_ex==layer.r_in ||
            throw(ArgumentError("bonded layers must meet without a gap or overlap"))
    end
    resistance=[layer.rho/((one(T)*π)*(layer.r_ex-layer.r_in)*
        (layer.r_ex+layer.r_in)) for layer in layers]
    if iszero(s)
        outer=complex(inv(sum(inv,resistance)),zero(T))
        inside=iszero(first(layers).r_in) ? zero(outer) : outer
        return (inner=inside,outer=outer,mutual=inside)
    end
    length(layers)==1 && return _bonded_work(layers,s)
    # A poorly conducting screen may have a surface resistance many orders
    # above the final metal impedance. Guard the exact subtraction, not its
    # physical formula. Julia >=1.12 scopes this precision to the calling task.
    bits=precision(T)
    contrast=max(0,ceil(Int,log2(maximum(resistance)/minimum(resistance))))
    extra=maximum(layers) do layer
        skin=abs(s)*vacuum_permeability(layer.r_ex)*layer.mu_r*
            layer.r_ex^2/layer.rho
        max(0,ceil(Int,-log2(skin)))+
            max(0,ceil(Int,log2(layer.r_ex/(layer.r_ex-layer.r_in))))
    end
    wide=setprecision(BigFloat,bits+contrast+extra+48) do
        converted=[(r_in=BigFloat(layer.r_in),r_ex=BigFloat(layer.r_ex),
            rho=BigFloat(layer.rho),mu_r=BigFloat(layer.mu_r)) for layer in layers]
        _bonded_work(converted,Complex{BigFloat}(s))
    end
    narrow(value)=T===BigFloat ?
        complex(BigFloat(real(value);precision=bits),BigFloat(imag(value);precision=bits)) :
        Complex{T}(value)
    result=map(narrow,wide)
    all(isfinite,result) || throw(ErrorException("bonded surface evaluation was not finite"))
    return result
end
