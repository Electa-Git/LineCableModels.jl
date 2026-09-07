# Pettersson (1994), equations (10), (11), (14), and (15).
# The same source-medium selection and branch convention serve Z and P.
function pettersson_images(state,pair)
    s=state.jω
    iszero(s) && throw(DomainError(s,"Pettersson images require nonzero frequency"))
    h1,h2=abs.(pair.heights)
    on_interface=iszero(h1)&&iszero(h2)
    if !on_interface
        iszero(h1)||iszero(h2) ?
            throw(ArgumentError("interface/off-interface pairs require a common modal prescription")) : nothing
        pair.layers[1]==pair.layers[2] || throw(ArgumentError(
            "Pettersson's closed images require both wires in the same half-space"))
    end
    medium=on_interface ? 1 : pair.layers[1]
    medium in (1,2) || throw(ArgumentError("Pettersson images require two half-spaces"))
    other=3-medium
    kappa=state.sigma[medium]+s*state.epsilon[medium]
    ratio=(state.sigma[other]+s*state.epsilon[other])/kappa
    beta=spectral_root(state.gamma_medium_squared[other]-
        state.gamma_medium_squared[medium],s)
    iszero(beta) && throw(DomainError(beta,
        "identical half-spaces do not define a finite image-reference distance"))
    self=pair.row==pair.column
    x=self ? zero(h1) : pair.separation
    direct=self ? pair.separation : hypot(x,h1-h2)
    direct>zero(direct) || throw(DomainError(direct,"wire distance must be positive"))
    if on_interface
        unit=complex(zero(h1),copysign(one(h1),imag(s)))
        shiftP=sqrt(2one(h1))*(1+unit)/beta
        shiftQ=unit*sqrt(2one(h1))*(ratio+1)/(ratio+unit)/beta
        # Equation (15) neglects the unit term for |beta*a| << 1.
        self && abs(beta*direct)>=1 && throw(DomainError(beta*direct,
            "the on-interface self approximation requires |beta*a| much smaller than one"))
        imageP=self ? shiftP : sqrt(shiftP^2+x^2)
        imageQ=self ? shiftQ : sqrt(shiftQ^2+x^2)
        mirror=direct
        ideal=zero(s)
    else
        H=h1+h2
        mirror=self ? H : hypot(x,H)
        imageP=self ? H+2/beta : sqrt((H+2/beta)^2+x^2)
        imageQ=self ? H+(ratio+1)/beta : sqrt((H+(ratio+1)/beta)^2+x^2)
        ideal=log(mirror/direct)
    end
    # The printed Q image has Im(dQ)<0 in air, >0 in ground (positive f).
    sense=(medium==1 ? -one(h1) : one(h1))*sign(imag(s))
    imag(imageQ)*sense<0 && (imageQ=-imageQ)
    magnetic=ideal+log(imageP/mirror)
    electric=ideal+2/(ratio+1)*log(imageQ/mirror)
    return (;magnetic,electric,kappa)
end
