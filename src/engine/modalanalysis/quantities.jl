"""Return the retained modal-to-phase voltage basis (phase × mode × frequency)."""
Tv(parameters::LineParameters{T,U,D}) where {T,U,D<:ModalDomain} =
    parameters.domain.operators.Tv

"""Return the retained modal-to-phase current basis (phase × mode × frequency)."""
Ti(parameters::LineParameters{T,U,D}) where {T,U,D<:ModalDomain} =
    parameters.domain.operators.Ti

"""Return the retained propagation roots (mode × frequency)."""
gamma(parameters::LineParameters{T,U,D}) where {T,U,D<:ModalDomain} =
    parameters.domain.gamma

function _characteristic(parameters::LineParameters{T,U,D}, quantity::Val,
        ::Type{ModalDomain}) where {T,U,D<:ModalDomain}
    coefficients = quantity isa Val{:Zc} ? parameters.Z.values : parameters.Y.values
    roots = gamma(parameters)
    n, _, nf = size(coefficients)
    S = promote_type(T,eltype(roots))
    result = Array{S,2}(undef,n,nf)
    for frequency in 1:nf, mode in 1:n
        result[mode,frequency] = coefficients[mode,mode,frequency] / roots[mode,frequency]
    end
    return result
end

function _characteristic(parameters::LineParameters{T,U,D}, quantity::Val,
        ::Type{PhaseDomain}) where {T,U,D<:ModalDomain}
    diagonal = _characteristic(parameters,quantity,ModalDomain)
    maps = operators(parameters)
    n,nf = size(diagonal)
    S = promote_type(eltype(diagonal),eltype(maps.Tv),eltype(maps.Ti))
    result = Array{S,3}(undef,n,n,nf)
    for frequency in 1:nf
        if quantity isa Val{:Zc}
            left = @view maps.Tv[:,:,frequency]
            right = @view maps.Ti[:,:,frequency]
        else
            left = @view maps.Ti[:,:,frequency]
            right = @view maps.Tv[:,:,frequency]
        end
        @views result[:,:,frequency] .= (left * Diagonal(diagonal[:,frequency])) / right
    end
    return result
end

"""Characteristic impedance of the selected diagonal modal approximation."""
Zc(parameters::LineParameters{T,U,D}, domain::Type{<:Engine.LineParamsDomain}=ModalDomain) where {T,U,D<:ModalDomain} =
    _characteristic(parameters,Val(:Zc),domain)

"""Characteristic admittance of the selected diagonal modal approximation."""
Yc(parameters::LineParameters{T,U,D}, domain::Type{<:Engine.LineParamsDomain}=ModalDomain) where {T,U,D<:ModalDomain} =
    _characteristic(parameters,Val(:Yc),domain)
