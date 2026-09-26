"""
    PropagationParameters(modal; line_length=line_length(modal))

Bind per-unit-length modal parameters to one finite segment.
The source normalization length and active segment length are independent.
"""
struct PropagationParameters{P<:LineParameters,L<:Real,D<:ComputationDetails} <: AbstractCoreResult
    parameters::P
    line_length::L
    details::D

    function PropagationParameters(parameters::P,line_length::L,details::D) where
            {P<:LineParameters,L<:Real,D<:ComputationDetails}
        Engine.domain(parameters)===ModalDomain || throw(ArgumentError(
            "PropagationParameters require modal-domain coefficients"))
        basis(parameters)===:pul || throw(ArgumentError(
            "PropagationParameters require per-unit-length coefficients"))
        _segment_length(line_length)
        haskey(details.data,:gridpoint) && haskey(details.data,:source_gridpoint) ||
            throw(ArgumentError("segment details must retain gridpoint ancestry"))
        retained=get(details.data,:segment,nothing)
        retained isa NamedTuple && haskey(retained,:normalization_length) &&
            haskey(retained,:line_length) ||
            throw(ArgumentError("segment details must retain the completed length description"))
        isequal(retained.line_length,line_length) ||
            throw(ArgumentError("segment details must describe the stored line_length"))
        return new{P,L,D}(parameters,line_length,details)
    end
end

function _segment_length(value)
    value isa Real && !(value isa Bool) && isfinite(nominal(value)) &&
        nominal(value) >= 0 ||
        throw(ArgumentError("segment line_length must be finite and nonnegative"))
    return value
end

function per_unit_length(parameters::LineParameters{T,U,D,:pul}) where {T,U,D<:ModalDomain}
    return parameters
end

function per_unit_length(parameters::LineParameters{T,U,D,:total}) where {T,U,D<:ModalDomain}
    ell0 = line_length(parameters)
    ell0 isa Real && !(ell0 isa Bool) && isfinite(nominal(ell0)) &&
        nominal(ell0) > 0 || throw(ArgumentError(
            "total modal coefficients require a known positive source normalization length"))
    modal = ModalDomain(parameters.domain.operators, gamma(parameters) ./ ell0)
    return LineParameters(modal,
        SeriesImpedance(parameters.Z.values ./ ell0; basis=:pul),
        ShuntAdmittance(parameters.Y.values ./ ell0; basis=:pul),
        parameters.f,parameters.details)
end

function PropagationParameters(parameters::LineParameters{T,U,D};
        line_length=LineCableModels.line_length(parameters),combine::Symbol=:product) where {T,U,D<:ModalDomain}
    return parameterize(PropagationParameters,_bind_segment,(parameters,line_length);combine)
end

function _bind_segment(parameters::LineParameters{T,U,D},line_length) where {T,U,D<:ModalDomain}
    line_length === nothing && throw(ArgumentError(
        "segment line_length is required when the source has no declared length"))
    segment_length = _segment_length(line_length)
    per_unit_length_parameters = per_unit_length(parameters)
    source_length = LineCableModels.line_length(parameters)
    source_gridpoint = get(parameters.details.data,:gridpoint,nothing)
    segment_gridpoint = isequal(segment_length,source_length) ? source_gridpoint : Grammar.gridpoint_id()
    source_record=(; (key=>value for (key,value) in pairs(parameters.details.data)
        if key ∉ (:timing,:comparison_unsupported))...)
    segment_details = merge(source_record,
        (gridpoint=segment_gridpoint, source_gridpoint,
         segment=(line_length=segment_length,normalization_length=source_length,)))
    completed_details = ComputationDetails(segment_details)
    return PropagationParameters(per_unit_length_parameters,segment_length,completed_details)
end

function PropagationParameters(source::PropagationParameters;
        line_length=source.line_length,combine::Symbol=:product)
    return parameterize(PropagationParameters,_bind_segment,(source,line_length);combine)
end

function _bind_segment(source::PropagationParameters,line_length)
    segment_length = _segment_length(line_length)
    if isequal(segment_length,source.line_length)
        return source
    end
    source_record=(; (key=>value for (key,value) in pairs(source.details.data)
        if key ∉ (:timing,:comparison_unsupported))...)
    segment_details = merge(source_record,
        (gridpoint=Grammar.gridpoint_id(),source_gridpoint=get(source.details.data,:gridpoint,nothing),
         segment=merge(source.details.data.segment,(line_length=segment_length,))))
    return PropagationParameters(source.parameters,segment_length,ComputationDetails(segment_details))
end

line_length(source::PropagationParameters) = source.line_length
basis(::PropagationParameters) = :pul
frequencies(source::PropagationParameters) = frequencies(source.parameters)
details(source::PropagationParameters) = source.details
Tv(source::PropagationParameters) = Tv(source.parameters)
Ti(source::PropagationParameters) = Ti(source.parameters)
gamma(source::PropagationParameters) = gamma(source.parameters)

"""Return modal attenuation constants in inverse metres, ordered mode × frequency."""
alpha(source::PropagationParameters) = real.(gamma(source))

"""Return modal phase constants in radians per metre, ordered mode × frequency."""
beta(source::PropagationParameters) = imag.(gamma(source))

"""Return phase velocities in metres per second from `2πf / beta`, ordered mode × frequency."""
function velocity(source::PropagationParameters)
    constants=beta(source)
    f=frequencies(source)
    return map(CartesianIndices(constants)) do index
        value=constants[index]
        T=promote_type(typeof(float(real(nominal(value)))),typeof(float(nominal(f[index[2]]))))
        (T(2)*T(π)*f[index[2]])/value
    end
end
Zc(source::PropagationParameters, domain::Type{<:Engine.LineParamsDomain}=ModalDomain) =
    Zc(source.parameters,domain)
Yc(source::PropagationParameters, domain::Type{<:Engine.LineParamsDomain}=ModalDomain) =
    Yc(source.parameters,domain)

function Base.getindex(source::PropagationParameters,
        selected::Union{Integer,AbstractRange{<:Integer},AbstractVector{<:Integer},Colon})
    parameters=source.parameters[selected]
    retained=ComputationDetails(merge(source.details.data,
        (modal=get(parameters.details.data,:modal,nothing),)))
    return PropagationParameters(parameters,source.line_length,retained)
end

"""Forward-wave factors over the segment's bound physical length."""
H(source::PropagationParameters) = H(source,ModalDomain)

function H(source::PropagationParameters, ::Type{ModalDomain};field=nothing)
    field === nothing || throw(ArgumentError("modal H does not accept a field"))
    return exp.(-gamma(source).*source.line_length)
end

function H(source::PropagationParameters, ::Type{PhaseDomain};field)
    field in (:voltage,:current) || throw(ArgumentError("field must be :voltage or :current"))
    return _phase_H(source,Val(field))
end

function _phase_H(source::PropagationParameters, ::Val{:voltage})
    return _similarity_H(H(source),Tv(source))
end
function _phase_H(source::PropagationParameters, ::Val{:current})
    return _similarity_H(H(source),Ti(source))
end

function _similarity_H(factors, maps)
    n,nf = size(factors)
    S=promote_type(eltype(factors),eltype(maps))
    result=Array{S,3}(undef,n,n,nf)
    for frequency in 1:nf
        basis_matrix=@view maps[:,:,frequency]
        @views result[:,:,frequency] .=
            (basis_matrix * Diagonal(factors[:,frequency])) / basis_matrix
    end
    return result
end
