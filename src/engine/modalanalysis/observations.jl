const _phase_Zc = Base.Fix2(Zc,(domain=PhaseDomain,))
const _phase_Yc = Base.Fix2(Yc,(domain=PhaseDomain,))
const _voltage_H = Base.Fix2(H,(domain=PhaseDomain,field=:voltage))
const _current_H = Base.Fix2(H,(domain=PhaseDomain,field=:current))

"""
$(TYPEDEF)

The name of the voltage propagation function, `H` bound to the phase domain and the
voltage field. The binding stores the field as a value, so the voltage and the current
functions share one selector type.
"""
abstract type VoltagePropagationFunction end

"""
$(TYPEDEF)

The name of the current propagation function, `H` bound to the phase domain and the
current field. See [`VoltagePropagationFunction`](@ref).
"""
abstract type CurrentPropagationFunction end

# Scientific units belong to the scalar quantity, independently of presentation. A
# selector's type names its quantity.
Units.quantity(selector::Union{typeof(gamma),typeof(alpha),typeof(beta),typeof(velocity),
    typeof(Zc),typeof(Yc),typeof(H),typeof(Tv),typeof(Ti),typeof(_phase_Zc),typeof(_phase_Yc)}) =
    Units.Quantity{typeof(selector), Nothing}()
Units.quantity(value::Base.Fix2{typeof(H)}) = value.x.field===:voltage ?
    Units.Quantity{VoltagePropagationFunction, Nothing}() :
    Units.Quantity{CurrentPropagationFunction, Nothing}()

Units.native_unit(::Units.Quantity{typeof(gamma), Nothing}) =
    Units.units(:base,:dimensionless;per=(:base,:meter))
Units.display_unit(::Units.Quantity{typeof(gamma), Nothing}) =
    Units.units(:base,:dimensionless;per=(:kilo,:meter))
Units.label(::Units.Quantity{typeof(gamma), Nothing}) = "Propagation constant"
Units.symbol(::Units.Quantity{typeof(gamma), Nothing}) = "γ"

Units.native_unit(::Units.Quantity{typeof(alpha), Nothing}) = Units.units(:base,:neper;per=(:base,:meter))
Units.display_unit(::Units.Quantity{typeof(alpha), Nothing}) = Units.units(:base,:neper;per=(:kilo,:meter))
Units.label(::Units.Quantity{typeof(alpha), Nothing}) = "Attenuation constant"
Units.symbol(::Units.Quantity{typeof(alpha), Nothing}) = "α"

Units.native_unit(::Units.Quantity{typeof(beta), Nothing}) = Units.units(:base,:radian;per=(:base,:meter))
Units.display_unit(::Units.Quantity{typeof(beta), Nothing}) = Units.units(:base,:radian;per=(:kilo,:meter))
Units.label(::Units.Quantity{typeof(beta), Nothing}) = "Phase constant"
Units.symbol(::Units.Quantity{typeof(beta), Nothing}) = "β"

Units.native_unit(::Units.Quantity{typeof(velocity), Nothing}) = Units.units(:base,:meter;per=(:base,:second))
Units.display_unit(::Units.Quantity{typeof(velocity), Nothing}) = Units.native_unit(Units.Quantity{typeof(velocity), Nothing}())
Units.label(::Units.Quantity{typeof(velocity), Nothing}) = "Phase velocity"
Units.symbol(::Units.Quantity{typeof(velocity), Nothing}) = "vₚ"

Units.native_unit(::Units.Quantity{typeof(Zc), Nothing}) = Units.units(:base,:ohm)
Units.display_unit(::Units.Quantity{typeof(Zc), Nothing}) = Units.units(:base,:ohm)
Units.label(::Units.Quantity{typeof(Zc), Nothing}) = "Characteristic impedance"
Units.symbol(::Units.Quantity{typeof(Zc), Nothing}) = "Zc"

Units.native_unit(::Units.Quantity{typeof(Yc), Nothing}) = Units.units(:base,:siemens)
Units.display_unit(::Units.Quantity{typeof(Yc), Nothing}) = Units.units(:base,:siemens)
Units.label(::Units.Quantity{typeof(Yc), Nothing}) = "Characteristic admittance"
Units.symbol(::Units.Quantity{typeof(Yc), Nothing}) = "Yc"

Units.native_unit(::Units.Quantity{typeof(H), Nothing}) = Units.units(:base,:dimensionless)
Units.display_unit(::Units.Quantity{typeof(H), Nothing}) = Units.units(:base,:dimensionless)
Units.label(::Units.Quantity{typeof(H), Nothing}) = "Propagation function"
Units.symbol(::Units.Quantity{typeof(H), Nothing}) = "H"

Units.native_unit(::Units.Quantity{typeof(Tv), Nothing}) = Units.units(:base,:dimensionless)
Units.display_unit(::Units.Quantity{typeof(Tv), Nothing}) = Units.units(:base,:dimensionless)
Units.label(::Units.Quantity{typeof(Tv), Nothing}) = "Modal-to-phase voltage transformation"
Units.symbol(::Units.Quantity{typeof(Tv), Nothing}) = "Tv"

Units.native_unit(::Units.Quantity{typeof(Ti), Nothing}) = Units.units(:base,:dimensionless)
Units.display_unit(::Units.Quantity{typeof(Ti), Nothing}) = Units.units(:base,:dimensionless)
Units.label(::Units.Quantity{typeof(Ti), Nothing}) = "Modal-to-phase current transformation"
Units.symbol(::Units.Quantity{typeof(Ti), Nothing}) = "Ti"

# The phase-domain variants keep the unit of the modal quantity they represent.
Units.native_unit(::Units.Quantity{typeof(_phase_Zc), Nothing}) =
    Units.native_unit(Units.Quantity{typeof(Zc), Nothing}())
Units.display_unit(::Units.Quantity{typeof(_phase_Zc), Nothing}) =
    Units.display_unit(Units.Quantity{typeof(Zc), Nothing}())
Units.label(::Units.Quantity{typeof(_phase_Zc), Nothing}) = "Phase-domain characteristic impedance"
Units.symbol(::Units.Quantity{typeof(_phase_Zc), Nothing}) = "Zc,phase"

Units.native_unit(::Units.Quantity{typeof(_phase_Yc), Nothing}) =
    Units.native_unit(Units.Quantity{typeof(Yc), Nothing}())
Units.display_unit(::Units.Quantity{typeof(_phase_Yc), Nothing}) =
    Units.display_unit(Units.Quantity{typeof(Yc), Nothing}())
Units.label(::Units.Quantity{typeof(_phase_Yc), Nothing}) = "Phase-domain characteristic admittance"
Units.symbol(::Units.Quantity{typeof(_phase_Yc), Nothing}) = "Yc,phase"

Units.native_unit(::Units.Quantity{VoltagePropagationFunction, Nothing}) =
    Units.native_unit(Units.Quantity{typeof(H), Nothing}())
Units.display_unit(::Units.Quantity{VoltagePropagationFunction, Nothing}) =
    Units.display_unit(Units.Quantity{typeof(H), Nothing}())
Units.label(::Units.Quantity{VoltagePropagationFunction, Nothing}) = "Voltage propagation function"
Units.symbol(::Units.Quantity{VoltagePropagationFunction, Nothing}) = "Hᵥ"

Units.native_unit(::Units.Quantity{CurrentPropagationFunction, Nothing}) =
    Units.native_unit(Units.Quantity{typeof(H), Nothing}())
Units.display_unit(::Units.Quantity{CurrentPropagationFunction, Nothing}) =
    Units.display_unit(Units.Quantity{typeof(H), Nothing}())
Units.label(::Units.Quantity{CurrentPropagationFunction, Nothing}) = "Current propagation function"
Units.symbol(::Units.Quantity{CurrentPropagationFunction, Nothing}) = "Hᵢ"

# Named parts. The real and imaginary parts of γ are the attenuation and phase constants,
# with units of their own. Those of Zc and Yc take the unit of the Units rule and carry
# their own names. Every other component follows the Units rules.
Units.quantity(::typeof(gamma),::typeof(real)) = Units.quantity(alpha)
Units.quantity(::typeof(gamma),::typeof(imag)) = Units.quantity(beta)
Units.label(::Units.Quantity{typeof(Zc), typeof(real)}) = "Characteristic resistance"
Units.label(::Units.Quantity{typeof(Zc), typeof(imag)}) = "Characteristic reactance"
Units.label(::Units.Quantity{typeof(Yc), typeof(real)}) = "Characteristic conductance"
Units.label(::Units.Quantity{typeof(Yc), typeof(imag)}) = "Characteristic susceptance"
Units.symbol(::Units.Quantity{typeof(Zc), typeof(real)}) = "R꜀"
Units.symbol(::Units.Quantity{typeof(Zc), typeof(imag)}) = "X꜀"
Units.symbol(::Units.Quantity{typeof(Yc), typeof(real)}) = "G꜀"
Units.symbol(::Units.Quantity{typeof(Yc), typeof(imag)}) = "B꜀"


# Finite representation bindings include only domain and, for H, field.
function _validate_phase_selector(selector::Base.Fix2)
    selector.f in (H,Zc,Yc) || throw(ArgumentError("unsupported finite representation selector"))
    value=selector.x
    value isa NamedTuple || throw(ArgumentError("representation binding must be a NamedTuple"))
    if selector.f===H
        keys(value)==(:domain,:field) && value.domain===PhaseDomain &&
            value.field in (:voltage,:current) || throw(ArgumentError("invalid phase H binding"))
    else
        keys(value)==(:domain,) && value.domain===PhaseDomain ||
            throw(ArgumentError("invalid phase characteristic binding"))
    end
    return selector
end

const _ModalLineParameters = LineParameters{T,U,D} where {T<:Complex,U<:Real,D<:ModalDomain}
for selector in (gamma,Zc,Yc,Tv,Ti)
    @eval function observe(source::_ModalLineParameters,::typeof($selector),indices...)
        values=$selector(source)
        return isempty(indices) ? values : getindex(values,indices...)
    end
end
for selector in (gamma,alpha,beta,velocity,Zc,Yc,H,Tv,Ti)
    @eval function observe(source::PropagationParameters,::typeof($selector),indices...)
        values=$selector(source)
        return isempty(indices) ? values : getindex(values,indices...)
    end
end
function observe(source::Union{_ModalLineParameters,PropagationParameters},
        selector::Base.Fix2{typeof(Zc)},indices...)
    _validate_phase_selector(selector)
    values=Zc(source,PhaseDomain)
    return isempty(indices) ? values : getindex(values,indices...)
end
function observe(source::Union{_ModalLineParameters,PropagationParameters},
        selector::Base.Fix2{typeof(Yc)},indices...)
    _validate_phase_selector(selector)
    values=Yc(source,PhaseDomain)
    return isempty(indices) ? values : getindex(values,indices...)
end
function observe(source::PropagationParameters,
        selector::Base.Fix2{typeof(H)},indices...)
    _validate_phase_selector(selector)
    values=H(source,PhaseDomain;field=selector.x.field)
    return isempty(indices) ? values : getindex(values,indices...)
end
const _ModalTransform=Union{typeof(abs),typeof(angle),typeof(real),typeof(imag)}
for selector in (gamma,Zc,Yc,Tv,Ti), source_type in (:_ModalLineParameters,:PropagationParameters)
    @eval observe(source::$source_type,::typeof($selector),
        transform::_ModalTransform,indices...) = transform.(observe(source,$selector,indices...))
end
@eval observe(source::PropagationParameters,::typeof(H),
    transform::_ModalTransform,indices...) = transform.(observe(source,H,indices...))
for selector in (Zc,Yc), source_type in (:_ModalLineParameters,:PropagationParameters)
    @eval observe(source::$source_type,bound::Base.Fix2{typeof($selector)},
        transform::_ModalTransform,indices...) = transform.(observe(source,bound,indices...))
end
observe(source::PropagationParameters,bound::Base.Fix2{typeof(H)},
    transform::_ModalTransform,indices...) = transform.(observe(source,bound,indices...))

_modal_transform_requests(selectors) = Tuple((selector,transform) for selector in selectors
    for transform in (abs,angle,real,imag))

function observables(::Type{<:LineParameters{T,U,D}}) where {T<:Complex,U<:Real,D<:ModalDomain}
    selectors=(gamma,Zc,Yc,Tv,Ti,_phase_Zc,_phase_Yc)
    return (observables(LineParameters)...,selectors...,
        _modal_transform_requests(selectors)...)
end
function observables(::Type{<:PropagationParameters})
    selectors=(gamma,Zc,Yc,H,Tv,Ti,_phase_Zc,_phase_Yc,_voltage_H,_current_H)
    return (selectors...,alpha,beta,velocity,_modal_transform_requests(selectors)...)
end

Commons.normalize_observation_selector(::typeof(alpha)) = (gamma,real)
Commons.normalize_observation_selector(::typeof(beta)) = (gamma,imag)
Commons.normalize_observation_selector(::Val{:alpha}) = Commons.normalize_observation_selector(alpha)
Commons.normalize_observation_selector(::Val{:beta}) = Commons.normalize_observation_selector(beta)

function _normalize_modal_request(source,request)
    identity=request_identity(request)
    prefix=identity isa Tuple ? identity : (identity,)
    selector=first(prefix)
    selector in observables(typeof(source)) ||
        selector in (gamma,alpha,beta,velocity,Zc,Yc,H,Tv,Ti) ||
        throw(ArgumentError("unsupported modal quantity"))
    selector isa Base.Fix2 && _validate_phase_selector(selector)
    length(prefix)<=2 && (length(prefix)==1 || prefix[2] in (abs,angle,real,imag)) ||
        throw(ArgumentError("unsupported modal transform"))
    rank=selector in (Tv,Ti,_phase_Zc,_phase_Yc,_voltage_H,_current_H) ? 3 : 2
    indices=request_indices(request)
    isempty(indices) && (indices=ntuple(_->Colon(),rank))
    length(indices)==rank || throw(DimensionMismatch("modal request requires $rank selectors"))
    dimensions=rank==2 ? size(gamma(source)) : size(Tv(source))
    map(observation_indices,indices,dimensions)
    return (prefix...,indices...)
end

function Commons.observation_requests(source::_ModalLineParameters,requests::Tuple;complete_pairs::Bool=false)
    selected=isempty(requests) ? (Engine.Z,Engine.Y,gamma,Zc,Yc,Tv,Ti) : requests
    primary=Tuple(filter(request -> begin
        identity=request_identity(request)
        first(identity isa Tuple ? identity : (identity,)) in
            (Engine.Z,Engine.Y,Engine.R,Engine.X,Engine.L,Engine.G,Engine.B,Engine.C,real,imag,abs,angle)
    end,selected))
    derived=Tuple(filter(request -> request ∉ primary,selected))
    normalized_primary=isempty(primary) ? () :
        invoke(Commons.observation_requests,
            Tuple{LineParameters,Tuple},
            source,primary;complete_pairs).retained
    normalized_derived=_modal_component_requests(source,derived;complete_pairs)
    retained=(normalized_primary...,normalized_derived...)
    allunique(retained) || throw(ArgumentError("observation requests must be distinct"))
    return (retained=retained,displayed=selected)
end

function Commons.observation_requests(source::PropagationParameters,requests::Tuple;complete_pairs::Bool=false)
    selected=isempty(requests) ? (gamma,Zc,Yc,H,Tv,Ti,velocity) : requests
    normalized=_modal_component_requests(source,selected;complete_pairs)
    allunique(normalized) || throw(ArgumentError("observation requests must be distinct"))
    return (retained=normalized,displayed=selected)
end

function _modal_component_requests(source,selected;complete_pairs)
    normalized=Tuple[]
    complex_selectors=(gamma,Zc,Yc,H,Tv,Ti,_phase_Zc,_phase_Yc,_voltage_H,_current_H)
    for request in selected
        item=_normalize_modal_request(source,request)
        identity=Commons.normalize_observation_selector(request_identity(item))
        prefix=identity isa Tuple ? identity : (identity,)
        item=(prefix...,request_indices(item)...)
        if identity in complex_selectors
            append!(normalized,((identity,real,request_indices(item)...),
                (identity,imag,request_indices(item)...)))
        else
            push!(normalized,item)
        end
    end
    allunique(normalized) || throw(ArgumentError("observation requests must be distinct"))
    for item in copy(normalized)
        identity=request_identity(item)
        identity isa Tuple && length(identity)==2 && last(identity) in (real,imag,abs,angle) || continue
        selector=first(identity)
        selector in complex_selectors || continue
        other=last(identity)===real ? imag :
            last(identity)===imag ? real : last(identity)===abs ? angle : abs
        companion=(selector,other,request_indices(item)...)
        companion in normalized && continue
        complete_pairs || throw(ArgumentError("$(nameof(selector isa Base.Fix2 ? selector.f : selector)) requires a complete component pair"))
        push!(normalized,companion)
    end
    return Tuple(normalized)
end

function _modal_coordinates(source,request,values)
    identity=request_identity(request)
    selector=first(identity isa Tuple ? identity : (identity,))
    dims=size(values)
    indices=request_indices(request)
    selected=map(observation_indices,indices,dims)
    n=size(Tv(source),1)
    modes=string.(1:n)
    phase=get(details(source isa PropagationParameters ? source.parameters : source).data,
        :phase_coordinates,nothing)
    phase=phase===nothing ? string.(1:n) : string.(phase)
    mixed=selector in (Tv,Ti)
    phase_matrix=selector in (_phase_Zc,_phase_Yc,_voltage_H,_current_H)
    labels=mixed || phase_matrix ? phase : modes
    if length(dims)==2
        return (kind=:vector,indices,axis=:mode,axis_label="Mode",
            positions=selected[1],
            samples=selected[2],frequencies=copy(frequencies(source)[selected[2]]),
            frequency_unit=Units.units(:base,:hertz),extent=(n,length(frequencies(source))),
            labels,domain=:ModalDomain)
    end
    record=(kind=:matrix,indices,rows=selected[1],columns=selected[2],samples=selected[3],
        frequencies=copy(frequencies(source)[selected[3]]),
        frequency_unit=Units.units(:base,:hertz),extent=(n,n,length(frequencies(source))),
        labels,domain=mixed || phase_matrix ? :PhaseDomain : :ModalDomain)
    return mixed ? merge(record,(column_labels=modes,column_domain=:ModalDomain)) : record
end

function _modal_observation_quantity(source,request;unit=nothing,clip=true,atol=nothing,frequencies=nothing)
    frequencies===nothing || frequencies==LineCableModels.frequencies(source) ||
        throw(ArgumentError("derived modal observation uses stored frequencies"))
    atol===nothing || throw(ArgumentError("derived modal quantities have no engineering cutoff"))
    identity=request_identity(request)
    prefix=identity isa Tuple ? identity : (identity,)
    selector=first(prefix)
    transform=length(prefix)==2 ? prefix[2] : identity
    raw=observe(source,selector)
    indices=request_indices(request)
    original=getindex(raw,indices...)
    values=length(prefix)==2 ? transform.(original) : original
    available=Engine.resolution_available.(original) .& Engine.resolution_available.(values)
    origin=length(prefix)==2 && transform===angle ? iszero.(nominal.(original)) : false
    magnitude_origin=length(prefix)==2 && transform===abs ?
        iszero.(nominal.(original)) : false
    available=available .& .!origin
    reasons=broadcast(available,origin,magnitude_origin) do valid,at_zero,at_magnitude_origin
        valid ? nothing : at_zero ? :undefined_phase :
            at_magnitude_origin ? :undefined_first_order_magnitude : :nonfinite_value
    end
    resolved=broadcast((value,valid)->valid ? value : missing,values,available)
    q=Commons.request_quantity(request)
    native=Units.native_unit(q,basis(source))
    target=unit===nothing ? Units.display_unit(q,basis(source)) : unit
    T=typeof(float(real(nominal(zero(eltype(raw))))))
    factor=Units.scale_factor(native,target,T)
    assumptions=Engine.observation_assumptions(source,selector)
    undefined=any(x -> x===:undefined_first_order_magnitude,
        reasons isa AbstractArray ? reasons : (reasons,))
    components=undefined ? (nominal_magnitude=Commons.detach(abs.(nominal.(original))),
        real=Commons.detach(real.(original)),imaginary=Commons.detach(imag.(original)),
        unit=Units.native_unit(selector isa Base.Fix2 ? selector.f : selector,basis(source))) : nothing
    return (request,quantity=q,family=Symbol(nameof(selector isa Base.Fix2 ? selector.f : selector)),
        statistic=:value,values=Commons.detach(resolved,factor),unit=target,
        basis=basis(source),coordinates=_modal_coordinates(source,request,raw),
        assumptions,thresholds=nothing,available,engineering_zero=false,clipped=false,
        missing_reason=reasons,unavailable_components=components)
end

function Engine.observation_assumptions(
        source::Union{_ModalLineParameters,PropagationParameters},selector)
    retained=details(source).data
    upstream=get(retained,:selections,nothing)
    modal=get(retained,:modal,nothing)
    return upstream===nothing || modal===nothing ? nothing :
        (upstream=(Z=get(upstream,:Z,nothing),Y=get(upstream,:Y,nothing)),
         modal=get(modal,:effective,nothing),rotate=get(modal,:rotate,nothing))
end

function Commons.observation_quantity(source::_ModalLineParameters,request;kwargs...)
    identity=request_identity(request)
    selector=first(identity isa Tuple ? identity : (identity,))
    if selector in (Engine.Z,Engine.Y,Engine.R,Engine.X,Engine.L,Engine.G,Engine.B,Engine.C)
        return invoke(Commons.observation_quantity,
            Tuple{LineParameters,Any},
            source,request;kwargs...)
    end
    return _modal_observation_quantity(source,request;kwargs...)
end
Commons.observation_quantity(source::PropagationParameters,request;kwargs...) =
    _modal_observation_quantity(source,request;kwargs...)

function Commons.observation_gridpoint(source::PropagationParameters)
    original=Commons.observation_gridpoint(source.parameters)
    inputs=original.inputs===nothing ? (segment=(line_length=source.line_length,),) :
        merge(original.inputs,(segment=(line_length=source.line_length,),))
    return merge(original,(id=get(source.details.data,:gridpoint,nothing),
        source_gridpoint=get(source.details.data,:source_gridpoint,nothing),inputs,
        transformation=get(source.details.data,:modal,nothing)))
end
