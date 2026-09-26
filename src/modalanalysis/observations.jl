# Scientific units belong to the scalar quantity, independently of presentation.
Units.quantity(::typeof(gamma)) = Units.Quantity{:propagation_constant}()
Units.quantity(::typeof(alpha)) = Units.Quantity{:attenuation_constant}()
Units.quantity(::typeof(beta)) = Units.Quantity{:phase_constant}()
Units.quantity(::typeof(velocity)) = Units.Quantity{:phase_velocity}()
Units.quantity(::typeof(Zc)) = Units.Quantity{:characteristic_impedance}()
Units.quantity(::typeof(Yc)) = Units.Quantity{:characteristic_admittance}()
Units.quantity(::typeof(H)) = Units.Quantity{:forward_response}()
Units.quantity(::typeof(Tv)) = Units.Quantity{:voltage_basis}()
Units.quantity(::typeof(Ti)) = Units.Quantity{:current_basis}()
Units.quantity(value::Base.Fix2{typeof(H)}) = value.x.field===:voltage ?
    Units.Quantity{:voltage_propagation_function}() :
    Units.Quantity{:current_propagation_function}()
Units.quantity(value::Base.Fix2{typeof(Zc)}) = Units.Quantity{:phase_characteristic_impedance}()
Units.quantity(value::Base.Fix2{typeof(Yc)}) = Units.Quantity{:phase_characteristic_admittance}()

function Units.quantity(value::Base.Fix2{typeof(H)},transform::Function)
    name=value.x.field===:voltage ? :voltage_propagation_function :
        :current_propagation_function
    return _modal_bound_quantity(name,transform)
end
Units.quantity(::Base.Fix2{typeof(Zc)},transform::Function) =
    _modal_bound_quantity(:phase_characteristic_impedance,transform)
Units.quantity(::Base.Fix2{typeof(Yc)},transform::Function) =
    _modal_bound_quantity(:phase_characteristic_admittance,transform)

_modal_bound_quantity(name::Symbol,::typeof(abs)) = Units.Quantity{(name,:magnitude)}()
_modal_bound_quantity(name::Symbol,::typeof(angle)) = Units.Quantity{(name,:phase_angle)}()
_modal_bound_quantity(name::Symbol,::typeof(real)) = Units.Quantity{(name,:real)}()
_modal_bound_quantity(name::Symbol,::typeof(imag)) = Units.Quantity{(name,:imag)}()

for (selector, name) in ((gamma,:propagation_constant),(Zc,:characteristic_impedance),
        (Yc,:characteristic_admittance),(H,:forward_response),(Tv,:voltage_basis),
        (Ti,:current_basis))
    @eval begin
        Units.quantity(::typeof($selector),::typeof(abs)) = Units.Quantity{($(QuoteNode(name)),:magnitude)}()
        Units.quantity(::typeof($selector),::typeof(angle)) = Units.Quantity{($(QuoteNode(name)),:phase_angle)}()
        Units.native_unit(::Units.Quantity{($(QuoteNode(name)),:magnitude)}) = Units.native_unit(Units.Quantity{$(QuoteNode(name))}())
        Units.display_unit(::Units.Quantity{($(QuoteNode(name)),:magnitude)}) = Units.display_unit(Units.Quantity{$(QuoteNode(name))}())
        Units.native_unit(::Units.Quantity{($(QuoteNode(name)),:phase_angle)}) = Units.units(:base,:radian)
        Units.display_unit(::Units.Quantity{($(QuoteNode(name)),:phase_angle)}) = Units.units(:base,:degree)
        Units.label(::Units.Quantity{($(QuoteNode(name)),:magnitude)}) = Units.label(Units.Quantity{$(QuoteNode(name))}()) * " magnitude"
        Units.label(::Units.Quantity{($(QuoteNode(name)),:phase_angle)}) = Units.label(Units.Quantity{$(QuoteNode(name))}()) * " angle"
        Units.symbol(::Units.Quantity{($(QuoteNode(name)),:magnitude)}) = "|" * Units.symbol(Units.Quantity{$(QuoteNode(name))}()) * "|"
        Units.symbol(::Units.Quantity{($(QuoteNode(name)),:phase_angle)}) = "∠" * Units.symbol(Units.Quantity{$(QuoteNode(name))}())
    end
end
Units.quantity(::typeof(gamma),::typeof(real)) = Units.quantity(alpha)
Units.quantity(::typeof(gamma),::typeof(imag)) = Units.quantity(beta)
Units.quantity(::typeof(Zc),::typeof(real)) = Units.Quantity{:characteristic_resistance}()
Units.quantity(::typeof(Zc),::typeof(imag)) = Units.Quantity{:characteristic_reactance}()
Units.quantity(::typeof(Yc),::typeof(real)) = Units.Quantity{:characteristic_conductance}()
Units.quantity(::typeof(Yc),::typeof(imag)) = Units.Quantity{:characteristic_susceptance}()
Units.quantity(::typeof(H),::typeof(real)) = Units.Quantity{:forward_response_real}()
Units.quantity(::typeof(H),::typeof(imag)) = Units.Quantity{:forward_response_imag}()
Units.quantity(::typeof(Tv),::typeof(real)) = Units.Quantity{:voltage_basis_real}()
Units.quantity(::typeof(Tv),::typeof(imag)) = Units.Quantity{:voltage_basis_imag}()
Units.quantity(::typeof(Ti),::typeof(real)) = Units.Quantity{:current_basis_real}()
Units.quantity(::typeof(Ti),::typeof(imag)) = Units.Quantity{:current_basis_imag}()

Units.native_unit(::Units.Quantity{:attenuation_constant}) = Units.units(:base,:neper;per=(:base,:meter))
Units.display_unit(::Units.Quantity{:attenuation_constant}) = Units.units(:base,:neper;per=(:kilo,:meter))
Units.label(::Units.Quantity{:attenuation_constant}) = "Attenuation constant"
Units.symbol(::Units.Quantity{:attenuation_constant}) = "α"
Units.native_unit(::Units.Quantity{:phase_constant}) = Units.units(:base,:radian;per=(:base,:meter))
Units.display_unit(::Units.Quantity{:phase_constant}) = Units.units(:base,:radian;per=(:kilo,:meter))
Units.label(::Units.Quantity{:phase_constant}) = "Phase constant"
Units.symbol(::Units.Quantity{:phase_constant}) = "β"
Units.native_unit(::Units.Quantity{:phase_velocity}) = Units.units(:base,:meter;per=(:base,:second))
Units.display_unit(::Units.Quantity{:phase_velocity}) = Units.native_unit(Units.Quantity{:phase_velocity}())
Units.label(::Units.Quantity{:phase_velocity}) = "Phase velocity"
Units.symbol(::Units.Quantity{:phase_velocity}) = "vₚ"
for (name,base,scientific,sym) in (
        (:characteristic_resistance,:characteristic_impedance,"Characteristic resistance","R꜀"),
        (:characteristic_reactance,:characteristic_impedance,"Characteristic reactance","X꜀"),
        (:characteristic_conductance,:characteristic_admittance,"Characteristic conductance","G꜀"),
        (:characteristic_susceptance,:characteristic_admittance,"Characteristic susceptance","B꜀"),
        (:forward_response_real,:forward_response,"Real propagation function","Re(H)"),
        (:forward_response_imag,:forward_response,"Imaginary propagation function","Im(H)"),
        (:voltage_basis_real,:voltage_basis,"Real modal-to-phase voltage transformation","Re(Tv)"),
        (:voltage_basis_imag,:voltage_basis,"Imaginary modal-to-phase voltage transformation","Im(Tv)"),
        (:current_basis_real,:current_basis,"Real modal-to-phase current transformation","Re(Ti)"),
        (:current_basis_imag,:current_basis,"Imaginary modal-to-phase current transformation","Im(Ti)"))
    @eval begin
        Units.native_unit(::Units.Quantity{$(QuoteNode(name))}) = Units.native_unit(Units.Quantity{$(QuoteNode(base))}())
        Units.display_unit(::Units.Quantity{$(QuoteNode(name))}) = Units.display_unit(Units.Quantity{$(QuoteNode(base))}())
        Units.label(::Units.Quantity{$(QuoteNode(name))}) = $scientific
        Units.symbol(::Units.Quantity{$(QuoteNode(name))}) = $sym
    end
end
Units.native_unit(::Units.Quantity{:propagation_constant}) =
    Units.units(:base,:dimensionless;per=(:base,:meter))
Units.display_unit(::Units.Quantity{:propagation_constant}) =
    Units.units(:base,:dimensionless;per=(:kilo,:meter))
Units.native_unit(::Units.Quantity{:characteristic_impedance}) = Units.units(:base,:ohm)
Units.display_unit(::Units.Quantity{:characteristic_impedance}) = Units.units(:base,:ohm)
Units.native_unit(::Units.Quantity{:characteristic_admittance}) = Units.units(:base,:siemens)
Units.display_unit(::Units.Quantity{:characteristic_admittance}) = Units.units(:base,:siemens)
for name in (:forward_response,:voltage_basis,:current_basis)
    @eval begin
        Units.native_unit(::Units.Quantity{$(QuoteNode(name))}) = Units.units(:base,:dimensionless)
        Units.display_unit(::Units.Quantity{$(QuoteNode(name))}) = Units.units(:base,:dimensionless)
    end
end
Units.label(::Units.Quantity{:propagation_constant}) = "Propagation constant"
Units.label(::Units.Quantity{:characteristic_impedance}) = "Characteristic impedance"
Units.label(::Units.Quantity{:characteristic_admittance}) = "Characteristic admittance"
Units.label(::Units.Quantity{:forward_response}) = "Propagation function"
Units.label(::Units.Quantity{:voltage_basis}) = "Modal-to-phase voltage transformation"
Units.label(::Units.Quantity{:current_basis}) = "Modal-to-phase current transformation"
Units.symbol(::Units.Quantity{:propagation_constant}) = "γ"
Units.symbol(::Units.Quantity{:characteristic_impedance}) = "Zc"
Units.symbol(::Units.Quantity{:characteristic_admittance}) = "Yc"
Units.symbol(::Units.Quantity{:forward_response}) = "H"
Units.symbol(::Units.Quantity{:voltage_basis}) = "Tv"
Units.symbol(::Units.Quantity{:current_basis}) = "Ti"

for (name,parent,scientific,sym) in (
        (:phase_characteristic_impedance,:characteristic_impedance,"Phase-domain characteristic impedance","Zc,phase"),
        (:phase_characteristic_admittance,:characteristic_admittance,"Phase-domain characteristic admittance","Yc,phase"),
        (:voltage_propagation_function,:forward_response,"Voltage propagation function","Hᵥ"),
        (:current_propagation_function,:forward_response,"Current propagation function","Hᵢ"))
    @eval begin
        Units.native_unit(::Units.Quantity{$(QuoteNode(name))}) = Units.native_unit(Units.Quantity{$(QuoteNode(parent))}())
        Units.display_unit(::Units.Quantity{$(QuoteNode(name))}) = Units.display_unit(Units.Quantity{$(QuoteNode(parent))}())
        Units.label(::Units.Quantity{$(QuoteNode(name))}) = $scientific
        Units.symbol(::Units.Quantity{$(QuoteNode(name))}) = $sym
        Units.native_unit(::Units.Quantity{($(QuoteNode(name)),:magnitude)}) = Units.native_unit(Units.Quantity{$(QuoteNode(name))}())
        Units.display_unit(::Units.Quantity{($(QuoteNode(name)),:magnitude)}) = Units.display_unit(Units.Quantity{$(QuoteNode(name))}())
        Units.label(::Units.Quantity{($(QuoteNode(name)),:magnitude)}) = $scientific * " magnitude"
        Units.symbol(::Units.Quantity{($(QuoteNode(name)),:magnitude)}) = "|" * $sym * "|"
        Units.native_unit(::Units.Quantity{($(QuoteNode(name)),:phase_angle)}) = Units.units(:base,:radian)
        Units.display_unit(::Units.Quantity{($(QuoteNode(name)),:phase_angle)}) = Units.units(:base,:degree)
        Units.label(::Units.Quantity{($(QuoteNode(name)),:phase_angle)}) = $scientific * " angle"
        Units.symbol(::Units.Quantity{($(QuoteNode(name)),:phase_angle)}) = "∠" * $sym
        Units.native_unit(::Units.Quantity{($(QuoteNode(name)),:real)}) = Units.native_unit(Units.Quantity{$(QuoteNode(name))}())
        Units.display_unit(::Units.Quantity{($(QuoteNode(name)),:real)}) = Units.display_unit(Units.Quantity{$(QuoteNode(name))}())
        Units.label(::Units.Quantity{($(QuoteNode(name)),:real)}) = "Real " * lowercase($scientific)
        Units.symbol(::Units.Quantity{($(QuoteNode(name)),:real)}) = "Re(" * $sym * ")"
        Units.native_unit(::Units.Quantity{($(QuoteNode(name)),:imag)}) = Units.native_unit(Units.Quantity{$(QuoteNode(name))}())
        Units.display_unit(::Units.Quantity{($(QuoteNode(name)),:imag)}) = Units.display_unit(Units.Quantity{$(QuoteNode(name))}())
        Units.label(::Units.Quantity{($(QuoteNode(name)),:imag)}) = "Imaginary " * lowercase($scientific)
        Units.symbol(::Units.Quantity{($(QuoteNode(name)),:imag)}) = "Im(" * $sym * ")"
    end
end

const _phase_Zc = Base.Fix2(Zc,(domain=PhaseDomain,))
const _phase_Yc = Base.Fix2(Yc,(domain=PhaseDomain,))
const _voltage_H = Base.Fix2(H,(domain=PhaseDomain,field=:voltage))
const _current_H = Base.Fix2(H,(domain=PhaseDomain,field=:current))

# Finite representation bindings carry only domain and, for H, field.
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

Grammar.normalize_observation_selector(::typeof(alpha)) = (gamma,real)
Grammar.normalize_observation_selector(::typeof(beta)) = (gamma,imag)
Grammar.normalize_observation_selector(::Val{:alpha}) = Grammar.normalize_observation_selector(alpha)
Grammar.normalize_observation_selector(::Val{:beta}) = Grammar.normalize_observation_selector(beta)

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

function Grammar.observation_requests(source::_ModalLineParameters,requests::Tuple;complete_pairs::Bool=false)
    selected=isempty(requests) ? (Engine.Z,Engine.Y,gamma,Zc,Yc,Tv,Ti) : requests
    primary=Tuple(filter(request -> begin
        identity=request_identity(request)
        first(identity isa Tuple ? identity : (identity,)) in
            (Engine.Z,Engine.Y,Engine.R,Engine.X,Engine.L,Engine.G,Engine.B,Engine.C,real,imag,abs,angle)
    end,selected))
    derived=Tuple(filter(request -> request ∉ primary,selected))
    normalized_primary=isempty(primary) ? () :
        invoke(Grammar.observation_requests,
            Tuple{LineParameters,Tuple},
            source,primary;complete_pairs).retained
    normalized_derived=_modal_component_requests(source,derived;complete_pairs)
    retained=(normalized_primary...,normalized_derived...)
    allunique(retained) || throw(ArgumentError("observation requests must be distinct"))
    return (retained=retained,displayed=selected)
end

function Grammar.observation_requests(source::PropagationParameters,requests::Tuple;complete_pairs::Bool=false)
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
        identity=Grammar.normalize_observation_selector(request_identity(item))
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
    q=Grammar.request_quantity(request)
    native=Units.native_unit(q,basis(source))
    target=unit===nothing ? Units.display_unit(q,basis(source)) : unit
    T=typeof(float(real(nominal(zero(eltype(raw))))))
    factor=Units.scale_factor(native,target,T)
    assumptions=Engine.observation_assumptions(source,selector)
    undefined=any(x -> x===:undefined_first_order_magnitude,
        reasons isa AbstractArray ? reasons : (reasons,))
    components=undefined ? (nominal_magnitude=Grammar.detach(abs.(nominal.(original))),
        real=Grammar.detach(real.(original)),imaginary=Grammar.detach(imag.(original)),
        unit=Units.native_unit(selector isa Base.Fix2 ? selector.f : selector,basis(source))) : nothing
    return (request,quantity=q,family=Symbol(nameof(selector isa Base.Fix2 ? selector.f : selector)),
        statistic=:value,values=Grammar.detach(resolved,factor),unit=target,
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
         modal=get(modal,:effective,nothing))
end

function Grammar.observation_quantity(source::_ModalLineParameters,request;kwargs...)
    identity=request_identity(request)
    selector=first(identity isa Tuple ? identity : (identity,))
    if selector in (Engine.Z,Engine.Y,Engine.R,Engine.X,Engine.L,Engine.G,Engine.B,Engine.C)
        return invoke(Grammar.observation_quantity,
            Tuple{LineParameters,Any},
            source,request;kwargs...)
    end
    return _modal_observation_quantity(source,request;kwargs...)
end
Grammar.observation_quantity(source::PropagationParameters,request;kwargs...) =
    _modal_observation_quantity(source,request;kwargs...)

function Grammar.observation_gridpoint(source::PropagationParameters)
    original=Grammar.observation_gridpoint(source.parameters)
    inputs=original.inputs===nothing ? (segment=(line_length=source.line_length,),) :
        merge(original.inputs,(segment=(line_length=source.line_length,),))
    return merge(original,(id=get(source.details.data,:gridpoint,nothing),
        source_gridpoint=get(source.details.data,:source_gridpoint,nothing),inputs,
        transformation=get(source.details.data,:modal,nothing)))
end
