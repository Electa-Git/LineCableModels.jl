const _ObservedLineSource = Union{LineParameters,SeriesImpedance,ShuntAdmittance}
_line_families(::LineParameters) = (Z,Y)
_line_families(::SeriesImpedance) = (Z,)
_line_families(::ShuntAdmittance) = (Y,)
_primary_family(selector) = selector in (Z,R,X,L) ? Z : selector in (Y,G,B,C) ? Y : nothing
observation_assumptions(::Union{SeriesImpedance,ShuntAdmittance},selector) = nothing
function observation_assumptions(source::AbstractCoreResult,selector)
    family=_primary_family(selector)
    family===nothing && return nothing
    return get(get(details(source).data,:selections,(;)),Symbol(nameof(family)),nothing)
end
_primary_key(identity) = identity isa Function ? nameof(identity) :
    first(identity) in (Z,Y) && any(in((abs,angle)),Base.tail(identity)) ?
    Symbol(nameof(first(identity)),last(filter(in((abs,angle)),Base.tail(identity)))===abs ? :_abs : :_angle) : nameof(first(identity))
_pair_keys(::typeof(Z)) = ((:R,:X),(:Z_abs,:Z_angle),(:R,:L))
_pair_keys(::typeof(Y)) = ((:G,:B),(:Y_abs,:Y_angle),(:G,:C))
_pair_selector(key::Symbol) = key===:Z_abs ? (Z,abs) : key===:Z_angle ? (Z,angle) :
    key===:Y_abs ? (Y,abs) : key===:Y_angle ? (Y,angle) : (getfield(@__MODULE__,key),)

function _line_request(source,request)
    identity=request_identity(request)
    prefix=identity isa Tuple ? identity : (identity,)
    indices=request_indices(request)
    first(prefix) in (real,imag,abs,angle) && throw(ArgumentError("use a physical quantity or a family-qualified transform"))
    if length(prefix)>=2 && prefix[2] in (real,imag)
        component=first(prefix)===Z ? (prefix[2]===real ? R : X) :
            first(prefix)===Y ? (prefix[2]===real ? G : B) :
            throw(ArgumentError("real/imag transforms require Z or Y"))
        prefix=(component,prefix[3:end]...)
    end
    family=_primary_family(first(prefix))
    family in _line_families(source) || throw(ArgumentError("quantity does not belong to an available primary family"))
    diagonal=diag in prefix
    rank=diagonal ? 2 : 3
    isempty(indices) && (indices=ntuple(_ -> Colon(),rank))
    !diagonal && length(indices)==2 && (indices=(indices...,Colon()))
    length(indices)==rank || throw(DimensionMismatch("a primary request requires $rank coordinate selectors"))
    all(x -> x isa Colon || x isa Integer && !(x isa Bool) || x isa AbstractVector{<:Integer},indices) ||
        throw(ArgumentError("coordinates must be integer indices, ranges, vectors, or :"))
    # Polar diagonal access is obtained from the original complex diagonal.
    supported=Tuple(filter(!=(diag),prefix))
    identity_checked=length(supported)==1 ? only(supported) : supported
    identity_checked in observables(typeof(source)) || throw(ArgumentError("unsupported line representation"))
    return (prefix...,indices...)
end

function _line_observation_requests(source::_ObservedLineSource,requests::Tuple;complete_pairs::Bool=false)
    selected=isempty(requests) ? Tuple(_line_families(source)) : requests
    expanded=Tuple[]
    for item in selected
        if item in (real,imag,abs,angle)
            append!(expanded,[_line_request(source,(family,item)) for family in _line_families(source)])
            continue
        end
        normalized=_line_request(source,item)
        identity=request_identity(normalized)
        prefix=identity isa Tuple ? identity : (identity,)
        if first(prefix) in (Z,Y) && all(==(diag),Base.tail(prefix))
            for key in first(_pair_keys(first(prefix)))
                push!(expanded,(_pair_selector(key)...,Base.tail(prefix)...,request_indices(normalized)...))
            end
        else
            push!(expanded,normalized)
        end
    end
    allunique(expanded) || throw(ArgumentError("observation requests must be distinct"))
    retained=Tuple[]
    for family in _line_families(source)
        entries=filter(request -> _primary_family(first(request))===family,expanded)
        if isempty(entries)
            template=first(expanded)
            diagonal=diag in (request_identity(template) isa Tuple ? request_identity(template) : ())
            entries=[(_pair_selector(key)...,(diagonal ? (diag,) : ())...,request_indices(template)...)
                for key in first(_pair_keys(family))]
        end
        keys=Tuple(_primary_key(request_identity(request)) for request in entries)
        choices=filter(pair -> all(in(pair),keys),_pair_keys(family))
        if length(entries)==1 && complete_pairs
            chosen=first(choices)
            diagonal=diag in (request_identity(first(entries)) isa Tuple ? request_identity(first(entries)) : ())
            entries=[(_pair_selector(key)...,(diagonal ? (diag,) : ())...,request_indices(first(entries))...)
                for key in chosen]
        elseif length(entries)!=2 || length(unique(keys))!=2 || length(choices)!=1
            throw(ArgumentError("$(nameof(family)) requires one complete pair: $(_pair_keys(family)); received $keys"))
        end
        dimensions=size(source isa LineParameters ? source.Z : source)
        coordinate(request) = begin
            identity=request_identity(request)
            diagonal=identity isa Tuple && diag in identity
            dims=diagonal ? (dimensions[1],dimensions[3]) : dimensions
            map(observation_indices,request_indices(request),dims)
        end
        coordinate(first(entries))==coordinate(last(entries)) ||
            throw(DimensionMismatch("both quantities of a primary pair must select identical coordinates and rank"))
        chosen=only(filter(pair -> Set(pair)==Set(_primary_key(request_identity(q)) for q in entries),_pair_keys(family)))
        append!(retained,[only(filter(q -> _primary_key(request_identity(q))===key,entries)) for key in chosen])
    end
    return (retained=Tuple(retained),displayed=Tuple(expanded))
end
Grammar.observation_requests(source::LineParameters,requests::Tuple;complete_pairs::Bool=false) =
    _line_observation_requests(source,requests;complete_pairs)
Grammar.observation_requests(source::Union{SeriesImpedance,ShuntAdmittance},requests::Tuple;complete_pairs::Bool=false) =
    _line_observation_requests(source,requests;complete_pairs)

"""
$(TYPEDSIGNATURES)

Select physical matrix coordinates for an indexed line-parameter observation.

# Arguments

- `source`: Line parameters, series impedance, or shunt admittance tensor.
- `request`: Normalized observable request with row/column/sample indices, or
  diagonal/sample indices for a diagonal request.
- `frequencies`: Supplied frequency samples \\[Hz\\], or `nothing`. Line parameters
  use their stored frequencies and reject a differing supplied vector.

# Returns

- A named tuple retaining the original indices, selected row/column/sample
  positions, frequencies \\[Hz\\], coordinate labels, full tensor extent, domain,
  and `:matrix` or `:diagonal` representation. Selection order is preserved.
"""
function line_coordinates(source,request,frequencies)
    identity=request_identity(request)
    diagonal=identity isa Tuple && diag in identity
    dimensions=size(source isa LineParameters ? source.Z : source)
    indices=request_indices(request)
    dims=diagonal ? (dimensions[1],dimensions[3]) : dimensions
    selected=map(observation_indices,indices,dims)
    rows=selected[1]
    columns=diagonal ? copy(rows) : selected[2]
    samples=last(selected)
    f=_resolution_frequencies(source,frequencies)
    f===nothing || length(f)==dimensions[3] || throw(DimensionMismatch("frequency count differs from tensor depth"))
    labels=source isa LineParameters ? get(details(source).data,:coordinates,nothing) : nothing
    labels=labels===nothing ? string.(1:dimensions[1]) : copy(labels)
    return (kind=diagonal ? :diagonal : :matrix,indices,rows,columns,samples,
        frequencies=f===nothing ? nothing : copy(f[samples]),frequency_unit=Units.units(:base,:hertz),
        extent=dimensions,labels,domain=source isa LineParameters ? nameof(domain(source)) : :unspecified)
end

_scaled_thresholds(value::Nothing,factor) = nothing
_scaled_thresholds(value::NamedTuple,factor) = map(x -> _scaled_thresholds(x,factor),value)
_scaled_thresholds(value,factor) = Grammar.detach(value,factor)

function _line_observation_quantity(source::_ObservedLineSource,request;
        unit=nothing,clip=true,atol=nothing,frequencies=nothing)
    identity=request_identity(request)
    prefix=identity isa Tuple ? identity : (identity,)
    selector=first(prefix)
    indices=request_indices(request)
    coordinates=line_coordinates(source,request,frequencies)
    f=_resolution_frequencies(source,frequencies)
    sampled_f=f===nothing ? nothing : f[last(indices)]
    polar=any(in((abs,angle)),prefix)
    phase=angle in prefix
    diagonal=diag in prefix
    original=polar ? (diagonal ? observe(source,selector,diag,indices...) : observe(source,selector,indices...)) :
        _line_observation_values(source,request;frequencies=f)
    resolution=observation_resolution(original,selector;atol,frequencies=sampled_f,
        result_basis=basis(source),line_length=line_length(source))
    values=polar ? (phase ? angle.(original) : abs.(original)) : original
    available=resolution.available .& resolution_available.(values)
    exact_origin=polar ? iszero.(nominal.(original)) : false
    phase && (available=available .& .!exact_origin)
    mask=resolution.unresolved===nothing ? false : resolution.unresolved
    resolved=broadcast(values,mask,available) do value,unresolved,valid
        !valid && return missing
        clip && phase && unresolved && return missing
        clip && unresolved ? value-nominal(value) : value
    end
    reasons=broadcast(available,mask,exact_origin) do valid,unresolved,origin
        !valid ? (polar && !phase && origin ? :undefined_first_order_magnitude : phase && origin ? :undefined_phase :
            selector in (L,C) ? :unavailable_proxy : :nonfinite_value) :
            clip && phase && unresolved ? :engineering_zero_phase : nothing
    end
    q=Grammar.request_quantity(request)
    native=Units.native_unit(q,basis(source))
    target=unit===nothing ? Units.display_unit(q,basis(source)) : unit
    T=typeof(float(real(nominal(zero(Base.nonmissingtype(eltype(original)))))))
    factor=Units.scale_factor(native,target,T)
    threshold_unit=resolution.unit
    threshold_factor=phase ? one(T) : factor
    undefined=polar && !phase && any(x -> x===:undefined_first_order_magnitude,
        reasons isa AbstractArray ? reasons : (reasons,))
    # Components are retained only for the undefined first-order polar value;
    # they are scientific values, never a reference back to the source tensor.
    components=undefined ? (nominal_magnitude=Grammar.detach(abs.(nominal.(original))),real=Grammar.detach(real.(original)),imaginary=Grammar.detach(imag.(original)),
        unit=Units.native_unit(selector,basis(source))) : nothing
    return (request,quantity=q,family=Symbol(nameof(_primary_family(selector))),statistic=:value,
        values=Grammar.detach(resolved isa AbstractArray && ndims(resolved)==0 ? only(resolved) : resolved,factor),unit=target,basis=basis(source),coordinates,
        assumptions=observation_assumptions(source,selector),
        thresholds=(kind=resolution.kind,values=_scaled_thresholds(resolution.atol,threshold_factor),
            unit=phase ? threshold_unit : target),available,engineering_zero=mask,clipped=clip,
        missing_reason=reasons,unavailable_components=components)
end
Grammar.observation_quantity(source::LineParameters,request;kwargs...) =
    _line_observation_quantity(source,request;kwargs...)
Grammar.observation_quantity(source::Union{SeriesImpedance,ShuntAdmittance},request;kwargs...) =
    _line_observation_quantity(source,request;kwargs...)

function Grammar.observation_requests(source::CableConstants,requests::Tuple;complete_pairs::Bool=false)
    selected=isempty(requests) ? (R,L,G,C) : requests
    Grammar.validate_observables(source,selected,())
    return (retained=selected,displayed=selected)
end

function Grammar.observation_quantity(source::CableConstants,request;
        unit=nothing,clip=true,atol=nothing,frequencies=nothing)
    frequencies===nothing || frequencies==[source.frequency] || throw(ArgumentError("frequency differs from cable constants"))
    indices=request_indices(request)
    length(indices)<=1 || throw(ArgumentError("cable constants select assembly indices"))
    index=isempty(indices) ? Colon() : only(indices)
    selected=observation_indices(index,length(source))
    selector=request_identity(request)
    values=getindex(observe(source,selector),index)
    resolution=observation_resolution(values,selector;atol,result_basis=:pul)
    resolved = if !clip
        values
    elseif resolution.unresolved === nothing
        resolution.available === false ? missing : values
    else
        broadcast(values, resolution.unresolved, resolution.available) do value, unresolved, available
            available || return missing
            return unresolved ? value - nominal(value) : value
        end
    end
    q=Units.quantity(selector)
    target=unit===nothing ? Units.display_unit(q,:pul) : unit
    T=typeof(float(nominal(zero(eltype(values)))))
    factor=Units.scale_factor(Units.native_unit(q,:pul),target,T)
    return (request,quantity=q,family=:constants,statistic=:value,values=Grammar.detach(resolved,factor),
        unit=target,basis=:pul,coordinates=(kind=:assemblies,indices=(index,),assemblies=selected,
            labels=copy(source.cores),frequencies=[source.frequency],frequency_unit=Units.units(:base,:hertz),extent=(length(source),1)),
        thresholds=(kind=resolution.kind,values=_scaled_thresholds(resolution.atol,factor),unit=target),
        available=resolution.available,engineering_zero=resolution.unresolved,clipped=clip,missing_reason=nothing)
end
