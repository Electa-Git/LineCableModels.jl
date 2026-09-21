"""
$(TYPEDEF)

Retain detached scientific products for one completed gridpoint. Collections
are ordinary vectors of these objects. No result source, parent collection,
construction closure, or flattened table is retained.

$(TYPEDFIELDS)
"""
struct ObservedResult
    "Original identity, physical inputs, selections, and uncertainty interpretation."
    gridpoint::NamedTuple
    "Requested numerical products, coordinates, units, and availability records."
    quantities::Vector{NamedTuple}
    "Completed comparisons with candidate and separate reference identities."
    errors::Vector{NamedTuple}
    "Completed execution and performance measurements in their original scopes."
    timings::NamedTuple

    function ObservedResult(gridpoint::NamedTuple,quantities::AbstractVector,
            errors::AbstractVector,timings::NamedTuple)
        foreach(_validate_observed_quantity,quantities)
        allunique(q.request for q in quantities) || throw(ArgumentError("retained requests must be distinct"))
        _observed_fields(gridpoint,(:id,),"gridpoint")
        id=gridpoint.id
        id===nothing || _validate_observed_id(id)
        haskey(timings,:candidate_id) && timings.candidate_id!=id &&
            throw(ArgumentError("recorded timings must identify this candidate"))
        all(row -> row isa NamedTuple && haskey(row,:candidate_id) &&
            haskey(row,:reference_id) && id!==nothing && row.candidate_id==id,errors) ||
            throw(ArgumentError("completed comparisons must identify this candidate and a separate reference"))
        foreach(_validate_observed_comparison,errors)
        allunique((row.reference_id,row.request,row.band,row.normalization) for row in errors) ||
            throw(ArgumentError("completed comparisons must be distinct"))
        return new(detach(gridpoint),NamedTuple[detach(q) for q in quantities],
            NamedTuple[detach(row) for row in errors],detach(timings))
    end
end

"""
$(TYPEDSIGNATURES)

Normalize retained requests at the observation boundary. Concrete scientific
owners enforce their complete representations. `complete_pairs=true` is used
by raw display conveniences; an observed-input consumer only selects retained
products. The return record separates `retained` from `displayed` requests.
"""
function observation_requests(source,requests::Tuple;complete_pairs::Bool=false)
    retained=isempty(requests) ? observables(typeof(source)) : requests
    validate_observables(source,retained,())
    return (;retained,displayed=retained)
end

function observation_requests(source::ObservedResult,requests::Tuple;complete_pairs::Bool=false)
    identities=request_identity.(getproperty.(source.quantities,:request))
    selected=isempty(requests) ? Tuple(count(==(identity),identities)==1 ? identity : product.request
        for (identity,product) in zip(identities,source.quantities)) : requests
    displayed=Any[]
    for request in selected
        identity=request_identity(request)
        family=identity isa Function ? nameof(identity) : nothing
        products=filter(product -> get(product,:family,nothing)===family &&
            get(product,:statistic,nothing)===:value,source.quantities)
        if isempty(products)
            observation_product(source,request)
            push!(displayed,request)
        else
            for product in products
                indices=request_indices(request)
                prefix=request_identity(product.request)
                selector=prefix isa Tuple ? prefix : (prefix,)
                selection=isempty(indices) ? prefix : (selector...,indices...)
                observation_product(source,selection)
                push!(displayed,selection)
            end
        end
    end
    result=Tuple(displayed)
    return (retained=result,displayed=result)
end

"""
$(TYPEDSIGNATURES)

Extract one owned numerical product with coordinates and units. Scientific
owners extend this operation to describe their tensor coordinates and apply
their numerical resolution. The fallback provides ordinary array coordinates.
"""
function observation_quantity(source,request;unit=nothing,clip=true,atol=nothing,frequencies=nothing)
    values=_observe_request(source,request)
    q=request_quantity(request)
    native=native_unit(q,basis(source))
    displayed=unit===nothing ? display_unit(q,basis(source)) : unit
    coordinate=(kind=:array,indices=request_indices(request),
        extent=values isa AbstractArray ? size(values) : (),)
    return (request,quantity=q,family=:other,statistic=:value,
        values=detach(values,scale_factor(native,displayed)),unit=displayed,
        basis=basis(source),coordinates=coordinate,thresholds=nothing,
        available=nothing,engineering_zero=nothing,clipped=false,missing_reason=nothing)
end

"""
$(TYPEDSIGNATURES)

Construct one detached observation from a completed primary result and already
completed comparison/timing records. This action never computes comparisons,
collects timings, or evaluates a problem.

# Keywords

- `comparisons=()`: Completed records identifying this candidate and its reference.
- `timings=(;)`: Recorded measurements; absent evidence stays absent.
- `gridpoint=nothing`: Explicit description for an external numerical source.
- `clip=true`, `atol=nothing`: Engineering recentering and native-unit cutoffs.
- `length_unit=:kilo`, `frequency_unit=:base`: Display-unit prefixes.
- `units=()`, `quantity_units=nothing`: Aligned or quantity-keyed unit overrides.
- `frequencies=nothing`: Frequency context \\[Hz\\] for standalone numerical tensors.
- `complete_pairs=false`: Complete a raw display selection through the shared normalizer.
"""
function ObservedResult(source,requests::Tuple=();comparisons=(),timings=(;),gridpoint=nothing,
        clip::Bool=true,atol=nothing,units::Tuple=(),length_unit::Symbol=:kilo,
        frequency_unit::Symbol=:base,quantity_units=nothing,frequencies=nothing,
        complete_pairs::Bool=false)
    selected=observation_requests(source,requests;complete_pairs).retained
    isempty(units) || length(units)==length(selected) || throw(DimensionMismatch(
        "one display unit is required for each retained quantity"))
    isempty(units) || quantity_units===nothing || throw(ArgumentError("use units or quantity_units, not both"))
    atol isa Real && length(unique(request_quantity.(selected)))>1 && throw(ArgumentError(
        "multiple quantities require component-keyed native-unit cutoffs"))
    description=gridpoint===nothing ? observation_gridpoint(source) : gridpoint
    quantities=map(eachindex(selected)) do index
        request=selected[index]
        q=request_quantity(request)
        unit=isempty(units) ? q isa Units.Quantity{:frequency} ? Units.units(frequency_unit,:hertz) :
            display_unit(q,basis(source),_unit_override(quantity_units,request);length_prefix=length_unit) :
            display_unit(q,basis(source),units[index];length_prefix=length_unit)
        record=observation_quantity(source,request;unit,clip,atol,frequencies)
        if haskey(record.coordinates,:frequencies) && record.coordinates.frequencies!==nothing
            frequency_target=Units.units(frequency_unit,:hertz)
            f=record.coordinates.frequencies
            factor=isempty(f) ? 1 : scale_factor(Units.units(:base,:hertz),frequency_target,typeof(float(nominal(first(f)))))
            record=merge(record,(coordinates=merge(record.coordinates,
                (frequencies=detach(f,factor),frequency_unit=frequency_target)),))
        end
        record
    end
    return ObservedResult(description,quantities,collect(comparisons),timings)
end

function observation_request(observed::ObservedResult,request)
    identity=request_identity(request)
    any(q -> request_identity(q.request)==identity,observed.quantities) ||
        throw(ArgumentError("the requested quantity was not retained"))
    return (;identity,quantity=request_quantity(request),indices=request_indices(request))
end

"""
$(TYPEDSIGNATURES)

Read a retained quantity in its recorded unit. Selection uses original
coordinates; no source extraction, clipping, or absent-quantity derivation occurs.
"""
function observe(observed::ObservedResult,selectors...)
    request=length(selectors)==1 ? only(selectors) : selectors
    return detach(observation_product(observed,request).values)
end

function basis(observed::ObservedResult)
    bases=unique(q.basis for q in observed.quantities)
    length(bases)==1 || throw(ArgumentError("observation does not have one quantity basis"))
    return only(bases)
end

observation_gridpoint(observed::ObservedResult) = detach(observed.gridpoint)
function Base.show(io::IO,observed::ObservedResult)
    print(io,"ObservedResult(",length(observed.quantities)," quantities, ",length(observed.errors)," comparisons)")
end
Base.show(io::IO,::MIME"text/plain",observed::ObservedResult) = show(io,observed)
Base.summary(io::IO,observed::ObservedResult) = show(io,observed)

"""
$(TYPEDSIGNATURES)

Lift atomic observation over an ordinary collection. Completed comparisons and
point timings are joined by original identities before detachment, so filtering
and reordering cannot associate a candidate with another point's evidence.
"""
function observables(sources::Union{AbstractVector,Tuple,AbstractResultSpace},requests::Tuple=();
        comparisons=nothing,timings=nothing,kwargs...)
    return map(collect(sources)) do source
        id=get(observation_gridpoint(source),:id,nothing)
        errors=comparisons===nothing ? (source isa ObservedResult ? source.errors : ()) :
            filter(record -> record.candidate_id==id,comparisons)
        recorded=timings===nothing ? (source isa ObservedResult ? source.timings : (;)) : timings isa NamedTuple ? timings : begin
            matches=filter(record -> record.candidate_id==id,timings)
            length(matches)<=1 || throw(ArgumentError("multiple timing records for one candidate identity"))
            isempty(matches) ? (;) : only(matches)
        end
        ObservedResult(source,requests;comparisons=errors,timings=recorded,kwargs...)
    end
end

"""
$(TYPEDSIGNATURES)

Select a retained quantity record and optional original coordinates. Missing or
ambiguous requests fail; this operation never extracts or derives a quantity.
"""
function observation_product(observed::ObservedResult,request;unit=nothing,frequency_unit=nothing)
    matches=filter(q -> request_identity(q.request)==request_identity(request),observed.quantities)
    length(matches)>1 && (matches=filter(q -> q.request==request,matches))
    length(matches)==1 || throw(ArgumentError("requested retained product is absent or ambiguous"))
    product=only(matches)
    indices=request_indices(request)
    isempty(indices) || request==product.request || (product=_selected_product(product,request,indices))
    return _reexpress_product(product;unit,frequency_unit)
end

function _selected_product(product,request,indices)
    c=product.coordinates
    c.kind in (:matrix,:diagonal,:assemblies,:samples) || throw(ArgumentError("select this retained product by its complete request"))
    matrix=haskey(c,:rows)
    dimensions=matrix ? c.kind===:diagonal ? (c.rows,c.samples) : (c.rows,c.columns,c.samples) : (c.assemblies,)
    if c.kind===:samples
        dimensions=(dimensions...,c.trials)
        length(indices)==length(dimensions)-1 && (indices=(indices...,Colon()))
    end
    length(indices)==length(dimensions) || throw(DimensionMismatch("request rank differs from retained coordinates"))
    positions=map(indices,dimensions) do requested,retained
        wanted=requested isa Colon ? retained : requested isa Integer ? [requested] : collect(requested)
        selected=map(wanted) do index
            position=findfirst(==(index),retained)
            position===nothing && throw(ArgumentError("coordinate $index was not retained"))
            position
        end
        requested isa Integer ? only(selected) : selected
    end
    selected_dimensions=map(dimensions,positions) do dimension,position
        position isa Integer ? [dimension[position]] : dimension[position]
    end
    select_values(value)=value isa AbstractArray ? reshape(value,length.(dimensions)...)[positions...] : value
    sample_axis=matrix ? (c.kind===:diagonal ? 2 : 3) : nothing
    f=if c.frequencies===nothing || sample_axis===nothing
        c.frequencies
    else
        selected_samples=positions[sample_axis]
        c.frequencies[selected_samples isa Integer ? [selected_samples] : selected_samples]
    end
    coordinate=matrix ? merge(c,(indices,rows=first(selected_dimensions),
        columns=c.kind===:diagonal ? first(selected_dimensions) : selected_dimensions[2],
        samples=selected_dimensions[sample_axis],frequencies=f)) :
        merge(c,(indices,assemblies=first(selected_dimensions)))
    c.kind===:samples && (coordinate=merge(coordinate,(trials=last(selected_dimensions),)))
    thresholds=product.thresholds
    if thresholds!==nothing && sample_axis!==nothing
        select_cutoff(cutoff::AbstractVector)=cutoff[positions[sample_axis]]
        select_cutoff(cutoff::NamedTuple)=map(select_cutoff,cutoff)
        select_cutoff(cutoff)=cutoff
        thresholds=merge(thresholds,(values=select_cutoff(thresholds.values),))
    end
    components=get(product,:unavailable_components,nothing)
    components===nothing || (components=merge(components,(nominal_magnitude=select_values(components.nominal_magnitude),
        real=select_values(components.real),imaginary=select_values(components.imaginary))))
    result=merge(product,(request,values=select_values(product.values),coordinates=coordinate,
        available=select_values(product.available),engineering_zero=select_values(product.engineering_zero),
        missing_reason=select_values(product.missing_reason),thresholds))
    return haskey(product,:unavailable_components) ? merge(result,(unavailable_components=components,)) : result
end

# Numeric equality also checks the uncertainty graph. For a dependency-aware
# number, equal marginal deviations alone do not make the difference certain.
_same_observed_number(a::Number,b::Number) = isequal(nominal(a),nominal(b)) &&
    isequal(uncertainty(a),uncertainty(b)) && iszero(uncertainty(a-b))
_same_observed_number(a,b) = isequal(a,b)
_same_observed_values(a::AbstractArray,b::AbstractArray) = size(a)==size(b) && all(_same_observed_number.(a,b))
_same_observed_values(a::Tuple,b::Tuple) = length(a)==length(b) && all(_same_observed_values(x,y) for (x,y) in zip(a,b))
_same_observed_values(a::NamedTuple,b::NamedTuple) = keys(a)==keys(b) && all(_same_observed_values(x,y) for (x,y) in zip(values(a),values(b)))
_same_observed_values(a,b) = _same_observed_number(a,b)

_same_observed_dependencies(a::Number,b::Number) = iszero(uncertainty(a-b))
_same_observed_dependencies(a::AbstractArray,b::AbstractArray) = size(a)==size(b) && all(_same_observed_dependencies.(a,b))
_same_observed_dependencies(a::NamedTuple,b::NamedTuple) = keys(a)==keys(b) && all(_same_observed_dependencies(x,y) for (x,y) in zip(values(a),values(b)))
_same_observed_dependencies(a,b) = true

"""
$(TYPEDSIGNATURES)

Return display groups and every original member identity for one retained
request. Eligibility is established by physical-point identity, owner-selected
formulations and controls, statistical interpretation, coordinates, and applied
cutoffs. Exact numerical and uncertainty-dependency agreement only verifies an
already established equivalence. Unidentified assumptions remain separate.
"""
function observation_groups(observed;request,band=nothing,normalization=nothing,reference=nothing)
    points=observed isa ObservedResult ? [observed] : observed
    groups=NamedTuple[]
    keys=Any[]
    products=Any[]
    for (index,point) in enumerate(points)
        product=if band===nothing
            observation_product(point,request)
        else
            matches=filter(row -> row.request==request && isequal(row.band,band) &&
                row.normalization==normalization && (reference===nothing || row.reference_id==reference),point.errors)
            length(matches)==1 || throw(ArgumentError("requested completed comparison is absent or ambiguous"))
            row=only(matches)
            (quantity=row.quantity,statistic=row.statistic,coordinates=row.coordinates,
                basis=get(row.settings,:basis,nothing),unit=row.absolute_unit,
                thresholds=(reference=get(row.settings,:atol,nothing),candidate=get(row.settings,:candidate_atol,nothing)),
                available=get(row.settings,:status,nothing),engineering_zero=nothing,
                assumptions=get(row,:assumptions,nothing),values=(absolute=row.absolute,relative=row.relative),
                interpretation=get(row.settings,:estimators,nothing))
        end
        id=get(point.gridpoint,:id,nothing)
        assumptions=get(product,:assumptions,nothing)
        physical=id===nothing ? nothing : (id.source_id,id.problem_index)
        key=(physical,assumptions,product.quantity,product.statistic,
            get(point.gridpoint,:uncertainty,nothing),product.coordinates,product.basis,
            product.unit,product.thresholds,
            band,normalization,reference,get(product,:interpretation,nothing))
        semantic=physical===nothing || assumptions===nothing || ismissing(assumptions) ? Int[] :
            findall(i -> _same_observed_values(keys[i],key),eachindex(keys))
        # The uncertainty interpretation includes the actual dependency graph.
        # Equal scalar means alone do not establish that interpretation.
        same_interpretation=filter(i -> _same_observed_dependencies(products[i].values,product.values),semantic)
        for i in same_interpretation
            _same_observed_values(products[i].values,product.values) || throw(ArgumentError(
                "semantically equivalent observations have conflicting numerical values"))
        end
        matched=isempty(same_interpretation) ? nothing : first(same_interpretation)
        if matched===nothing
            push!(keys,key);push!(products,product)
            push!(groups,(representative=index,members=[index],identities=[id]))
        else
            push!(groups[matched].members,index)
            push!(groups[matched].identities,id)
        end
    end
    return groups
end

function _description_values!(output,path,value;name=path,unit="",indices=(),text=nothing)
    if value isa NamedTuple
        descriptions=get(value,:field_descriptions,(;))
        for (key,child) in pairs(value)
            key in (:frequencies,:field_descriptions,:kind,:system_id,:cable_id) && continue
            child_path=isempty(path) ? string(key) : path*"."*string(key)
            field=get(descriptions,key,nothing)
            child_name=(field===nothing ? string(key) : field.name)*join("[$i]" for i in indices)
            _description_values!(output,child_path,child;name=child_name,
                unit=field===nothing ? unit : field.unit,indices,
                text=field===nothing ? nothing : get(field,:text,nothing))
        end
    elseif value isa Union{AbstractArray,Tuple}
        for (index,child) in enumerate(value)
            suffix="["*string(index)*"]"
            _description_values!(output,path*suffix,child;name=name*suffix,unit,indices=(indices...,index),text)
        end
    elseif value isa Union{Number,AbstractString,Symbol}
        output[path]=(;value,name,unit,text)
    end
    return output
end

"""
Describe differences in captured physical inputs, active methods, and individual
controls. Scientific text comes from completion-time owner descriptions; no
formula or problem is reconstructed. Administrative gridpoint IDs are retained
in the observations but do not form automatic legend prefixes.
"""
function observation_labels(observed;request=nothing)
    points=observed isa ObservedResult ? [observed] : observed
    isempty(points) && return String[]
    descriptions=[begin
        values=_description_values!(Dict{String,Any}(),"",get(point.gridpoint,:inputs,nothing))
        uncertainty=get(point.gridpoint,:uncertainty,nothing)
        annotations=get(point.gridpoint,:uncertainty_descriptions,nothing)
        if uncertainty isa NamedTuple && !isempty(uncertainty)
            annotations!==nothing && all(key -> haskey(annotations,key),keys(uncertainty)) ||
                throw(ArgumentError("retained uncertainty descriptions are absent; construct observations through the current UQ owner before plotting"))
            uncertainty=merge(uncertainty,(field_descriptions=annotations,))
        end
        _description_values!(values,"uncertainty",uncertainty)
        values
    end for point in points]
    paths=sort(unique(collect(Iterators.flatten(keys(record) for record in descriptions))))
    varying=filter(paths) do path
        haskey(first(descriptions),path) || return true
        original=first(descriptions)[path]
        any(record -> !haskey(record,path) ||
            !isequal(record[path].value,original.value) || record[path].unit!=original.unit,descriptions)
    end
    fields=map(points) do point
        retained=get(point.gridpoint,:formulation_fields,(;))
        family=request===nothing || isempty(retained) ? :all : Units.family(request_quantity(request))===Val(:series) ? :Z : :Y
        entries=get(retained,family,())
        all(field -> haskey(field,:meaning) && haskey(field,:control_fields),entries) ||
            throw(ArgumentError("retained formulation descriptions lack individual control meanings; capture descriptions with the current completion owner before plotting"))
        entries
    end
    # Only active peers establish a varying selection. A missing backend slot
    # is not an equal selection: the backend summary identifies that difference.
    peers(field)=[other for entries in fields for other in entries if other.meaning==field.meaning]
    varies(field)=any(other -> !isequal(other.selection.identifier,field.selection.identifier),peers(field))
    return map(eachindex(points)) do index
        parts=String[]
        methods=filter(varies,fields[index])
        for field in methods
            # When differing child methods already name the calculation, use
            # the backend itself instead of repeating its method summary.
            text=isempty(field.meaning) && all(other -> isempty(other.meaning),methods) ? field.summary : field.value
            if !isempty(field.name) && count(other -> !isempty(other.meaning),methods)>1
                text=field.name*"="*text
            end
            push!(parts,text)
        end
        for field in fields[index], control in field.control_fields
            active=[other for peer in peers(field) for other in peer.control_fields
                if other.scope==control.scope]
            any(other -> !isequal(other.value,control.value),active) || continue
            text=control.text
            isempty(field.name) || (text=field.name*": "*text)
            push!(parts,text)
        end
        for path in varying
            haskey(descriptions[index],path) || continue
            field=descriptions[index][path]
            push!(parts,field.text===nothing ?
                field.name*"="*string(field.value)*(isempty(field.unit) ? "" : " "*field.unit) : field.text)
        end
        if isempty(parts)
            roots=filter(field -> isempty(field.meaning),fields[index])
            if !isempty(roots)
                append!(parts,(field.summary for field in roots))
            else
                source_name=get(points[index].gridpoint,:name,nothing)
                push!(parts,source_name===nothing ? "Result $index" : string(source_name))
            end
        end
        join(parts,", ")
    end
end
