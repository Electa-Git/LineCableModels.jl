# Re-expression uses recorded units. It never enters source acquisition or
# changes the retained numerical eligibility decision.
_observed_scale(value::Missing,from,to) = value
_observed_scale(value::Nothing,from,to) = value
_observed_scale(value::AbstractArray,from,to) = map(x -> _observed_scale(x,from,to),value)
_observed_scale(value::NamedTuple,from,to) = map(x -> _observed_scale(x,from,to),value)
_observed_scale(value::Tuple,from,to) = map(x -> _observed_scale(x,from,to),value)
function _observed_scale(value::Number,from,to)
    from==to && return value
    T=typeof(float(real(nominal(value))))
    convert_value()=value*scale_factor(from,to,T)
    return T===BigFloat ? setprecision(convert_value,BigFloat,precision(real(nominal(value)))) : convert_value()
end

function _reexpress_product(product;unit=nothing,frequency_unit=nothing)
    target=unit===nothing ? product.unit : unit
    scale_factor(product.unit,target) # Check dimensions even if all values are missing.
    result=product
    if target!=product.unit
        values=if product.coordinates.kind===:histogram
            merge(product.values,(
                lower=_observed_scale(product.values.lower,product.unit,target),
                upper=_observed_scale(product.values.upper,product.unit,target),
                density=_observed_scale(product.values.density,inv(product.unit),inv(target))))
        else
            _observed_scale(product.values,product.unit,target)
        end
        thresholds=product.thresholds
        if thresholds!==nothing && thresholds.unit==product.unit
            thresholds=merge(thresholds,(values=_observed_scale(thresholds.values,product.unit,target),unit=target))
        end
        result=merge(product,(values,unit=target,thresholds))
        if product.coordinates.kind===:histogram
            d=product.distribution
            cdf(c)=c===nothing ? nothing : merge(c,(x=_observed_scale(c.x,product.unit,target),))
            distribution=merge(d,(edges=_observed_scale(d.edges,product.unit,target),
                empirical_cdf=cdf(d.empirical_cdf),model_cdf=cdf(d.model_cdf),
                qq=d.qq===nothing ? nothing : _observed_scale(d.qq,product.unit,target)))
            result=merge(result,(distribution,ordinate_units=merge(product.ordinate_units,(density=inv(target),))))
        end
    end
    c=result.coordinates
    if frequency_unit!==nothing && haskey(c,:frequencies) && c.frequencies!==nothing
        target_frequency=frequency_unit isa Symbol ? Units.units(frequency_unit,:hertz) : frequency_unit
        scale_factor(c.frequency_unit,target_frequency)
        result=merge(result,(coordinates=merge(c,(
            frequencies=_observed_scale(c.frequencies,c.frequency_unit,target_frequency),
            frequency_unit=target_frequency)),))
    end
    return result
end

function _retained_unit(product,override,length_unit)
    override isa UnitExpr && return override
    override===nothing || override isa Symbol || throw(ArgumentError("unit override must be a unit expression or metric prefix"))
    numerator=product.unit.numerator
    if override!==nothing && !isempty(numerator)
        numerator=(Units.Unit(first(numerator).name,override),Base.tail(numerator)...)
    end
    denominator=map(product.unit.denominator) do unit
        length_unit!==nothing && unit.name===:meter ? Units.Unit(:meter,length_unit) : unit
    end
    return UnitExpr(numerator,denominator)
end

"""
$(TYPEDSIGNATURES)

Select or re-express an existing observation without acquiring a source.
Omitted unit options preserve recorded units. Explicit units convert from each
product's recorded unit; frequency conversion preserves sample identities.
Coordinates, availability, uncertainty dependencies, comparisons, and timing
associations remain retained. New clipping, thresholds, or frequency samples
require a new observation of the primary result.
"""
function ObservedResult(source::ObservedResult,requests::Tuple=();
        comparisons=nothing,timings=nothing,gridpoint=nothing,clip=nothing,atol=nothing,
        units::Tuple=(),length_unit::Union{Nothing,Symbol}=nothing,
        frequency_unit=nothing,quantity_units=nothing,frequencies=nothing,
        complete_pairs::Bool=false)
    atol===nothing && frequencies===nothing || throw(ArgumentError(
        "retained observations cannot apply new cutoffs or recover frequency samples"))
    gridpoint===nothing || isequal(gridpoint,source.gridpoint) || throw(ArgumentError("retained gridpoint identity cannot be replaced"))
    comparisons===nothing || isequal(collect(comparisons),source.errors) || throw(ArgumentError("retained comparisons cannot be replaced"))
    timings===nothing || isequal(timings,source.timings) || throw(ArgumentError("retained timing associations cannot be replaced"))
    selected=observation_requests(source,requests;complete_pairs).retained
    isempty(units) || length(units)==length(selected) || throw(DimensionMismatch("units must align with retained requests"))
    isempty(units) || quantity_units===nothing || throw(ArgumentError("use units or quantity_units, not both"))
    products=map(eachindex(selected)) do index
        request=selected[index]
        product=observation_product(source,request)
        clip===nothing || clip===product.clipped || throw(ArgumentError("retained clipping cannot be changed"))
        override=isempty(units) ? _unit_override(quantity_units,request) : units[index]
        unit=_retained_unit(product,override,length_unit)
        _reexpress_product(product;unit,frequency_unit)
    end
    return ObservedResult(source.gridpoint,products,source.errors,source.timings)
end

"""
$(TYPEDSIGNATURES)

Interpret a retained request across observations. Convert compatible
quantity and frequency units to the first product's units (or explicit targets),
and order coefficients by the first product's original coordinates. Every trace
keeps its own samples. Missing coefficients and incompatible units fail without
interpolation, numerical acquisition, or changes to scientific eligibility.
With `band`, select original sample identities from completed comparison
records. `reference_id` disambiguates the recorded reference; a separately
included reference uses the same unambiguous saved selection. Missing or
conflicting records fail before returning any products.
"""
function observation_product(points::Union{Tuple,AbstractVector{<:ObservedResult}},request;
        unit=nothing,frequency_unit=nothing,band=nothing,reference_id=nothing)
    isempty(points) && throw(ArgumentError("at least one observation is required"))
    first_product=observation_product(first(points),request)
    target=something(unit,first_product.unit)
    coordinate=first_product.coordinates
    frequency_target=frequency_unit===nothing ? get(coordinate,:frequency_unit,nothing) : frequency_unit
    selections=band===nothing ? nothing : _retained_band_samples(points,request,band,reference_id)
    products=map(eachindex(points)) do index
        point=points[index]
        product=observation_product(point,request)
        c=product.coordinates
        c.kind==coordinate.kind || throw(ArgumentError("overlaid products require the same coordinate kind"))
        if c.kind in (:matrix,:diagonal,:vector)
            aligned=c.kind===:vector ? Set(c.positions)==Set(coordinate.positions) :
                Set(c.rows)==Set(coordinate.rows) && Set(c.columns)==Set(coordinate.columns)
            aligned ||
                throw(DimensionMismatch("overlaid products must retain the requested original coefficients"))
            reordered=c.kind===:vector ? c.positions!=coordinate.positions :
                c.rows!=coordinate.rows || c.columns!=coordinate.columns
            if reordered
                identity=request_identity(product.request)
                prefix=identity isa Tuple ? identity : (identity,)
                indices=c.kind===:matrix ? (coordinate.rows,coordinate.columns,Colon()) :
                    c.kind===:vector ? (coordinate.positions,Colon()) : (coordinate.rows,Colon())
                product=_selected_product(product,(prefix...,indices...),indices)
            end
        end
        if selections!==nothing
            c=product.coordinates
            c.kind in (:matrix,:diagonal,:vector) || throw(ArgumentError(
                "retained band selection requires frequency coordinates"))
            samples=filter(in(selections[index]),c.samples)
            identity=request_identity(product.request)
            prefix=identity isa Tuple ? identity : (identity,)
            indices=c.kind===:matrix ? (c.rows,c.columns,samples) :
                c.kind===:vector ? (c.positions,samples) : (c.rows,samples)
            product=_selected_product(product,(prefix...,indices...),indices)
        end
        detach(_reexpress_product(product;unit=target,frequency_unit=frequency_target))
    end
    selections===nothing || any(p -> !isempty(p.coordinates.samples),products) ||
        throw(ArgumentError("comparison band $(repr(band)) contains no retained samples in the supplied observations"))
    return products
end

# Comparison records, not plotting-side frequency rules, define a retained band.
function _retained_band_samples(points,request,band,reference_id)
    matching=map(points) do point
        id=get(point.gridpoint,:id,nothing)
        filter(point.errors) do row
            isequal(row.candidate_id,id) && isequal(row.band,band) &&
                request_identity(row.request)==request_identity(request) &&
                (reference_id===nothing || isequal(row.reference_id,reference_id))
        end
    end
    records=collect(Iterators.flatten(matching))
    isempty(records) && throw(ArgumentError("no completed comparison retains band $(repr(band)) for this request and reference"))
    references=unique(row.reference_id for row in records)
    length(references)==1 || throw(ArgumentError("retained band has multiple references; supply reference_id"))
    retained_reference=only(references)
    definition(row)=(;
        indices=row.settings.indices,
        requested_bounds=get(row.settings,:requested_bounds,nothing),
        actual_bounds=get(row.settings,:actual_bounds,nothing),
        fundamental=get(row.settings,:fundamental,nothing),
        harmonics=get(row.settings,:harmonics,nothing))
    first_definition=definition(first(records))
    all(row -> isequal(definition(row),first_definition),records) || throw(ArgumentError(
        "completed comparisons disagree on the saved definition or sample selection of band $(repr(band))"))
    return map(eachindex(points)) do index
        id=get(points[index].gridpoint,:id,nothing)
        !isempty(matching[index]) || isequal(id,retained_reference) || throw(ArgumentError(
            "observation $(repr(id)) has no completed comparison for band $(repr(band)) and the selected reference"))
        copy(first_definition.indices)
    end
end

# Direct owner dispatch obeys the same retained-input requirements as construction.
function observation_quantity(source::ObservedResult,request;unit=nothing,clip=nothing,atol=nothing,frequencies=nothing)
    atol===nothing && frequencies===nothing || throw(ArgumentError("retained quantities cannot apply new cutoffs or samples"))
    product=observation_product(source,request;unit)
    clip===nothing || clip===product.clipped || throw(ArgumentError("retained clipping cannot be changed"))
    return detach(product)
end
