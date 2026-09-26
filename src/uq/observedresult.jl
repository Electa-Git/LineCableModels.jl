# Acquisition is owned by UQ because primary elements do not own empirical
# statistics, samples, histogram models, or sampling-workload information.
_is_uq_product(request) = request_identity(request) isa Tuple &&
    first(request_identity(request)) in (statistics,samples,histograms)

function Grammar.observation_requests(source::Union{MonteCarloResult,LinearErrorResult},
        requests::Tuple;point::Integer=firstindex(source),complete_pairs::Bool=false)
    primary=Tuple(filter(!_is_uq_product,requests))
    products=Tuple(filter(_is_uq_product,requests))
    normalized=isempty(requests) || !isempty(primary) ?
        Grammar.observation_requests(source[point],primary;complete_pairs) : (retained=(),displayed=())
    expanded=Tuple[]
    for request in products
        identity=request_identity(request)
        observation_request(source,request)
        if first(identity)===statistics && length(identity)==2
            for transform in (Statistics.mean,Statistics.std,minimum,
                    Base.Fix2(Statistics.quantile,0.05),Statistics.median,
                    Base.Fix2(Statistics.quantile,0.95),maximum)
                push!(expanded,(identity...,transform,request_indices(request)...))
            end
        else
            push!(expanded,request)
        end
    end
    selected=(normalized.retained...,expanded...)
    allunique(selected) || throw(ArgumentError("observation requests must be distinct"))
    return (retained=selected,displayed=(normalized.displayed...,expanded...))
end

function _uq_coordinates(core::Engine.LineParameters,selector,indices)
    rank=3
    selected=isempty(indices) ? ntuple(_ -> Colon(),rank) : indices
    length(selected)==rank || throw(DimensionMismatch("statistical products require row, column and sample indices"))
    return Engine.line_coordinates(core,(selector,selected...),nothing)
end
function _uq_coordinates(core::Engine.CableConstants,selector,indices)
    length(indices)<=1 || throw(DimensionMismatch("cable statistics select assembly indices"))
    selected=isempty(indices) ? Colon() : only(indices)
    return (kind=:assemblies,indices=(selected,),assemblies=observation_indices(selected,length(core)),
        labels=copy(core.cores),frequencies=[core.frequency],frequency_unit=Units.units(:base,:hertz),extent=(length(core),1))
end

function Grammar.observation_quantity(source::Union{MonteCarloResult,LinearErrorResult},
        point::Integer,request;unit=nothing,clip=true,atol=nothing,frequencies=nothing)
    core=source[point]
    _is_uq_product(request) || return Grammar.observation_quantity(core,request;unit,clip,atol,frequencies)
    identity=request_identity(request)
    product,selector=identity[1:2]
    indices=request_indices(request)
    q=request_quantity(request)
    target=unit===nothing ? Units.display_unit(q,basis(core)) : unit
    distribution=nothing
    if product===histograms
        rank=core isa Engine.LineParameters ? 3 : 1
        length(indices) in (rank,rank+1) || throw(DimensionMismatch("a histogram selects exactly one marginal and optional bin count"))
        marginal=indices[1:rank]
        all(i -> i isa Integer && !(i isa Bool),marginal) || throw(ArgumentError("histograms require integer marginal coordinates"))
        bins=length(indices)==rank ? nothing : last(indices)
        model=observe(source,histograms,selector,point,marginal...,bins)
        T=typeof(float(first(model.edges)))
        factor=Units.scale_factor(Units.native_unit(q,basis(core)),target,T)
        model=detach(model,factor)
        retained_samples=samples(source)===nothing ? nothing :
            vec(detach(observe(source,samples,selector,point,marginal...,Colon()),factor))
        edges=model.edges
        probabilities=model.density.*diff(edges)
        counts=retained_samples===nothing ? nothing : begin
            count=zeros(Int,length(model.density))
            for value in retained_samples
                index=value==last(edges) ? length(count) : clamp(searchsortedlast(edges,value),1,length(count))
                count[index]+=1
            end
            count
        end
        empirical=retained_samples===nothing ? nothing :
            (x=sort(retained_samples),y=collect(1:length(retained_samples))./length(retained_samples))
        model_cdf=(x=copy(edges),y=cumulative_probability.(Ref(model),edges))
        qq=retained_samples===nothing ? nothing : quantile_pairs(model,retained_samples)
        distribution=(edges=copy(edges),empirical_cdf=empirical,model_cdf,qq)
        values=(lower=edges[1:end-1],upper=edges[2:end],density=copy(model.density),
            probability=probabilities,count=counts===nothing ? fill(missing,length(probabilities)) : counts)
        coordinates=merge(_uq_coordinates(core,selector,marginal),(kind=:histogram,))
        return (request,quantity=q,family=:statistics,statistic=:histogram,values,unit=target,basis=basis(core),coordinates,
            thresholds=nothing,available=true,engineering_zero=false,clipped=false,missing_reason=nothing,
            assumptions=nothing,distribution,ordinate_units=(density=Units.UnitExpr(target.denominator,target.numerator),
                probability=Units.units(:base,:dimensionless),count=Units.units(:base,:dimensionless)))
    end
    acquisition_indices=indices
    if product===samples && !isempty(indices)
        rank=core isa Engine.LineParameters ? 3 : 1
        length(indices)==rank && (acquisition_indices=(indices...,Colon()))
    end
    values=observe(source,identity...,point,acquisition_indices...)
    T=typeof(float(real(nominal(zero(eltype(values))))))
    factor=Units.scale_factor(Units.native_unit(q,basis(core)),target,T)
    coordinate_indices=product===samples && !isempty(indices) ? indices[1:(core isa Engine.LineParameters ? 3 : 1)] : indices
    coordinates=_uq_coordinates(core,selector,coordinate_indices)
    if product===samples
        rank=core isa Engine.LineParameters ? 3 : 1
        trial_index=length(indices)>rank ? indices[rank+1] : Colon()
        coordinates=merge(coordinates,(kind=:samples,indices=(coordinates.indices...,trial_index),trials=observation_indices(trial_index,trial_count(source,point))))
    end
    # These are selected estimators or individual trials, not replacements for
    # primary uncertain numbers. In particular, a retained std is never clipped.
    available=Engine.resolution_available.(values)
    retained=map((value,valid) -> valid ? value : missing,values,available)
    statistic=product===samples ? :samples : last(identity) isa Base.Fix2 ?
        Symbol("quantile_",last(identity).x) : nameof(last(identity))
    return (request,quantity=q,family=:statistics,statistic,values=detach(retained,factor),
        unit=target,basis=basis(core),coordinates,thresholds=nothing,available,
        engineering_zero=false,clipped=false,missing_reason=map(valid -> valid ? nothing : :nonfinite_value,available),
        assumptions=nothing,distribution)
end

"""
$(TYPEDSIGNATURES)

Detach one UQ point using its owner's point-indexed accessors. Product requests
omit the point index because `point` already identifies the atomic observation.
For example, `(statistics,R,mean)` retains an estimator and
`(histograms,R,1,2,3,20)` retains a marginal distribution with up to twenty bins.
Primary requests follow the ordinary primary-result pair rules.
"""
function Grammar.ObservedResult(source::Union{MonteCarloResult,LinearErrorResult},point::Integer,
        requests::Tuple=();comparisons=(),timings=(;),gridpoint=nothing,clip::Bool=true,
        atol=nothing,units::Tuple=(),length_unit::Symbol=:kilo,frequency_unit::Symbol=:base,
        quantity_units=nothing,frequencies=nothing,complete_pairs::Bool=false)
    selected=Grammar.observation_requests(source,requests;point,complete_pairs).retained
    isempty(units) || length(units)==length(selected) || throw(DimensionMismatch("units must align with retained requests"))
    targets = if isempty(units)
        Grammar.unit_targets(selected, basis(source);
            length_prefix=length_unit, overrides=quantity_units)
    else
        map(selected, units) do request, unit
            Units.display_unit(request_quantity(request), basis(source), unit;
                length_prefix=length_unit)
        end
    end
    description=gridpoint===nothing ? Grammar.observation_gridpoint(source[point]) : gridpoint
    sampling=source isa MonteCarloResult ? merge(confidence(source,point),
        (frequencies=source[point] isa Engine.LineParameters ? detach(Engine.frequencies(source[point])) :
            [source[point].frequency],basis=basis(source))) : nothing
    description=merge(description,(sampling,
        uncertainty_descriptions=_uncertainty_descriptions(source isa MonteCarloResult ? MonteCarlo : LinearError),
        uncertainty=(estimator=source isa MonteCarloResult ? :empirical : :first_order,
        representation=source isa MonteCarloResult ? :marginal_mean_std : :dependency_preserving)))
    quantities=map(eachindex(selected)) do index
        request=selected[index]
        unit=targets[index]
        record=Grammar.observation_quantity(source,point,request;unit,clip,atol,frequencies)
        if record.coordinates.frequencies!==nothing
            f=record.coordinates.frequencies
            target=Units.units(frequency_unit,:hertz)
            factor=Units.scale_factor(Units.units(:base,:hertz),target,typeof(float(nominal(first(f)))))
            record=merge(record,(coordinates=merge(record.coordinates,(frequencies=detach(f,factor),frequency_unit=target)),))
        end
        record
    end
    return Grammar.ObservedResult(description,quantities,collect(comparisons),timings)
end

function observables(source::Union{MonteCarloResult,LinearErrorResult},requests::Tuple=();
        comparisons=(),timings=(;),kwargs...)
    return [begin
        id=Grammar.observation_gridpoint(source[point]).id
        errors=filter(record -> record.candidate_id==id,comparisons)
        recorded=timings isa NamedTuple ? timings : begin
            matched=filter(record -> record.candidate_id==id,timings)
            length(matched)<=1 || throw(ArgumentError("duplicate point timings"))
            isempty(matched) ? (;) : only(matched)
        end
        Grammar.ObservedResult(source,point,requests;comparisons=errors,timings=recorded,kwargs...)
    end for point in eachindex(source)]
end
