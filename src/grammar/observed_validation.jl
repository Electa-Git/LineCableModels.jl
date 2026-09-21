# These checks describe the ordinary records consumed by tables, plots, and
# archives. There is no registry or second product representation.
function _observed_fields(record,fields,what)
    record isa NamedTuple && all(key -> haskey(record,key),fields) ||
        throw(ArgumentError("$what requires fields $(fields)"))
end
_observed_length(value::AbstractArray)=length(value)
_observed_length(value)=1

function _observed_indices(indices,extent,name)
    indices isa AbstractVector && allunique(indices) && all(i -> i isa Integer && !(i isa Bool) && 1<=i<=extent,indices) ||
        throw(ArgumentError("$name must contain distinct original indices within its extent"))
end

function _validate_observed_quantity(product)
    _observed_fields(product,(:request,:quantity,:family,:statistic,:values,:unit,:basis,
        :coordinates,:thresholds,:available,:engineering_zero,:clipped,:missing_reason),"observed quantity")
    request_quantity(product.request)==product.quantity || throw(ArgumentError("request and quantity disagree"))
    product.family isa Symbol && product.statistic isa Symbol && product.clipped isa Bool ||
        throw(ArgumentError("family/statistic must be symbols and clipped must be Boolean"))
    scale_factor(native_unit(product.quantity,product.basis),product.unit)
    c=product.coordinates
    _observed_fields(c,(:kind,:indices,:extent),"quantity coordinates")
    c.extent isa Tuple && all(n -> n isa Integer && n>=0,c.extent) || throw(ArgumentError("invalid coordinate extent"))
    c.indices isa Tuple || throw(ArgumentError("coordinate indices must be a tuple"))
    declared=request_indices(product.request)
    c.kind===:samples && length(declared)==length(c.indices)-1 &&
        (declared=(declared...,Colon()))
    if !isempty(declared) && c.kind!==:array
        length(declared)>=length(c.indices) && declared[1:length(c.indices)]==c.indices ||
            throw(ArgumentError("request selectors disagree with retained coordinate selectors"))
    end
    dims=if c.kind in (:matrix,:diagonal) || c.kind in (:samples,:histogram) && haskey(c,:rows)
        _observed_fields(c,(:rows,:columns,:samples,:frequencies,:frequency_unit,:labels,:domain),"matrix coordinates")
        length(c.extent)==3 || throw(DimensionMismatch("matrix extent requires three axes"))
        _observed_indices(c.rows,c.extent[1],"rows")
        _observed_indices(c.columns,c.extent[2],"columns")
        _observed_indices(c.samples,c.extent[3],"samples")
        length(c.labels)>=max(c.extent[1],c.extent[2]) || throw(DimensionMismatch("matrix labels do not cover the extent"))
        c.kind===:diagonal && c.rows!=c.columns && throw(ArgumentError("diagonal coordinates must agree"))
        c.frequencies===nothing || length(c.frequencies)==length(c.samples) || throw(DimensionMismatch("frequency and sample counts differ"))
        c.kind===:diagonal ? (length(c.rows),length(c.samples)) : (length(c.rows),length(c.columns),length(c.samples))
    elseif c.kind===:assemblies || c.kind in (:samples,:histogram) && haskey(c,:assemblies)
        _observed_fields(c,(:assemblies,:labels,:frequencies,:frequency_unit),"assembly coordinates")
        length(c.extent)==2 && c.extent[2]==1 || throw(DimensionMismatch("assembly extent requires one operating frequency"))
        _observed_indices(c.assemblies,c.extent[1],"assemblies")
        length(c.labels)==c.extent[1] && allunique(c.labels) || throw(ArgumentError("assembly labels must identify every assembly uniquely"))
        length(c.frequencies)==1 || throw(DimensionMismatch("assemblies require one operating frequency"))
        (length(c.assemblies),)
    elseif c.kind===:array
        c.extent
    else
        throw(ArgumentError("unsupported observed coordinate kind $(c.kind)"))
    end
    if haskey(c,:frequency_unit)
        scale_factor(c.frequency_unit,Units.units(:base,:hertz))
        c.frequencies===nothing || all(f -> f isa Real && isfinite(nominal(f)) && nominal(f)>=0,c.frequencies) ||
            throw(ArgumentError("frequency coordinates must be finite and nonnegative"))
    end
    if c.kind===:samples
        _observed_fields(c,(:trials,),"sample coordinates")
        allunique(c.trials) && all(t -> t isa Integer && !(t isa Bool) && t>0,c.trials) ||
            throw(ArgumentError("trial coordinates must be distinct positive integers"))
        dims=(dims...,length(c.trials))
    end
    if c.kind===:histogram
        _observed_fields(product.values,(:lower,:upper,:density,:probability,:count),"histogram values")
        n=length(product.values.lower)
        all(v -> v isa AbstractVector && length(v)==n,Base.values(product.values)) || throw(DimensionMismatch("histogram column lengths differ"))
        _observed_fields(product,(:distribution,:ordinate_units),"histogram product")
        _observed_fields(product.distribution,(:edges,:empirical_cdf,:model_cdf,:qq),"histogram distribution")
        length(product.distribution.edges)==n+1 || throw(DimensionMismatch("histogram edges and bins differ"))
        for cdf in (product.distribution.empirical_cdf,product.distribution.model_cdf)
            cdf===nothing && continue
            _observed_fields(cdf,(:x,:y),"CDF coordinates")
            cdf.x isa AbstractVector && cdf.y isa AbstractVector && length(cdf.x)==length(cdf.y) ||
                throw(DimensionMismatch("CDF coordinate lengths differ"))
        end
        qq=product.distribution.qq
        if qq!==nothing
            _observed_fields(qq,(:model,:sample,:reference),"Q–Q coordinates")
            qq.model isa AbstractVector && qq.sample isa AbstractVector && length(qq.model)==length(qq.sample) ||
                throw(DimensionMismatch("Q–Q model and sample coordinate lengths differ"))
            qq.reference isa Tuple && length(qq.reference)==2 || throw(ArgumentError("Q–Q identity line requires two endpoints"))
        end
        _observed_fields(product.ordinate_units,(:density,:probability,:count),"histogram ordinate units")
        scale_factor(inv(product.unit),product.ordinate_units.density)
        for unit in (product.ordinate_units.probability,product.ordinate_units.count)
            scale_factor(Units.units(:base,:dimensionless),unit)
        end
        product.available isa Bool && product.engineering_zero isa Bool && product.missing_reason===nothing ||
            throw(ArgumentError("histogram availability must be scalar"))
        product.thresholds===nothing || throw(ArgumentError("histograms cannot carry primary clipping thresholds"))
        return nothing
    end
    values=product.values
    values isa Union{Number,Missing,AbstractArray} || throw(ArgumentError("numerical products require scalar or array values"))
    _observed_length(values)==prod(dims) || throw(DimensionMismatch("value count differs from retained coordinates"))
    actual=values isa AbstractArray ? size(values) : ()
    reduced=length(c.indices)==length(dims) ? Tuple(n for (n,i) in zip(dims,c.indices) if !(i isa Integer)) : dims
    actual in (dims,reduced) || throw(DimensionMismatch("value shape differs from retained coordinate axes"))
    for (field,allowed) in ((:available,x -> x isa Bool),(:engineering_zero,x -> x isa Bool),
            (:missing_reason,x -> x===nothing || x isa Symbol))
        value=getproperty(product,field)
        value===nothing && continue
        if value isa AbstractArray
            size(value)==actual || (length(value)==1 && isempty(actual)) || throw(DimensionMismatch("$field shape differs from values"))
            all(allowed,value) || throw(ArgumentError("invalid $field entries"))
        else
            allowed(value) || throw(ArgumentError("invalid $field"))
        end
    end
    if product.thresholds!==nothing
        _observed_fields(product.thresholds,(:kind,:values,:unit),"quantity thresholds")
        # Phase eligibility retains physical Cartesian thresholds, not angles.
        threshold_unit=product.thresholds.unit
        phase=request_identity(product.request) isa Tuple && angle in request_identity(product.request)
        threshold_quantity=phase ? quantity(first(request_identity(product.request))) : product.quantity
        scale_factor(native_unit(threshold_quantity,product.basis),threshold_unit)
        samples=haskey(c,:samples) ? length(c.samples) : nothing
        _validate_observed_cutoff(product.thresholds.values,samples)
    end
    return nothing
end

function _validate_observed_cutoff(value,samples)
    value===nothing && return nothing
    if value isa NamedTuple
        foreach(v -> _validate_observed_cutoff(v,samples),values(value))
    elseif value isa AbstractVector
        samples===nothing || length(value)==samples || throw(DimensionMismatch("cutoff count differs from samples"))
        foreach(v -> _validate_observed_cutoff(v,nothing),value)
    else
        value isa Real && isfinite(nominal(value)) && nominal(value)>=0 || throw(ArgumentError("cutoffs must be finite and nonnegative"))
    end
    return nothing
end

function _validate_observed_comparison(row)
    _observed_fields(row,(:candidate_id,:reference_id,:request,:quantity,:statistic,:band,
        :normalization,:absolute,:relative,:absolute_unit,:relative_unit,:coordinates,:settings,:maxima),"completed comparison")
    for id in (row.candidate_id,row.reference_id)
        _validate_observed_id(id)
    end
    _observed_fields(row.settings,(:basis,:indices,:sample_count,:status,:normalization_reason),"comparison settings")
    _observed_fields(row.maxima,(:absolute,:relative),"comparison maxima")
    request_quantity(row.request)==row.quantity || throw(ArgumentError("comparison request and quantity disagree"))
    row.absolute isa AbstractMatrix && row.relative isa AbstractMatrix && size(row.absolute)==size(row.relative) ||
        throw(DimensionMismatch("comparison matrices must have matching extents"))
    length(row.coordinates)==size(row.absolute,1)==size(row.absolute,2) ||
        throw(DimensionMismatch("comparison coordinates and matrices disagree"))
    size(row.settings.status)==size(row.absolute)==size(row.settings.normalization_reason) ||
        throw(DimensionMismatch("comparison status and reason matrices must match the error matrices"))
    row.settings.sample_count==length(row.settings.indices) || throw(DimensionMismatch("comparison sample count disagrees with indices"))
    allunique(row.settings.indices) && all(i -> i isa Integer && !(i isa Bool) && i>0,row.settings.indices) ||
        throw(ArgumentError("comparison sample indices must be distinct positive integers"))
    for maximum in (row.maxima.absolute,row.maxima.relative)
        _observed_fields(maximum,(:value,:index),"comparison maximum")
    end
    all(ismissing.(row.absolute).==ismissing.(row.relative)) || throw(ArgumentError("absolute and relative errors must share eligibility"))
    scale_factor(native_unit(row.quantity,row.settings.basis),row.absolute_unit)
    scale_factor(Units.units(:base,:dimensionless),row.relative_unit)
    return nothing
end

function _validate_observed_id(id)
    _observed_fields(id,(:source_id,:problem_index,:formulation_index),"gridpoint identity")
    id.source_id===nothing && throw(ArgumentError("identified gridpoints require a source identity"))
    all(i -> i isa Integer && !(i isa Bool) && i>0,(id.problem_index,id.formulation_index)) ||
        throw(ArgumentError("point and formulation indices must be positive integers"))
    return nothing
end
