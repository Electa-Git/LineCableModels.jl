"""
$(TYPEDEF)

Publish an [`Engine.CableConstants`](@ref) result as separate R/L/C/G tables,
each with one operating-frequency row and a named column per assembly.

$(TYPEDFIELDS)
"""
struct CableConstantsTableDefinition <: AbstractReportDefinition
    "Whether detached display residue is replaced with exact zero."
    clip::Bool
end
CableConstantsTableDefinition() = CableConstantsTableDefinition(true)

"""
$(TYPEDEF)

Define requests and display units for separate line-parameter quantity tables.
Each full matrix table has one frequency column and every ordered coefficient.

$(TYPEDFIELDS)
"""
struct LineParametersTableDefinition{Q <: Tuple, U} <: AbstractReportDefinition
    "Explicit observable requests in quantity-table order."
    requests::Q
    "SI prefix used to display frequency."
    frequency_unit::Symbol
    "Length prefix used for per-length quantities."
    length_unit::Symbol
    "Optional display-unit overrides resolved from the requests."
    quantity_units::U
    "Whether detached display residue is replaced with exact zero."
    clip::Bool
end

"""
$(TYPEDEF)

Select and publish per-term comparisons of scalar or formulation-space results.
Retained products supply detailed tables and a compact summary.

$(TYPEDFIELDS)
"""
struct BenchmarkTableDefinition{S <: NamedTuple, I, O <: NamedTuple} <: AbstractReportDefinition
    "Normalized quantities, frequency bands and numerical comparison controls."
    settings::S
    "Explicitly supplied comparison controls, used when selecting retained products."
    requested::Tuple{Vararg{Symbol}}
    "Whether detached display residue is replaced with exact zero."
    clip::Bool
    "Optional explicit illustration request; nothing produces no figure."
    illustration::I
    "Options for the explicitly requested illustration."
    plot_options::O
end

"""
$(TYPEDSIGNATURES)

Request per-term comparisons grouped by formulation and frequency band.
Quantities use the observation grammar. The default bands are the entire range,
near DC, harmonic, narrowband and wideband. `fundamental` is in Hz. No figure is
created unless `illustration` is explicitly supplied. A benchmark report retains one scalar reference separately from its candidates.
The explicit comparison operation also supports declared collection pairings.
"""
function BenchmarkTableDefinition(; clip::Bool=false, illustration=nothing, plot_options=(;), kwargs...)
    haskey(kwargs,:requests) && any(key -> haskey(kwargs,key),(:quantities,:statistics)) &&
        throw(ArgumentError("use requests or quantities/statistics, not both"))
    statistics = Tuple(get(kwargs,:statistics,(:value,)))
    quantities = Tuple(get(kwargs,:quantities,statistics == (:value,) ? (Z,Y,R,L,G,C) : (R,L,C,G)))
    requests = get(kwargs,:requests,nothing)
    if requests === nothing
        selectors = map(quantities) do value
            value isa Function && return value
            matches = filter(selector -> nameof(selector) == value, (Z,Y,R,X,L,G,B,C))
            length(matches) == 1 || throw(ArgumentError("unknown physical quantity $value"))
            only(matches)
        end
        transforms = map(statistics) do entry
            entry === :value && return nothing
            entry isa Function && return entry
            matches = filter(selector -> nameof(selector) == entry,
                (Statistics.mean,Statistics.std,Statistics.median,minimum,maximum))
            length(matches) == 1 || throw(ArgumentError("unknown statistic $entry"))
            only(matches)
        end
        requests = Tuple(transform === nothing ? selector : (UQ.statistics,selector,transform)
            for selector in selectors for transform in transforms)
    end
    controls = (; (key=>value for (key,value) in kwargs if !(key in (:requests,:quantities,:statistics)))...)
    result = BenchmarkTableDefinition(Tuple(requests);clip,illustration,plot_options,controls...)
    requested = Tuple(unique(key in (:quantities,:statistics) ? :requests : key for key in keys(kwargs)))
    return BenchmarkTableDefinition(result.settings,requested,clip,illustration,plot_options)
end

"""
$(TYPEDSIGNATURES)

Select scientific requests for per-term RMS comparisons. Statistical requests
use `(statistics, quantity, statistic)` function tuples. Bands apply equally to
deterministic and statistical products. Absolute limits use native physical
units; relative errors are dimensionless. Plotting remains explicitly requested.
"""
function BenchmarkTableDefinition(requests::Tuple; clip::Bool=false, illustration=nothing,
        plot_options=(;), kwargs...)
    defaults=(bands=(:all,:dc,:harmonic,:narrow,:wide),
        normalizations=(:reference_rms,), atol=nothing, fundamental=50.0, harmonics=50,
        unsupported=(;), pairing=nothing)
    isempty(setdiff(keys(kwargs),keys(defaults))) || throw(ArgumentError("unknown benchmark comparison controls"))
    supplied=merge(defaults,(;kwargs...))
    settings=merge((;requests),supplied,(
        bands=Tuple(supplied.bands),
        normalizations=Tuple(supplied.normalizations)))
    return validate(BenchmarkTableDefinition(settings,(:requests,keys(kwargs)...),clip,illustration,plot_options))
end

function BenchmarkTableDefinition(clip::Bool; kwargs...)
    return BenchmarkTableDefinition(; clip, kwargs...)
end

"""Validate comparison requests before numerical execution."""
function validate(definition::BenchmarkTableDefinition)
    settings=definition.settings
    !isempty(settings.requests) && allunique(settings.requests) ||
        throw(ArgumentError("benchmark requests must be nonempty and distinct"))
    for request in settings.requests
        request_quantity(request)
        isempty(request_indices(request)) || throw(ArgumentError("benchmark point selection uses pairing"))
    end
    if settings.pairing !== nothing
        settings.pairing isa Union{Tuple,AbstractVector} && !isempty(settings.pairing) &&
            all(pair -> pair isa Tuple{Integer,Integer} && all(index -> !(index isa Bool) && index>0,pair),settings.pairing) &&
            sort(last.(collect(settings.pairing))) == collect(1:length(settings.pairing)) ||
            throw(ArgumentError("pairing must list positive reference/candidate indices with each candidate exactly once"))
    end
    !isempty(settings.bands) && allunique(settings.bands) ||
        throw(ArgumentError("benchmark needs distinct frequency bands"))
    !isempty(settings.normalizations) && allunique(settings.normalizations) ||
        throw(ArgumentError("benchmark needs distinct RMS normalizations"))
    for band in settings.bands, normalization in settings.normalizations
        validate(Engine.compare; band, normalization, atol=settings.atol,
            fundamental=settings.fundamental, harmonics=settings.harmonics,
            unsupported=settings.unsupported)
    end
    return definition
end


LineParametersTableDefinition(requests::Tuple=();frequency_unit::Symbol=:base,
    length_unit::Symbol=:kilo,quantity_units=nothing,clip::Bool=true) =
    LineParametersTableDefinition(requests,frequency_unit,length_unit,quantity_units,clip)

_bound_name(selector::Base.Fix2,statistic) =
    selector.f in (LineCableModels.ModalAnalysis.H,LineCableModels.ModalAnalysis.Zc,
        LineCableModels.ModalAnalysis.Yc) ?
        string(nameof(selector.f),"_phase",
            selector.f===LineCableModels.ModalAnalysis.H ? "_$(selector.x.field)" : "") :
        string(statistic)

_quantity_name(product) = begin
    identity=request_identity(product.request)
    identity isa Function ? (identity isa Base.Fix2 ? Symbol(_bound_name(identity,product.statistic)) : nameof(identity)) :
        Symbol(join([entry isa Base.Fix2 ? _bound_name(entry,product.statistic) : string(nameof(entry)) for entry in identity],"_"))
end

function _quantity_table(product; gridpoint_id=nothing)
    coordinates=product.coordinates
    values=product.values
    scalar=values isa Number || ismissing(values)
    if coordinates.kind in (:matrix,:diagonal,:vector)
        diagonal=coordinates.kind===:diagonal
        vector=coordinates.kind===:vector
        dimensions=vector ? (length(coordinates.positions),length(coordinates.samples)) :
            diagonal ? (length(coordinates.rows),length(coordinates.samples)) :
            (length(coordinates.rows),length(coordinates.columns),length(coordinates.samples))
        shaped=reshape(scalar ? [values] : values,dimensions...)
        f=coordinates.frequencies
        table=f===nothing ? DataFrame(sample=coordinates.samples) : DataFrame(frequency=f)
        columns=Pair{Symbol,Any}[]
        if vector
            for (i,position) in enumerate(coordinates.positions)
                name=Symbol(coordinates.axis_label," ",coordinates.labels[position])
                table[!,name]=copy(shaped[i,:])
                push!(columns,name=>(quantity=product.quantity,unit=product.unit,
                    axis=coordinates.axis,position,label=coordinates.labels[position]))
            end
        elseif diagonal
            for (i,row) in enumerate(coordinates.rows)
                name=Symbol("[",row,",",row,"]")
                table[!,name]=copy(shaped[i,:])
                push!(columns,name=>(quantity=product.quantity,unit=product.unit,row,column=row))
            end
        else
            # Full matrices are deliberately row-major, including both off-diagonals.
            for (i,row) in enumerate(coordinates.rows), (j,column) in enumerate(coordinates.columns)
                name=get(coordinates,:column_domain,nothing)===:ModalDomain ?
                    Symbol("Conductor ",coordinates.labels[row],", Mode ",coordinates.column_labels[column]) :
                    Symbol("[",row,",",column,"]")
                table[!,name]=copy(shaped[i,j,:])
                push!(columns,name=>(quantity=product.quantity,unit=product.unit,row,column,
                    row_label=coordinates.labels[row],
                    column_label=get(coordinates,:column_labels,coordinates.labels)[column]))
            end
        end
        first_column=f===nothing ? (:sample=>(quantity=nothing,unit=nothing)) :
            (:frequency=>(quantity=Units.Quantity{:frequency}(),unit=coordinates.frequency_unit))
        metadata!(table,"observation_columns",(; (first_column,columns...)...);style=:note)
    elseif coordinates.kind===:assemblies
        table=DataFrame(frequency=coordinates.frequencies)
        columns=Pair{Symbol,Any}[:frequency=>(quantity=Units.Quantity{:frequency}(),unit=coordinates.frequency_unit)]
        for (index,assembly) in enumerate(coordinates.assemblies)
            name=Symbol(coordinates.labels[assembly])
            name===:frequency && throw(ArgumentError("assembly name conflicts with the frequency column"))
            table[!,name]=[scalar ? values : values[index]]
            push!(columns,name=>(quantity=product.quantity,unit=product.unit,assembly))
        end
        metadata!(table,"observation_columns",(;columns...);style=:note)
    elseif coordinates.kind===:samples
        if haskey(coordinates,:rows)
            dims=(length(coordinates.rows),length(coordinates.columns),length(coordinates.samples),length(coordinates.trials))
            shaped=reshape(scalar ? [values] : values,dims...)
            table=DataFrame([(frequency=coordinates.frequencies[k],row=coordinates.rows[i],column=coordinates.columns[j],
                trial=coordinates.trials[t],value=shaped[i,j,k,t]) for t in 1:dims[4] for k in 1:dims[3] for i in 1:dims[1] for j in 1:dims[2]])
        else
            shaped=reshape(scalar ? [values] : values,length(coordinates.assemblies),length(coordinates.trials))
            table=DataFrame([(assembly=coordinates.assemblies[i],trial=coordinates.trials[t],value=shaped[i,t])
                for t in eachindex(coordinates.trials) for i in eachindex(coordinates.assemblies)])
        end
        metadata!(table,"observation_columns",(value=(quantity=product.quantity,unit=product.unit),);style=:note)
    elseif values isa NamedTuple
        table=DataFrame(values)
        metadata!(table,"observation_columns",(;);style=:note)
    else
        vector=vec(scalar ? [values] : values)
        table=DataFrame(index=collect(eachindex(vector)),value=copy(vector))
        metadata!(table,"observation_columns",(value=(quantity=product.quantity,unit=product.unit),);style=:note)
    end
    coordinate_columns=coordinates.kind===:samples ? Tuple(filter(!=(:value),propertynames(table))) :
        coordinates.kind in (:matrix,:diagonal,:vector,:assemblies,:array) ? (first(propertynames(table)),) : ()
    metadata!(table,"coordinate_columns",coordinate_columns;style=:note)
    metadata!(table,"coordinates",Grammar.detach(coordinates);style=:note)
    metadata!(table,"quantity",product.quantity;style=:note)
    metadata!(table,"request",Grammar.detach(product.request);style=:note)
    metadata!(table,"statistic",Grammar.detach(product.statistic);style=:note)
    metadata!(table,"family",product.family;style=:note)
    metadata!(table,"unit",product.unit;style=:note)
    metadata!(table,"basis",product.basis;style=:note)
    metadata!(table,"missing_reason",Grammar.detach(product.missing_reason);style=:note)
    gridpoint_id===nothing || metadata!(table,"gridpoint_id",gridpoint_id;style=:note)
    return table
end

"""
$(TYPEDSIGNATURES)

Build one table per retained quantity, grouped by the owning physical family.
Full matrices retain every coefficient in row-major order. Each row represents
one retained frequency or sample coordinate.
"""
function _quantity_tables(products; gridpoint_id=nothing)
    allunique((q.family,_quantity_name(q)) for q in products) || throw(ArgumentError(
        "multiple retained products share a quantity name; select a complete request with tabulate(observed, request)"))
    families=unique(q.family for q in products)
    return (;(family=>(;(_quantity_name(q)=>_quantity_table(q;gridpoint_id) for q in products if q.family==family)...)
        for family in families)...)
end
tabulate(observed::ObservedResult) = _quantity_tables(observed.quantities;gridpoint_id=observed.gridpoint.id)

"""Build one table from a retained quantity request."""
tabulate(observed::ObservedResult,request) = _quantity_table(Grammar.observation_product(observed,request);gridpoint_id=observed.gridpoint.id)
tabulate(observed::AbstractVector{<:ObservedResult}) = map(tabulate,observed)

select(::CableConstantsTableDefinition,observed::ObservedResult;reference=nothing) = observed.quantities
function select(definition::LineParametersTableDefinition,observed::ObservedResult;reference=nothing)
    return [Grammar.observation_product(observed,request)
        for request in Grammar.observation_requests(observed,definition.requests).retained]
end
function tabulate(::Union{TableReportDefinition,CableConstantsTableDefinition,LineParametersTableDefinition},
        observed,selected;reference=nothing)
    observed isa ObservedResult && return _quantity_tables(selected;gridpoint_id=observed.gridpoint.id)
    return map((point,products) -> _quantity_tables(products;gridpoint_id=point.gridpoint.id),
        observed,selected)
end
function report(definition::CableConstantsTableDefinition,source::Engine.CableConstants;kwargs...)
    return report(definition,ObservedResult(source;clip=definition.clip,kwargs...))
end
function report(definition::LineParametersTableDefinition,source::Engine.LineParameters;kwargs...)
    return report(definition,ObservedResult(source,definition.requests;complete_pairs=true,clip=definition.clip,
        frequency_unit=definition.frequency_unit,length_unit=definition.length_unit,
        quantity_units=definition.quantity_units,kwargs...))
end

"""
$(TYPEDSIGNATURES)

Reject aggregate conversion because an observation contains separate physical
quantities. Select a quantity table with `ReportBuilder.tabulate(observed, R)`
or a leaf such as `ReportBuilder.tabulate(observed).Z.R`.
"""
function DataFrame(::ObservedResult)
    throw(ArgumentError("ObservedResult contains separate quantity tables; use ReportBuilder.tabulate(observed, R) or ReportBuilder.tabulate(observed).Z.R"))
end
