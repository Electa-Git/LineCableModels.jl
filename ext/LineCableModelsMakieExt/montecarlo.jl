function _observed_distribution(result::LineCableModels.MonteCarloResult,request;
        point::Integer=1,bins=nothing,kwargs...)
    selector=request isa Function ? request : first(request)
    indices=request isa Function ? (1,) : Base.tail(request)
    retained=(LineCableModels.histograms,selector,indices...,bins)
    return Grammar.ObservedResult(result,point,(retained,);kwargs...)
end

function _distribution_product(observed,request)
    products=filter(q -> q.coordinates.kind===:histogram,observed.quantities)
    if request!==nothing
        products=filter(q -> q.request==request || request_identity(q.request)==request_identity(request) ||
            request isa Function && request_identity(q.request)[2]===request,products)
    end
    length(products)==1 || throw(ArgumentError("select one retained histogram marginal"))
    return only(products)
end

function _distribution_plot(observed::Grammar.ObservedResult,request,kind;
        normalization=:none,qqline=:identity,title=nothing,fig_size=(800,400),
        backend=nothing,display_plot=true,controls=true,export_theme=:default,open_export=true,kwargs...)
    for key in (:point,:bins,:clip,:atol,:units,:length_unit,:quantity_units,:frequency_unit,:freq_unit,:frequencies,:complete_pairs)
        haskey(kwargs,key) && throw(ArgumentError("$key belongs to observation construction; this statistical view consumes retained products"))
    end
    product=_distribution_product(observed,request)
    distribution=product.distribution
    if kind===:qq
        qqline in (:identity,:none) || throw(ArgumentError("qqline must be :identity or :none"))
        pairs=distribution.qq
        pairs===nothing && throw(ArgumentError("Q–Q coordinates were not retained"))
        x=pairs.sample;y=pairs.model
        ordinate=(values=y,quantity=product.quantity,unit=product.unit)
    elseif kind===:empirical
        empirical=distribution.empirical_cdf
        empirical===nothing && throw(ArgumentError("empirical CDF coordinates were not retained"))
        x=empirical.x;y=empirical.y
        ordinate=(values=y,quantity=Units.Quantity{:cumulative_probability}(),unit=Units.UnitExpr())
    elseif kind===:cdf
        x=distribution.model_cdf.x;y=distribution.model_cdf.y
        ordinate=(values=y,quantity=Units.Quantity{:cumulative_probability}(),unit=Units.UnitExpr())
    else
        normalization in (:none,:probability,:pdf) || throw(ArgumentError("normalization must be :none, :probability or :pdf"))
        field=kind===:density || normalization===:pdf ? :density : normalization===:probability ? :probability : :count
        y=getproperty(product.values,field)
        any(ismissing,y) && throw(ArgumentError("histogram counts were not retained"))
        x=distribution.edges[1:end-1]
        quantity=field===:density ? Units.Quantity{:probability_density}() : field===:count ?
            Units.Quantity{:sample_count}() : Units.Quantity{:probability}()
        unit=field===:density ? product.ordinate_units.density : Units.UnitExpr()
        ordinate=(values=y,quantity,unit)
    end
    abscissa=(values=x,quantity=product.quantity,unit=product.unit)
    heading=title===nothing ? "$(Units.symbol(product.quantity)) $(kind)" : String(title)
    return _addon_statistical_plot(abscissa,ordinate;title=heading,fig_size,backend,
            display_plot,controls,export_theme,open_export,kwargs...) do axis,groups,order,labels,series
        plot=if kind===:histogram
            edges=distribution.edges
            Makie.barplot!(axis,(edges[1:end-1].+edges[2:end])./2,y;width=diff(edges),gap=0,label=heading)
        elseif kind===:density
            stairs!(axis,distribution.edges,vcat(y,last(y));step=:post,label=heading)
        elseif kind===:empirical
            stairs!(axis,x,y;step=:post,label=heading)
        elseif kind===:qq
            scatter!(axis,x,y;label=heading)
        else
            lines!(axis,x,y;label=heading)
        end
        groups[:distribution]=Any[plot];push!(order,:distribution);labels[:distribution]=heading
        push!(series,(xdata=x,ydata=y,plots=Any[plot]))
        if kind===:qq && qqline===:identity
            limits=collect(distribution.qq.reference)
            identity=lines!(axis,limits,limits;color=:black,linestyle=:dash,label="identity")
            groups[:identity]=Any[identity];push!(order,:identity);labels[:identity]="identity"
        end
    end
end

Makie.hist(observed::Grammar.ObservedResult,request=nothing;kwargs...) = _distribution_plot(observed,request,:histogram;kwargs...)
Makie.stairs(observed::Grammar.ObservedResult,request=nothing;kwargs...) = _distribution_plot(observed,request,:density;kwargs...)
Makie.ecdfplot(observed::Grammar.ObservedResult,request=nothing;kwargs...) = _distribution_plot(observed,request,:empirical;kwargs...)
Makie.lines(observed::Grammar.ObservedResult,request=nothing;kwargs...) = _distribution_plot(observed,request,:cdf;kwargs...)
Makie.qqplot(observed::Grammar.ObservedResult,request=nothing;kwargs...) = _distribution_plot(observed,request,:qq;kwargs...)

function Makie.hist(source::LineCableModels.MonteCarloResult,request=LineCableModels.R;
        point=1,bins=nothing,kwargs...)
    acquisition,presentation=_plot_observation_options(kwargs)
    observed=_observed_distribution(source,request;point,bins,acquisition...)
    return Makie.hist(observed;presentation...)
end
function Makie.stairs(source::LineCableModels.MonteCarloResult,request=LineCableModels.R;
        point=1,bins=nothing,kwargs...)
    acquisition,presentation=_plot_observation_options(kwargs)
    observed=_observed_distribution(source,request;point,bins,acquisition...)
    return Makie.stairs(observed;presentation...)
end
function Makie.ecdfplot(source::LineCableModels.MonteCarloResult,request=LineCableModels.R;
        point=1,bins=nothing,kwargs...)
    acquisition,presentation=_plot_observation_options(kwargs)
    observed=_observed_distribution(source,request;point,bins,acquisition...)
    return Makie.ecdfplot(observed;presentation...)
end
function Makie.lines(source::LineCableModels.MonteCarloResult,request=LineCableModels.R;
        point=1,bins=nothing,kwargs...)
    acquisition,presentation=_plot_observation_options(kwargs)
    observed=_observed_distribution(source,request;point,bins,acquisition...)
    return Makie.lines(observed;presentation...)
end
function Makie.qqplot(source::LineCableModels.MonteCarloResult,request=LineCableModels.R;
        point=1,bins=nothing,kwargs...)
    acquisition,presentation=_plot_observation_options(kwargs)
    observed=_observed_distribution(source,request;point,bins,acquisition...)
    return Makie.qqplot(observed;presentation...)
end
