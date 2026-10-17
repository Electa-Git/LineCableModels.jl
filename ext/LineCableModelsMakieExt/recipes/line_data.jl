# Renderers select detached records. Scientific normalization and extraction
# belong to ObservedResult construction.
function _line_plot_data(records::Tuple)
    for q in records
        kind=q.coordinates.kind
        kind in (:matrix,:diagonal,:vector,:assemblies,:array) || throw(ArgumentError(
            "ordinary plot does not select an interpretation for retained $kind products; use a retained mean/std request or an explicit statistical plotting verb"))
        kind===:array && ndims(vcat(q.values))>1 && throw(ArgumentError(
            "generic arrays of rank greater than one need owner-provided coordinate interpretation"))
    end
    coordinates=map(records) do product
        c=product.coordinates
        c.kind in (:matrix,:diagonal,:vector) ?
            (c.kind===:vector ? c.positions : c.rows,
                c.kind===:matrix ? c.columns : [1],c.samples) :
            ([1],[1],collect(1:length(vcat(product.values))))
    end
    xdata=map(records) do product
        c=product.coordinates
        if c.kind===:assemblies
            return (values=string.(c.labels[c.assemblies]),quantity=Units.Quantity{:dimensionless}(),
                unit=Units.units(:base,:dimensionless),label="Assembly")
        elseif c.kind===:array
            n=length(vcat(product.values))
            indices=get(c,:indices,())
            selected=length(indices)==1 ? only(indices) : Colon()
            original=selected isa Integer ? [selected] : selected isa AbstractVector || selected isa AbstractRange ? collect(selected) : nothing
            original===nothing || length(original)==n || throw(DimensionMismatch("retained element indices do not match values"))
            return (values=original===nothing ? collect(1:n) : original,
                quantity=Units.Quantity{:dimensionless}(),unit=Units.units(:base,:dimensionless),
                label=original===nothing ? "Retained element position" : "Original element index")
        end
        for (values_key,unit_key,quantity) in (
                (:frequencies,:frequency_unit,Units.Quantity{:frequency}()),
                (:times,:time_unit,Units.Quantity{:time}()))
            haskey(c,values_key) && c[values_key]!==nothing || continue
            return (values=c[values_key],quantity,unit=c[unit_key])
        end
        throw(ArgumentError(
            "a plot requires the frequency or time axis retained during observation"))
    end
    observations=map(records,coordinates) do product,indices
        # One value is an array of one.
        merge(product,(values=reshape(vcat(product.values),length.(indices)...),))
    end
    resolutions=map(q -> (clip=q.clipped,kind=q.thresholds===nothing ? :unassessed : q.thresholds.kind),records)
    return (;xdata,observations,coordinates,resolutions)
end
