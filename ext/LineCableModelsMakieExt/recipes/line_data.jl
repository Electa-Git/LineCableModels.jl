# Renderers select detached records. Scientific normalization and extraction
# belong to ObservedResult construction.
function _prepare_line_observations(records::Tuple)
    for q in records
        kind=q.coordinates.kind
        kind in (:matrix,:diagonal,:assemblies,:array) || throw(ArgumentError(
            "ordinary plot does not select an interpretation for retained $kind products; use a retained mean/std request or an explicit statistical plotting verb"))
        kind===:array && q.values isa AbstractArray && ndims(q.values)>1 && throw(ArgumentError(
            "generic arrays of rank greater than one need owner-provided coordinate interpretation"))
    end
    coordinates=map(records) do product
        c=product.coordinates
        c.kind in (:matrix,:diagonal) ?
            (c.rows,c.kind===:diagonal ? [1] : c.columns,c.samples) :
            ([1],[1],collect(1:length(product.values isa Number || ismissing(product.values) ? [product.values] : product.values)))
    end
    frequencies=map(records) do product
        c=product.coordinates
        if c.kind===:assemblies
            return (values=string.(c.labels[c.assemblies]),quantity=Units.Quantity{:dimensionless}(),
                unit=Units.units(:base,:dimensionless),label="Assembly")
        elseif c.kind===:array
            n=product.values isa Number || ismissing(product.values) ? 1 : length(product.values)
            indices=get(c,:indices,())
            selected=length(indices)==1 ? only(indices) : Colon()
            original=selected isa Integer ? [selected] : selected isa AbstractVector || selected isa AbstractRange ? collect(selected) : nothing
            original===nothing || length(original)==n || throw(DimensionMismatch("retained element indices do not match values"))
            return (values=original===nothing ? collect(1:n) : original,
                quantity=Units.Quantity{:dimensionless}(),unit=Units.units(:base,:dimensionless),
                label=original===nothing ? "Retained element position" : "Original element index")
        end
        c.frequencies===nothing && throw(ArgumentError(
            "a frequency plot requires a frequency axis retained during observation"))
        (values=c.frequencies,quantity=Units.Quantity{:frequency}(),unit=c.frequency_unit)
    end
    observations=map(records,coordinates) do product,indices
        merge(product,(values=reshape(product.values isa Number || ismissing(product.values) ?
            [product.values] : product.values,length.(indices)...),))
    end
    resolutions=map(q -> (clip=q.clipped,kind=q.thresholds===nothing ? :unassessed : q.thresholds.kind),records)
    return (;frequencies,observations,coordinates,resolutions)
end
