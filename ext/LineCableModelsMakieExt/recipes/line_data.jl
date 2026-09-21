# Renderers select detached records. Scientific normalization and extraction
# belong to ObservedResult construction.
_line_request_family(request) = Units.family(request_quantity(request))
_diagonal_request(request) = request_identity(request) isa Tuple && diag in request_identity(request)
_family_parent(::Val{:series}) = Z
_family_parent(::Val{:shunt}) = Y

function _prepare_line_observations(records::Tuple)
    all(q -> q.coordinates.kind in (:matrix,:diagonal),records) ||
        throw(ArgumentError("matrix plots require retained matrix or diagonal quantities"))
    coordinates=map(records) do product
        c=product.coordinates
        (c.rows,c.kind===:diagonal ? [1] : c.columns,c.samples)
    end
    frequencies=map(records) do product
        c=product.coordinates
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
