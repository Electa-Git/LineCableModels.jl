# Renderers select detached records. Scientific normalization and extraction
# belong to ObservedResult construction.
_line_request_family(request) = Units.family(request_quantity(request))
_diagonal_request(request) = request_identity(request) isa Tuple && diag in request_identity(request)
_family_parent(::Val{:series}) = Z
_family_parent(::Val{:shunt}) = Y

function _prepare_line_observations(observed::Grammar.ObservedResult;ydata,kwargs...)
    records=map(request -> Grammar.observation_product(observed,request),ydata)
    all(q -> q.coordinates.kind in (:matrix,:diagonal),records) ||
        throw(ArgumentError("matrix plots require retained matrix or diagonal quantities"))
    coordinates=map(records) do product
        c=product.coordinates
        (c.rows,c.kind===:diagonal ? [1] : c.columns,c.samples)
    end
    first_coordinate=first(records).coordinates
    all(q -> isequal(q.coordinates.frequencies,first_coordinate.frequencies) &&
        q.coordinates.frequency_unit==first_coordinate.frequency_unit,records) ||
        throw(DimensionMismatch("quantities on one dashboard must have identical retained frequency coordinates"))
    first_coordinate.frequencies===nothing && throw(ArgumentError(
        "a frequency plot requires a frequency axis retained during observation"))
    frequency=(values=first_coordinate.frequencies,quantity=Units.Quantity{:frequency}(),unit=first_coordinate.frequency_unit)
    observations=map(records,coordinates) do product,indices
        merge(product,(values=reshape(product.values isa Number || ismissing(product.values) ?
            [product.values] : product.values,length.(indices)...),))
    end
    resolutions=map(q -> (clip=q.clipped,kind=q.thresholds===nothing ? :unassessed : q.thresholds.kind),records)
    return (;frequency,observations,coordinates,resolutions)
end
