"""
$(TYPEDEF)

Select retained UQ statistical products for quantity-wise tables.

$(TYPEDFIELDS)
"""
struct MonteCarloTableDefinition{U} <: AbstractReportDefinition
    "Length prefix used during raw-input observation construction."
    length_unit::Symbol
    "Optional quantity-unit overrides."
    quantity_units::U
    "Engineering recentering for primary quantities."
    clip::Bool
end
MonteCarloTableDefinition(length_unit::Symbol=:kilo,quantity_units=nothing) =
    MonteCarloTableDefinition(length_unit,quantity_units,true)

function report(definition::MonteCarloTableDefinition,source::AbstractUncertaintyResult)
    requests=Tuple((UQ.statistics,selector,statistic) for selector in (R,L,G,C) for statistic in (Statistics.mean,Statistics.std))
    observed=observables(source,requests;length_unit=definition.length_unit,
        quantity_units=definition.quantity_units,clip=definition.clip)
    return report(definition,observed)
end
function tabulate(::MonteCarloTableDefinition,observed;reference=nothing)
    return tabulate(observed)
end

function _sampling_tables(points,reference)
    sampling=NamedTuple[];precision=NamedTuple[];statistics=NamedTuple[]
    operands=reference===nothing ? ((:candidate,points),) : ((:candidate,points),(:reference,[reference]))
    for (role,observations) in operands,point in observations
        id=point.gridpoint.id
        for quantity in point.quantities
            quantity.family===:statistics || continue
            push!(statistics,(role,identity=id,request=quantity.request,statistic=quantity.statistic,
                unit=quantity.unit,values=Grammar.detach(quantity.values)))
        end
        record=get(point.gridpoint,:sampling,nothing)
        record===nothing && continue
        scalar=(; (key=>(value===nothing ? missing : value isa Union{Number,Bool,Symbol,AbstractString,Missing} ? value : string(value))
            for (key,value) in pairs(record) if !(key in (:mean_standard_error,:frequencies,:basis)))...)
        push!(sampling,merge(scalar,(role,point=id===nothing ? record.point : id.problem_index,
            formulation_index=id===nothing ? missing : id.formulation_index,
            method=only(Grammar.observation_labels([point])),std_sampling_precision=missing)))
        for (quantity,values) in pairs(record.mean_standard_error),index in CartesianIndices(values)
            coordinates=Tuple(index)
            unit=Units.native_unit(Units.quantity(getfield(Engine,quantity)),record.basis)
            push!(precision,(role,identity=id,point=id===nothing ? record.point : id.problem_index,
                formulation_index=id===nothing ? missing : id.formulation_index,
                method=only(Grammar.observation_labels([point])),quantity,index=coordinates,
                row=first(coordinates),column=length(coordinates)==3 ? coordinates[2] : missing,
                frequency_Hz=Grammar.nominal(record.frequencies[length(coordinates)==3 ? last(coordinates) : 1]),
                standard_error=values[index],unit=Units.label(unit)))
        end
    end
    return (statistics=DataFrame(statistics),sampling=DataFrame(sampling),mean_sampling_precision=DataFrame(precision))
end
