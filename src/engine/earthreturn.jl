function earth!(workspace::LineParametersWorkspace, frequency::Int,
        calculations::Tuple = workspace.plan.earth.calculations,
        materials::Tuple = workspace.buffers.earth.calculations)
    foreach(calculations, materials) do calculation, material
        earth!(calculation, material, workspace, frequency)
    end
    return workspace
end

# Build the calculation's Functor at the frequency, evaluate its parts into the destinations,
# convert them when the formula requires it, then publish each served quantity's pairs.
function earth!(calculation::NamedTuple, materials::NamedTuple, workspace, frequency::Int)
    formula = something(calculation.impedance, calculation.admittance).formula
    functor = Functor(formula, (; jω = workspace.input.jω[frequency], frequency,
        materials.rho, materials.epsilon, materials.mu, materials.thickness,
        calculation.media); workspace)
    earth!(functor.input.destinations, calculation, functor, workspace)
    physical = earth!(formula, functor, workspace)
    if calculation.impedance !== nothing && haskey(physical, :impedance)
        destination = workspace.buffers.Zearth
        for index in calculation.impedance.pairs
            pair = calculation.pairs[index].pair
            destination[pair.row, pair.column] = physical.impedance[pair.row, pair.column]
        end
    end
    if calculation.admittance !== nothing && haskey(physical, :admittance)
        destination = workspace.buffers.Pearth
        for index in calculation.admittance.pairs
            pair = calculation.pairs[index].pair
            destination[pair.row, pair.column] = physical.admittance[pair.row, pair.column]
        end
    end
    return workspace
end

# The parts of a formula that does not solve a whole system wrote the physical coefficients
# into `Zearth` or `Pearth`, so the conversion does not return any matrix.
earth!(::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation}, ::Functor, workspace) = (;)

"""
$(TYPEDSIGNATURES)

Evaluate the parts of an earth calculation at the frequency of `functor` and distribute their
scalar coefficients into the aligned matrices `destinations`. A part gives its expressions and
options. The calculation's pairs give each pair's source-target geometry and the earlier pair
with the same inputs, `reuse_from`.

A pair takes the value of that earlier pair only when their current media agree under
`same_physical_state`. `buffers.earth.pairs` records, at each frequency, the pair whose
computed value each pair took. Every logical integral and warning retains its receiving row
and source column.
"""
function earth!(destinations::Tuple{Vararg{AbstractMatrix}}, calculation::NamedTuple,
        functor::Functor, workspace)
    computed = workspace.buffers.earth.pairs
    length(computed.representatives) >= length(calculation.pairs) ||
        throw(DimensionMismatch("earth interaction scratch is too small"))
    empty!(computed.warnings)
    foreach(calculation.parts) do part
        earth!(destinations, part, calculation, functor, workspace)
    end
    empty!(computed.warnings)
    return destinations
end

function earth!(destinations, part, calculation, functor, workspace)
    computed = workspace.buffers.earth.pairs
    observations = workspace.buffers.observations
    materials = (functor.input.rho, functor.input.epsilon, functor.input.mu)
    for index in part.pairs
        entry = calculation.pairs[index]
        pair = entry.pair
        # Follow the earlier pairs with the same inputs until one whose media agree.
        current = index
        source = entry.reuse_from
        while source != current
            same = let source = source
                all(
                    values -> same_physical_state(@view(values[:, index]),
                        @view(values[:, source])),
                    materials)
            end
            same && break
            current = source
            source = calculation.pairs[source].reuse_from
        end
        if source == current
            computed.representatives[index] = index
            first_integral = observations === nothing ? 1 : length(observations)+1
            first_warning = length(computed.warnings)+1
            point = Functor(functor, (; pair, entry.physical,
                rho = @view(functor.input.rho[:, index]),
                epsilon = @view(functor.input.epsilon[:, index]),
                mu = @view(functor.input.mu[:, index]), part.options))
            validate(point.input, functor.formula)
            values = map(expression -> expression(point, workspace), part.expressions)
            all(value -> value isa Number && isfinite(value), values) || throw(DomainError(
                values, "earth coefficients must be finite scalars"))
            converted = map(value -> oftype(functor.input.jω, value), values)
            length(converted) == length(destinations) ||
                throw(DimensionMismatch("earth coefficients must match their destinations"))
            foreach(destinations, converted) do destination, value
                destination[pair.row, pair.column] = value
            end
            computed.integral_ranges[index] = first_integral:(observations === nothing ? 0 :
                                                           length(observations))
            computed.warning_ranges[index] = first_warning:length(computed.warnings)
        else
            representative = computed.representatives[source]
            computed.representatives[index] = representative
            previous_pair = calculation.pairs[representative].pair
            for destination in destinations
                destination[pair.row, pair.column] = destination[previous_pair.row, previous_pair.column]
            end
            context = (receiver = pair.row, source = pair.column)
            if observations !== nothing
                for position in computed.integral_ranges[representative]
                    record = observations[position]
                    integral_context = record.context isa NamedTuple ?
                                       merge(record.context, context) : record.context
                    push!(observations, merge(record, (context = integral_context,)))
                end
            end
            for position in computed.warning_ranges[representative]
                record = computed.warnings[position]
                integral_context = record.context isa NamedTuple ?
                                   merge(record.context, context) : record.context
                record_integral!(nothing, nothing, record.value, record.estimated_error,
                    record.controls, integral_context)
            end
        end
    end
    return destinations
end
