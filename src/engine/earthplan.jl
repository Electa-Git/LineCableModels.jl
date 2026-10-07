"""
$(TYPEDEF)

Hold the earth calculations of a line-parameter computation, fixed before the frequency loop.
The workspace stores it as `plan.earth`, beside its arrays in `buffers.earth`. An earth
calculation is a record with these fields.

- `impedance` and `admittance`: each is `(formula, pairs)`, the formula of that quantity and
  the conductor pairs whose values it publishes. It is `nothing` when the calculation does
  not serve the quantity.
- `earth`: the layered earth model or the reduction that its expressions consume.
- `reductions`: on a reduced earth, each is `(expression, options, pairs)`, one expression of
  the reduction, its normalized options and the conductor pairs it serves. It is `nothing` on
  the layered earth.
- `parts`: each is `(expressions, options, pairs)`, the expressions of one part of the
  formula, their normalized options and the conductor pairs they serve.
- `pairs`: each is `(pair, physical, reuse_from)`. It holds the pair on the decided earth and
  the physical pair. `reuse_from` is the earlier pair with the same inputs whose computed value
  this pair takes when their media agree, or its own index when the pair is computed.
- `media`: the medium of each conductor on the decided earth.

The plan of one slot serves only that slot's quantity. The merge of the impedance and the
admittance plans gives one calculation both quantities when it can serve them.

$(TYPEDFIELDS)
"""
struct EarthPlan
    "Earth calculations. Their record types differ, and the frequency loop dispatches on them once."
    calculations::Tuple
end

"""
$(TYPEDSIGNATURES)

Build the calculation of the earth formula `formula` for the physical conductor pairs
`physical[indices]`, on the decided `earth`: the layered earth `model`, or a reduction of it.
The constructor checks each pair and the geometric restrictions of each expression. It checks
that `formula` defines its expressions for the layers of `model`. On a reduction, it checks
that the reduction admits each expression and defines its own expression for each physical
pair. The distinct expressions, with the options that `formula` declares for each, become the
parts.
`geometry` holds the conductor radii and the physical layer of each conductor.

A formula that computes the whole system adds a method on its own type.
"""
function EarthPlan(formula::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation},
        earth, model::EarthModel, physical::AbstractVector{<:EarthPair}, indices,
        geometry::NamedTuple)
    earth === model && validate(model, formula)
    decided = [EarthPair(physical[index], earth) for index in indices]
    expressions = map(decided) do pair
        validate(pair)
        expression = Expression(formula, pair)
        validate(pair, expression)
        expression
    end
    projected = formulation_options(formula, expressions)
    for expression in projected.expressions
        validate(expression, model.layers)
    end
    # The reduction's own expressions, their options and the physical pairs each serves.
    reductions = if earth === model
        nothing
    else
        rule = EquivalentHomogeneous.rule(earth)
        foreach(expression -> validate(rule, expression), expressions)
        reduced = [Expression(rule, physical[index]) for index in indices]
        normalized = formulation_options(rule, reduced)
        foreach(validate, normalized.expressions)
        Tuple(map(normalized.expressions, normalized.options) do expression, options
            (; expression, options, pairs = findall(==(expression), reduced))
        end)
    end
    # The pairs that each distinct expression serves.
    served = [findall(==(expression), expressions) for expression in projected.expressions]
    # A pair takes the value of the nearest earlier pair of its part with the same inputs.
    reuse_from = collect(eachindex(decided))
    for members in served, ordinal in eachindex(members)
        position = members[ordinal]
        for earlier in (ordinal - 1):-1:1
            candidate = members[earlier]
            if same_physical_state(formula, decided[position], decided[candidate], geometry)
                reuse_from[position] = candidate
                break
            end
        end
    end
    parts = map(projected.expressions, projected.options, served) do expression, options, pairs
        (expressions = (expression,), options, pairs)
    end
    pairs = [(pair = decided[position], physical = physical[index],
                 reuse_from = reuse_from[position]) for (position, index) in enumerate(indices)]
    media = earth === model ? geometry.layers :
            [layer == 1 ? 1 : 2 for layer in geometry.layers]
    published = (; formula, pairs = collect(eachindex(decided)))
    impedance = formula isa EarthImpedanceFormulation ? published : nothing
    admittance = formula isa EarthAdmittanceFormulation ? published : nothing
    calculation = (; impedance, admittance, earth, reductions, parts = Tuple(parts), pairs,
        media)
    return EarthPlan((calculation,))
end

"""
$(TYPEDSIGNATURES)

Build the earth slot that holds the formula `formula`, for the physical conductor pairs
`physical` of an earth `model`. The slot decides its earth once: an explicit
`equivalent_earth` applies. Without one, a formula that does not admit any layer from 3 to N
consumes the `:default` reduction of a model with N > 2 layers, and every other formula sees
the layered earth.
"""
function EarthPlan(formula::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation},
        model::EarthModel, physical::AbstractVector{<:EarthPair}, geometry::NamedTuple)
    layers = length(model.layers)
    earth = if formula.equivalent_earth !== nothing
        formula.equivalent_earth
    elseif layers > 2 &&
           !any(k -> _admits_layer(Expression(formula, first(physical)), k), 3:layers)
        EquivalentHomogeneous.AbstractSequence(LineCableModels.formula(:default))
    else
        model
    end
    return EarthPlan(formula, earth, model, physical, eachindex(physical), geometry)
end

"""
$(TYPEDSIGNATURES)

Build the earth slot that holds a recipe of homogeneous formulas, for the physical conductor
pairs `physical` of an earth `model`. The slot decides its earth once: the explicit
reduction that its formulas agree on, or else the `:default` reduction on more than two
layers. Each pair then takes the recipe's formula for its layers on that earth, and each
distinct formula builds one calculation.

# Errors

- Throws `ArgumentError` when a formula of the recipe admits a layer above 2, or when the
  formulas give different equivalent earths.
"""
function EarthPlan(recipe::NamedTuple, model::EarthModel,
        physical::AbstractVector{<:EarthPair}, geometry::NamedTuple)
    selected = Tuple(leaf for leaf in recipe if leaf !== nothing)
    layers = length(model.layers)
    for leaf in selected
        expression = Expression(leaf, first(physical))
        any(k -> _admits_layer(expression, k), 3:max(3, layers)) && throw(ArgumentError(
            "formula :$(formula_id(leaf)) is multilayer and is used alone"))
    end
    explicit = unique(leaf.equivalent_earth for leaf in selected
        if leaf.equivalent_earth !== nothing)
    length(explicit) > 1 && throw(ArgumentError(
        "the formulas of an earth recipe give different equivalent earths"))
    earth = !isempty(explicit) ? only(explicit) :
            layers > 2 ? EquivalentHomogeneous.AbstractSequence(LineCableModels.formula(:default)) :
            model
    leaves = [Expression(recipe, EarthPair(pair, earth)).selection for pair in physical]
    calculations = NamedTuple[]
    for leaf in unique(leaves)
        indices = findall(value -> value === leaf, leaves)
        append!(calculations,
            EarthPlan(leaf, earth, model, physical, indices, geometry).calculations)
    end
    return EarthPlan(Tuple(calculations))
end

"""
$(TYPEDSIGNATURES)

Merge the earth plans of the impedance and the admittance slot. The first plan's calculations
compute the impedance and the second plan's the admittance. An impedance and an admittance
calculation become one calculation serving both when their formulas have the same physical
state and their parts the same options.
"""
function EarthPlan(impedance::EarthPlan, admittance::EarthPlan)
    calculations = NamedTuple[]
    remaining = collect(NamedTuple, admittance.calculations)
    for calculation in impedance.calculations
        index = findfirst(remaining) do candidate
            same_physical_state(calculation.impedance.formula, candidate.admittance.formula) &&
                same_physical_state(first(calculation.parts).options.data,
                    first(candidate.parts).options.data)
        end
        if index === nothing
            push!(calculations, calculation)
        else
            push!(calculations,
                merge(calculation, (admittance = remaining[index].admittance,)))
            deleteat!(remaining, index)
        end
    end
    append!(calculations, remaining)
    return EarthPlan(Tuple(calculations))
end

"""
$(TYPEDSIGNATURES)

Whether `formula` gives the conductor pairs `a` and `b` the same value when their media
agree. `geometry` holds the conductor radii. By default the pairs' destination indices take
part, so distinct pairs are computed separately. A formula whose arithmetic does not read
them compares the inputs it reads instead.
"""
function same_physical_state(
        ::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation},
        a::EarthPair, b::EarthPair, geometry::NamedTuple)
    return same_physical_state((a.row, a.column), (b.row, b.column))
end

"""
$(TYPEDSIGNATURES)

Extend `buffers` with `earth`, the arrays of the earth calculations of `earth`.

- `static` and `evaluated` hold the earth's layer properties, static and at each frequency.
  `evaluated` is `nothing` unless a calculation consumes the evaluated layers.
- `calculations` holds, for each calculation, the `rho`, `epsilon` and `mu` of the media that
  each of its conductor pairs sees, one column per pair, and the layer thicknesses of a
  layered earth.
- `pairs` holds, for the conductor pairs, the pair whose computed value each one took at the
  current frequency, the ranges of its integrals and warnings, and the warning records.

The integration warnings, when a formula integrates, share the warning records: the
formulas of the calculations provision their own arrays first, through the formulation's
methods.
"""
function initialize_buffers(earth::EarthPlan, ::Type{T}, input, plan, buffers) where {T}
    layers = input.earth.layers
    static = (rho = collect(getproperty.(layers, :rho)),
        eps_r = collect(getproperty.(layers, :eps_r)),
        mu_r = collect(getproperty.(layers, :mu_r)))
    needed = any(calculation -> !(calculation.earth isa EquivalentHomogeneous.BeforeFD),
        earth.calculations)
    evaluated = needed ?
                map(_ -> Matrix{eltype(input.freq)}(undef, length(layers), input.n_frequencies),
        static) : nothing
    Evaluated = NamedTuple{(:rho, :eps_r, :mu_r), NTuple{3, Matrix{eltype(input.freq)}}}
    # The layered earth keeps every physical layer, and its thicknesses when it has interior
    # layers. A reduction keeps air and one equivalent medium.
    calculations = map(earth.calculations) do calculation
        count = calculation.earth isa EarthModel ? length(layers) : 2
        columns = length(calculation.pairs)
        thickness = count > 2 ? T[layer.thickness for layer in layers] : nothing
        (rho = Matrix{T}(undef, count, columns), epsilon = Matrix{T}(undef, count, columns),
            mu = Matrix{T}(undef, count, columns), thickness)
    end
    pairs = (representatives = zeros(Int, input.n_cables^2),
        integral_ranges = Vector{UnitRange{Int}}(undef, input.n_cables^2),
        warning_ranges = Vector{UnitRange{Int}}(undef, input.n_cables^2),
        warnings = NamedTuple[])
    # The record type does not depend on the calculations, nor on whether the evaluated
    # layers are needed.
    Record = NamedTuple{(:static, :evaluated, :calculations, :pairs),
        Tuple{typeof(static), Union{Nothing, Evaluated}, Tuple, typeof(pairs)}}
    buffers = merge(buffers, (earth = Record((static, evaluated, calculations, pairs)),))
    haskey(buffers, :quadrature) || return buffers
    return merge(buffers, (quadrature = merge(buffers.quadrature,
        (warnings = buffers.earth.pairs.warnings,)),))
end
