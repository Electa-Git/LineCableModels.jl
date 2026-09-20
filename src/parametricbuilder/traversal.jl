function _formulations(inner::AbstractFormulation)
    return typeof(inner)[inner]
end

function _formulations(source::Gridspace{Target}) where {Target}
    formulations = Vector{Target}(undef, length(source))
    for (index, point) in enumerate(points(source))
        formulation = materialize(point)
        formulation isa AbstractFormulation || throw(ArgumentError(
            "formulation Gridspace produced $(typeof(formulation))",
        ))
        formulations[index] = formulation
    end
    return formulations
end

function _formulations(source::AbstractGrid)
    isempty(source) && return AbstractFormulation[]
    collected = collect(source)
    formulation_type = foldl(
        typejoin,
        (typeof(formulation) for formulation in collected)
    )
    formulation_type <: AbstractFormulation || throw(ArgumentError(
        "formulation Grid must contain completed AbstractFormulation values",
    ))
    formulations = Vector{formulation_type}(undef, length(collected))
    for (index, formulation) in enumerate(collected)
        formulation isa AbstractFormulation || throw(ArgumentError(
            "formulation Grid produced $(typeof(formulation))",
        ))
        formulations[index] = formulation
    end
    return formulations
end

"""
$(TYPEDSIGNATURES)

Evaluate established scalar problem/formulation dispatch for every completed
formulation in a collection. Owners may provide a more specific method to
share immutable lowering work. Every result must have one consistent concrete
type.
"""
function compute(
        problem::AbstractProblemDefinition,
        formulations::AbstractVector{<:AbstractFormulation};
        options::Union{NamedTuple,ComputationOptions} = ComputationOptions()
)
    options = options isa NamedTuple ? ComputationOptions(options) : options
    isempty(formulations) && throw(ArgumentError(
        "a formulation collection must contain at least one formulation",
    ))
    first_result = compute(problem, first(formulations); options)
    check_core_result(typeof(first_result))
    values = Vector{typeof(first_result)}(undef, length(formulations))
    values[1] = first_result
    for index in 2:length(formulations)
        value = compute(problem, formulations[index]; options)
        typeof(value) === eltype(values) || throw(ArgumentError(
            "batched computation produced inconsistent core result types",
        ))
        values[index] = value
    end
    return values
end

"""
$(TYPEDSIGNATURES)

Resolve the formulation source once, materialize every problem point once,
and evaluate their Cartesian product. Optional details are aligned with the
same problem-index-fastest storage order.

# Returns

- A named tuple containing `values`, `(problems, formulations)` axes, and
  `ComputationDetails`, empty or containing `(points=records,)`.
"""
function traverse(problem::ParametricProblem, formulation)
    return with_scan_progress() do receiver
        _traverse(problem,formulation,receiver)
    end
end

function _traverse(problem,formulation,receiver)
    source_id = Grammar.gridpoint_id().source_id
    point_count = length(problem.space)
    point_count > 0 || throw(ArgumentError(
        "higher-order problem space must contain at least one core problem",
    ))

    formulations = _formulations(formulation.inner)
    formulation_count = length(formulations)
    formulation_count > 0 || throw(ArgumentError(
        "higher-order formulation space must contain at least one formulation",
    ))

    point_source = points(problem.space)
    first_item = iterate(point_source)
    first_item === nothing && throw(DimensionMismatch(
        "problem-space iteration ended before its declared cardinality",
    ))
    first_point, state = first_item
    first_problem = materialize(first_point)
    receiver === nothing || report_progress(receiver,
        (kind=:scan, stage=:computing, completed=0, total=point_count*formulation_count,
            point=1, batch=formulation_count))
    first_batch = [Engine.retain_gridpoint(value,
        Grammar.gridpoint_id(;source_id,problem_index=1,formulation_index=index))
        for (index,value) in enumerate(compute(first_problem, formulations; options = problem.options))]
    length(first_batch) == formulation_count || throw(DimensionMismatch(
        "batched computation did not return one result per formulation",
    ))
    first_result = first(first_batch)
    check_core_result(typeof(first_result))
    values = Vector{typeof(first_result)}(
        undef,
        point_count * formulation_count
    )
    @inbounds for formulation_index in 1:formulation_count
        value = first_batch[formulation_index]
        typeof(value) === eltype(values) || throw(ArgumentError(
            "batched computation produced inconsistent core result types",
        ))
        values[1 + (formulation_index - 1) * point_count] = value
    end

    retained = if formulation.options.data.retain_details
        first_record = computation_details(
            typeof(first(formulations)),
            first_result
        )
        records = Vector{typeof(first_record)}(
            undef,
            point_count * formulation_count
        )
        @inbounds for formulation_index in 1:formulation_count
            record = computation_details(
                typeof(formulations[formulation_index]),
                first_batch[formulation_index]
            )
            typeof(record) === eltype(records) || throw(ArgumentError(
                "batched computation produced inconsistent details record types",
            ))
            records[1 + (formulation_index - 1) * point_count] = record
        end
        records
    else
        nothing
    end

    for index in 2:point_count
        receiver === nothing || report_progress(receiver,
            (kind=:scan, stage=:computing, completed=(index-1)*formulation_count, total=point_count*formulation_count,
                point=index, batch=formulation_count))
        item = iterate(point_source, state)
        item === nothing && throw(DimensionMismatch(
            "problem-space iteration ended before its declared cardinality",
        ))
        point, state = item
        resolved_problem = materialize(point)
        batch = [Engine.retain_gridpoint(value,
            Grammar.gridpoint_id(;source_id,problem_index=index,formulation_index=fi))
            for (fi,value) in enumerate(compute(resolved_problem, formulations; options = problem.options))]
        length(batch) == formulation_count || throw(DimensionMismatch(
            "batched computation did not return one result per formulation",
        ))
        @inbounds for formulation_index in 1:formulation_count
            core_result = batch[formulation_index]
            typeof(core_result) === eltype(values) || throw(ArgumentError(
                "higher-order computation produced inconsistent core result types",
            ))
            result_index = index + (formulation_index - 1) * point_count
            values[result_index] = core_result

            if retained !== nothing
                record = computation_details(
                    typeof(formulations[formulation_index]),
                    core_result
                )
                typeof(record) === eltype(retained) || throw(ArgumentError(
                    "higher-order computation produced inconsistent details record types",
                ))
                retained[result_index] = record
            end
        end
    end
    iterate(point_source, state) === nothing || throw(DimensionMismatch(
        "problem-space iteration exceeded its declared cardinality",
    ))
    receiver === nothing || report_progress(receiver,
        (kind=:scan, stage=:computed, completed=point_count*formulation_count, total=point_count*formulation_count,
            batch=formulation_count))

    retained_details = ComputationDetails(retained === nothing ? (;) : (points = retained,))
    axes = (
        problems = problem.space,
        formulations = formulations
    )
    return (; values, details = retained_details, axes)
end

"""
$(TYPEDSIGNATURES)

Compute one target-bearing scalar grid point. Core workflows may add a more
specific lowering route; the general compatibility path materializes exactly
that selected point and never the surrounding finite space.
"""
function compute(
        point::Gridpoint{Target},
        formulation;
        options::Union{NamedTuple,ComputationOptions} = ComputationOptions()
) where {Target <: AbstractProblemDefinition}
    options = options isa NamedTuple ? ComputationOptions(options) : options
    problem = materialize(point)::Target
    return compute(problem, formulation; options)
end

function compute(problem::ParametricProblem, formulation::Combinatorial)
    traversed = traverse(problem, formulation)
    return ParametricResult(
        formulation,
        traversed.values,
        traversed.axes,
        traversed.details
    )
end

"""
$(TYPEDSIGNATURES)

Evaluate a deterministic formulation `Gridspace` with default
[`Combinatorial`](@ref) settings. A scalar problem forms a singleton problem
axis; a problem `Gridspace` forms a Cartesian product with the formulations.

# Keywords

- `options=(;)`: Options supplied to each core computation.

# Returns

- A [`ParametricResult`](@ref) indexed by problem and formulation, with no
  retained traversal details. Select `Combinatorial` explicitly to customize
  its settings.
"""
function compute(
        problem::Union{AbstractProblemDefinition, Gridspace{<:AbstractProblemDefinition}},
        formulations::Gridspace{<:AbstractFormulation};
        options::Union{NamedTuple,ComputationOptions} = ComputationOptions()
)
    options = options isa NamedTuple ? ComputationOptions(options) : options
    return compute(ParametricProblem(problem, options), Combinatorial(formulations))
end

"""
$(TYPEDSIGNATURES)

Evaluate a deterministic formulation `Gridspace` with default
[`Combinatorial`](@ref) settings, retaining the existing `ParametricProblem`
and its core computation options. Return a [`ParametricResult`](@ref).
"""
function compute(
        problem::ParametricProblem,
        formulations::Gridspace{<:AbstractFormulation}
)
    return compute(problem, Combinatorial(formulations))
end
