"""
$(TYPEDEF)

Supertype for policies that project a completed result space into a finite
space of complete downstream problems.
"""
abstract type AbstractProjectionDefinition end

"""
$(SIGNATURES)

Validate whether a projection definition admits a result space.

Required methods reject unsupported definition/source pairs before later
projection stages run.
"""
function entitle end

"""
$(SIGNATURES)

Select and align the stored products required by a projection definition.
"""
function select end

"""
$(SIGNATURES)

Derive a finite iterable of representative states from selected result
products.
"""
function derive end

@required AbstractProjectionDefinition begin
    entitle(::AbstractProjectionDefinition, source)
    select(::AbstractProjectionDefinition, source)
    derive(::AbstractProjectionDefinition, selected)
    materialize(::AbstractProjectionDefinition, representative)
end

"""
$(TYPEDSIGNATURES)

Project a completed result space into a `Gridspace` of complete downstream
problems.

The fixed stage order is `entitle`, `select`, `derive`, `materialize`, and
`finish`. Projection definitions provide methods for the first four stages;
ParametricBuilder derives `finish` from the materialized problem vector.

# Arguments

- `definition`: Projection policy.
- `source`: Completed result space.

# Returns

- A `Gridspace{P}` preserving the order of the derived representatives, where
  `P` is one concrete downstream problem type.

# Errors

- `ArgumentError`: No representative is derived, a materialized value is a
  result-space envelope, or materialized problem types differ.
"""
@orchestrator AbstractProjectionDefinition function project(
        definition::D,
        source::S
) where {
        D <: AbstractProjectionDefinition,
        S <: AbstractResultSpace
}
    entitled = entitle(definition, source)
    selected = select(definition, entitled)
    representatives = derive(definition, selected)

    item = iterate(representatives)
    item === nothing && throw(ArgumentError(
        "projection must derive at least one representative",
    ))

    representative, state = item
    first_problem = materialize(definition, representative)
    problem_type = typeof(first_problem)
    isconcretetype(problem_type) || throw(ArgumentError(
        "projected problem type must be concrete; got $problem_type",
    ))
    problem_type <: AbstractResultSpace && throw(ArgumentError(
        "projection must materialize complete downstream problems, not result spaces",
    ))

    problems = problem_type[first_problem]
    while true
        item = iterate(representatives, state)
        item === nothing && break
        representative, state = item
        problem = materialize(definition, representative)
        typeof(problem) === problem_type || throw(ArgumentError(
            "projection produced inconsistent downstream problem types",
        ))
        push!(problems, problem)
    end

    return finish(definition, problems)
end

"""
$(TYPEDSIGNATURES)

Wrap ordered, materialized downstream problems in a finite `Gridspace`.
"""
function finish(
        ::AbstractProjectionDefinition,
        problems::Vector{P}
) where {P}
    return Gridspace{P}(identity, (Grid(problems),))
end
