struct BenchmarkCalculation{P, F, O <: NamedTuple}
    id::Symbol
    problem::P
    formulation::F
    options::O

    function BenchmarkCalculation(
            id::Symbol,
            problem::P,
            formulation::F,
            options::O
    ) where {P, F, O <: NamedTuple}
        occursin(r"^[a-z][a-z0-9_]*$", string(id)) || throw(ArgumentError(
            "calculation identifiers must be lowercase; got $(repr(id))",
        ))
        return new{P, F, O}(id, problem, formulation, options)
    end
end

struct BenchmarkDefinition{M, R <: BenchmarkCalculation,
    C <: BenchmarkCalculation, P <: NamedTuple, T}
    id::Symbol
    case_id::Symbol
    collection::Symbol
    source_file::String
    source_sha256::String
    model::M
    reference::R
    candidate::C
    comparison_settings::P
    tolerances::T

    function BenchmarkDefinition(
            id::Symbol,
            case_id::Symbol,
            collection::Symbol,
            source_file::String,
            source_sha256::String,
            model::M,
            reference::R,
            candidate::C,
            comparison_settings::P,
            tolerances::T
    ) where {M, R <: BenchmarkCalculation,
            C <: BenchmarkCalculation, P <: NamedTuple, T}
        occursin(r"^[a-z][a-z0-9_]*$", string(id)) || throw(ArgumentError(
            "benchmark identifiers must be lowercase; got $(repr(id))",
        ))
        occursin(r"^[a-z][a-z0-9_]*$", string(collection)) || throw(ArgumentError(
            "benchmark collections must be lowercase identifiers",
        ))
        case_id === model.id || throw(ArgumentError(
            "benchmark :$id names case :$case_id but loaded :$(model.id)",
        ))
        reference.id == candidate.id && throw(ArgumentError(
            "benchmark calculations must have distinct identifiers",
        ))
        isfile(source_file) || throw(ArgumentError(
            "benchmark source file is missing: $source_file",
        ))
        comparison_settings = BenchmarkTableDefinition(; comparison_settings...).settings
        return new{M, R, C, typeof(comparison_settings), T}(
            id,
            case_id,
            collection,
            source_file,
            source_sha256,
            model,
            reference,
            candidate,
            comparison_settings,
            tolerances
        )
    end
end

function benchmark_definition(
        id::Symbol,
        case_id::Symbol,
        collection::Symbol,
        source_file::AbstractString,
        model,
        reference::BenchmarkCalculation,
        candidate::BenchmarkCalculation,
        comparison_settings,
        tolerances
)
    path = realpath(source_file)
    return BenchmarkDefinition(
        id,
        case_id,
        collection,
        path,
        bytes2hex(sha256(read(path))),
        model,
        reference,
        candidate,
        comparison_settings,
        tolerances
    )
end

function BenchmarkCalculation(id::Symbol, problem, formulation; options::NamedTuple = (;))
    BenchmarkCalculation(id, problem, formulation, options)
end

"""
    benchmark_definition(id::Symbol; kwargs...)

Materialize a catalogue declaration from `gauntlet/benchmarks`, passing explicit
configuration and variations to its constructor. This action does not compute.
"""
function benchmark_definition(id::Symbol; kwargs...)
    paths=String[]
    for (directory, _, names) in walkdir(joinpath(GAUNTLET_ROOT, "benchmarks"))
        string(id)*".jl" in names && push!(paths, joinpath(directory, string(id)*".jl"))
    end
    length(paths) == 1 || throw(ArgumentError("unknown or ambiguous benchmark :$id"))
    declaration=Base.include(@__MODULE__, only(paths))
    return Base.invokelatest(declaration; kwargs...)
end
