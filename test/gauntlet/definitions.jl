struct BenchmarkCalculation{P, F, O <: NamedTuple}
    id::Symbol
    owner::Symbol
    problem::P
    formulation::F
    options::O

    function BenchmarkCalculation(
            id::Symbol,
            owner::Symbol,
            problem::P,
            formulation::F,
            options::O
    ) where {P, F, O <: NamedTuple}
        occursin(r"^[a-z][a-z0-9_]*$", string(id)) || throw(ArgumentError(
            "calculation identifiers must be lowercase; got $(repr(id))",
        ))
        owner in (:engine, :uq, :external) || throw(ArgumentError(
            "calculation owner must be :engine, :uq, or :external",
        ))
        return new{P, F, O}(id, owner, problem, formulation, options)
    end
end

function benchmark_calculation(
        id::Symbol,
        owner::Symbol,
        problem,
        formulation;
        options::NamedTuple = (;)
)
    return BenchmarkCalculation(id, owner, problem, formulation, options)
end

"""
    LineParametersPolicy(; quantities=(:Z, :Y), bands=(:all,),
        normalizations=(:reference_rms,), atol=nothing,
        fundamental=50.0, harmonics=50, unsupported=(;))

Select observables, stored-frequency bands, RMS normalizations and numerical-zero
tolerances for a benchmark. Units and band selection follow `Engine.compare`.
The ordered benchmark operands determine the reference, independently of backend.
"""
struct LineParametersPolicy{Q, B, N, A, U}
    quantities::Q
    bands::B
    normalizations::N
    atol::A
    fundamental::Float64
    harmonics::Int
    unsupported::U
end

function LineParametersPolicy(; quantities = (:Z, :Y), bands = (:all,),
        normalizations = (:reference_rms,), atol = nothing,
        fundamental = 50.0, harmonics = 50, unsupported = (;))
    !isempty(quantities) && all(q -> q in (:Z, :Y, :R, :L, :G, :C), quantities) ||
        throw(ArgumentError("benchmark quantities must select Z, Y, R, L, G, or C"))
    !isempty(bands) || throw(ArgumentError("benchmark needs at least one frequency band"))
    !isempty(normalizations) &&
    all(n -> n in (:reference_rms, :pointwise), normalizations) ||
        throw(ArgumentError("benchmark normalizations must select :reference_rms or :pointwise"))
    allunique(quantities) && allunique(bands) && allunique(normalizations) ||
        throw(ArgumentError("benchmark quantities, bands and normalizations must not contain duplicates"))
    return LineParametersPolicy(Tuple(quantities), Tuple(bands), Tuple(normalizations),
        atol, Float64(fundamental), Int(harmonics), unsupported)
end
struct UQMomentPolicy end

struct OwnedBenchmark{M, R <: BenchmarkCalculation,
    C <: BenchmarkCalculation, P, T}
    id::Symbol
    case_id::Symbol
    collection::Symbol
    source_file::String
    source_sha256::String
    model::M
    reference::R
    candidate::C
    comparison_policy::P
    tolerances::T

    function OwnedBenchmark(
            id::Symbol,
            case_id::Symbol,
            collection::Symbol,
            source_file::String,
            source_sha256::String,
            model::M,
            reference::R,
            candidate::C,
            comparison_policy::P,
            tolerances::T
    ) where {M, R <: BenchmarkCalculation,
            C <: BenchmarkCalculation, P, T}
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
        return new{M, R, C, P, T}(
            id,
            case_id,
            collection,
            source_file,
            source_sha256,
            model,
            reference,
            candidate,
            comparison_policy,
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
        comparison_policy,
        tolerances
)
    path = realpath(source_file)
    return OwnedBenchmark(
        id,
        case_id,
        collection,
        path,
        bytes2hex(sha256(read(path))),
        model,
        reference,
        candidate,
        comparison_policy,
        tolerances
    )
end
