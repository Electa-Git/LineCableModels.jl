"""$(TYPEDEF)

Base type for external origins of validation evidence.
"""
abstract type Datasource end

"""Identify PSCAD as a datasource of imported validation evidence."""
struct PSCAD <: Datasource end

"""Reserve finite-element solvers as a future validation datasource."""
struct FEM <: Datasource end

"""Resolve a registered datasource symbol through dispatch."""
datasource(name::Symbol) = datasource(Val(name))
datasource(::Val{:pscad}) = PSCAD()
datasource(::Val{:fem}) = FEM()
function datasource(::Val{name}) where {name}
    throw(ArgumentError("unregistered Gauntlet datasource :$name"))
end

"""
    ingest(datasource, source, destination; kwargs...)

Normalize external benchmark evidence into a Gauntlet dataset.
"""
function ingest end

"""
    load(datasource, record)

Reconstruct one native [`Case`](@ref) from a normalized datasource record.
"""
function load end

"""
    decode(datasource, record)

Decode one datasource record into normalized [`Reference`](@ref) evidence.
"""
function decode end

ingest(name::Symbol, args...; kwargs...) = ingest(datasource(name), args...; kwargs...)
load(name::Symbol, args...; kwargs...) = load(datasource(name), args...; kwargs...)
decode(name::Symbol, args...; kwargs...) = decode(datasource(name), args...; kwargs...)

function _fem_not_implemented()
    throw(ErrorException("FEM datasource is not implemented"))
end

ingest(::FEM, args...; kwargs...) = _fem_not_implemented()
load(::FEM, args...; kwargs...) = _fem_not_implemented()
decode(::FEM, args...; kwargs...) = _fem_not_implemented()

"""$(TYPEDEF)

Base type for physical case families.
"""
abstract type Family end
"""$(TYPEDEF)

Identify a coaxial-cable case.
"""
struct Coax <: Family end
"""$(TYPEDEF)

Identify an overhead-line case.
"""
struct Overhead <: Family end
"""$(TYPEDEF)

Identify a mixed overhead/cable case.
"""
struct Mixed <: Family end
"""$(TYPEDEF)

Identify a pipe-type cable case.
"""
struct Pipe <: Family end

"""$(TYPEDEF)

Base type for scientific reconstruction fidelity.
"""
abstract type Fidelity end
"""$(TYPEDEF)

Mark a reconstruction whose geometry and formulation equivalence are verified.
"""
struct Exact <: Fidelity end
"""$(TYPEDEF)

Mark a reconstruction containing one or more declared approximations.
"""
struct Approximate <: Fidelity end
"""$(TYPEDEF)

Mark a case rejected by the source solver before reference output existed.
"""
struct Rejected <: Fidelity end

"""$(TYPEDEF)

Base type for source conductor-reduction states.
"""
abstract type Reduction end
"""$(TYPEDEF)

Mark a source case that retains auxiliary conductors.
"""
struct Retained <: Reduction end
"""$(TYPEDEF)

Mark a source case whose declared auxiliary conductors are reduced.
"""
struct Reduced <: Reduction end
"""$(TYPEDEF)

Mark a source case with no optional conductor reduction.
"""
struct NoReduction <: Reduction end

"""$(TYPEDEF)

Base type for typed scientific and performance checks.
"""
abstract type Check end
"""$(TYPEDEF)

Compare a phase-domain matrix quantity `Q`, such as `:Z` or `:Y`.
"""
struct MatrixCheck{Q} <: Check end
"""$(TYPEDEF)

Compare a sequence or modal quantity `Q`.
"""
struct ModalCheck{Q} <: Check end
"""$(TYPEDEF)

Evaluate and compare an imported vector-fit quantity `Q`.
"""
struct FitCheck{Q} <: Check end
"""$(TYPEDEF)

Evaluate a physical diagnostic `Q`.
"""
struct PhysicalCheck{Q} <: Check end
"""$(TYPEDEF)

Compare bytes and allocations with a pinned performance baseline.
"""
struct PerformanceCheck <: Check end

"""$(TYPEDEF)

Base type for comparison outcomes.
"""
abstract type Verdict end
"""$(TYPEDEF)

Mark a gating comparison that satisfies its tolerance.
"""
struct Pass <: Verdict end
"""$(TYPEDEF)

Mark a gating comparison or execution that failed.
"""
struct Fail <: Verdict end
"""$(TYPEDEF)

Retain a non-gating numerical discrepancy for scientific inspection.
"""
struct Diagnostic <: Verdict end
"""$(TYPEDEF)

Mark a reference channel that the source did not emit.
"""
struct Unavailable <: Verdict end
"""$(TYPEDEF)

Mark a case that its source solver rejected before producing a reference.
"""
struct ReferenceRejected <: Verdict end

"""$(TYPEDEF)

Base type for suite fidelity policies.
"""
abstract type Policy end
"""Gate only cases whose native reconstruction is scientifically exact."""
struct ExactOnly <: Policy end
"""Gate every available comparison, including approximate reconstructions."""
struct AllCases <: Policy end

"""
$(TYPEDEF)

Identify one source matrix terminal and its native phase mapping.

$(TYPEDFIELDS)
"""
struct Port
    "Stable source-terminal identifier."
    id::String
    "One-based cable index."
    cable::Int
    "Cable-component identifier."
    component::String
    "One-based native phase index; zero denotes an eliminated terminal."
    phase::Int
end

"""
$(TYPEDEF)

Describe one deliberate source-to-native approximation.

$(TYPEDFIELDS)
"""
struct Assumption
    "Scientific subject affected by the approximation."
    subject::Symbol
    "Concise explanation of the approximation and its scope."
    detail::String
end

"""
$(TYPEDEF)

Preserve the origin and immutable identity of imported evidence.

$(TYPEDFIELDS)
"""
struct Provenance
    "Datasource identifier."
    datasource::Symbol
    "Datasource version."
    version::String
    "Harvest campaign identifier."
    campaign::String
    "Datasource component definition."
    definition::String
    "SHA-256 digest of the donor source."
    source_sha256::String
    "SHA-256 digest of the effective case specification."
    specification_sha256::String
    "SHA-256 case identifier."
    case_sha256::String
    "Reported source-tool execution time \\[s\\], when available."
    elapsed_seconds::Union{Nothing, Float64}
    "Datasource artifact names and SHA-256 digests."
    artifact_sha256::Dict{String, String}
    "SHA-256 identity shared by reduction variants of one physical case."
    cohort::String
end

"""
$(TYPEDEF)

Store one rational column of a characteristic-admittance fit.

$(TYPEDFIELDS)
"""
struct RationalColumn{T <: Real}
    "Complex poles \\[s⁻¹\\]."
    poles::Vector{Complex{T}}
    "Residues for every output row and pole \\[S/s\\]."
    residues::Matrix{Complex{T}}

    function RationalColumn(
            poles::AbstractVector{<:Complex{T}},
            residues::AbstractMatrix{<:Complex{T}}
    ) where {T <: Real}
        size(residues, 2) == length(poles) || throw(DimensionMismatch(
            "one residue column is required for every pole",
        ))
        return new{T}(collect(poles), Matrix(residues))
    end
end

"""
$(TYPEDEF)

Store one delayed rational group of a propagation-function fit.

$(TYPEDFIELDS)
"""
struct DelayGroup{T <: Real}
    "Propagation delay \\[s\\]."
    delay::T
    "Complex poles \\[s⁻¹\\]."
    poles::Vector{Complex{T}}
    "Residue matrix for every pole \\[s⁻¹\\]."
    residues::Array{Complex{T}, 3}

    function DelayGroup(
            delay::T,
            poles::AbstractVector{<:Complex{T}},
            residues::AbstractArray{<:Complex{T}, 3}
    ) where {T <: Real}
        size(residues, 3) == length(poles) || throw(DimensionMismatch(
            "one residue matrix is required for every propagation pole",
        ))
        return new{T}(delay, collect(poles), Array(residues))
    end
end

"""
$(TYPEDEF)

Represent imported characteristic-admittance and propagation fits.

$(TYPEDFIELDS)
"""
struct Fit{T <: Real}
    "Rational characteristic-admittance columns."
    columns::Vector{RationalColumn{T}}
    "Constant characteristic-admittance matrix \\[S\\]."
    constant::Matrix{T}
    "Delayed rational propagation groups."
    groups::Vector{DelayGroup{T}}
    "Validated evaluation interval \\[Hz\\]."
    frequency_range::Tuple{T, T}

    function Fit(
            columns::AbstractVector{<:RationalColumn{T}},
            constant::AbstractMatrix{T},
            groups::AbstractVector{<:DelayGroup{T}},
            frequency_range::Tuple{T, T}
    ) where {T <: Real}
        n = length(columns)
        size(constant) == (n, n) || throw(DimensionMismatch(
            "the fit constant must be square and match its column count",
        ))
        all(group -> size(group.residues)[1:2] == (n, n), groups) ||
            throw(DimensionMismatch("propagation residues must match the fit size"))
        first(frequency_range) > zero(T) || throw(DomainError(
            frequency_range, "fit frequencies must be positive"
        ))
        last(frequency_range) >= first(frequency_range) || throw(DomainError(
            frequency_range, "fit frequency range must be increasing"
        ))
        return new{T}(
            collect(columns), Matrix(constant), collect(groups), frequency_range
        )
    end
end

"""
$(TYPEDEF)

Hold imported modal transforms, propagation eigenvalues, and transfer channels.

$(TYPEDFIELDS)
"""
struct Modes{T <: Real}
    "Evaluation frequencies \\[Hz\\]."
    frequency::Vector{T}
    "Phase-to-modal transformation matrices \\[dimensionless\\]."
    transform::Union{Nothing, Array{Complex{T}, 3}}
    "Eigenvalues of the per-length ``YZ`` product \\[m⁻²\\]."
    propagation::Union{Nothing, Matrix{Complex{T}}}
    "Calculated characteristic-admittance matrices \\[S\\]."
    characteristic::Union{Nothing, Array{Complex{T}, 3}}
    "Fitted characteristic-admittance matrices \\[S\\]."
    fitted_characteristic::Union{Nothing, Array{Complex{T}, 3}}
    "Calculated phase-domain propagation matrices \\[dimensionless\\]."
    phase_propagation::Union{Nothing, Array{Complex{T}, 3}}
    "Calculated modal propagation channels \\[dimensionless\\]."
    modal_propagation::Union{Nothing, Matrix{Complex{T}}}
    "Fitted phase-domain propagation matrices \\[dimensionless\\]."
    fitted_phase_propagation::Union{Nothing, Array{Complex{T}, 3}}

    function Modes(
            frequency::AbstractVector{T};
            transform = nothing,
            propagation = nothing,
            characteristic = nothing,
            fitted_characteristic = nothing,
            phase_propagation = nothing,
            modal_propagation = nothing,
            fitted_phase_propagation = nothing
    ) where {T <: Real}
        f = collect(frequency)
        all(isfinite, f) || throw(DomainError(f, "modal frequencies must be finite"))
        issorted(f) || throw(ArgumentError("modal frequencies must be sorted"))
        return new{T}(
            f,
            transform,
            propagation,
            characteristic,
            fitted_characteristic,
            phase_propagation,
            modal_propagation,
            fitted_phase_propagation
        )
    end
end

"""
$(TYPEDEF)

Hold imported open- and short-circuit terminal responses.

$(TYPEDFIELDS)
"""
struct Terminal{T <: Real}
    "Evaluation frequencies \\[Hz\\]."
    frequency::Vector{T}
    "Open-circuit terminal response, when emitted."
    open::Union{Nothing, Matrix{Complex{T}}}
    "Short-circuit terminal response, when emitted."
    short::Union{Nothing, Matrix{Complex{T}}}
    "Datasource output families that contained headers but no samples."
    empty::Vector{Symbol}

    function Terminal(
            frequency::AbstractVector{T};
            open = nothing,
            short = nothing,
            empty = Symbol[]
    ) where {T <: Real}
        f = collect(frequency)
        all(value -> isfinite(value) && value > zero(T), f) || throw(DomainError(
            f, "terminal-response frequencies must be finite and positive"
        ))
        issorted(f) || throw(ArgumentError(
            "terminal-response frequencies must be sorted",
        ))
        for values in (open, short)
            values === nothing && continue
            size(values, 1) == length(f) || throw(DimensionMismatch(
                "terminal-response rows must match the frequency count",
            ))
        end
        return new{T}(f, open, short, unique!(collect(empty)))
    end
end

"""
$(TYPEDEF)

Collect phase, sequence, modal, fitted, and terminal reference quantities.

$(TYPEDFIELDS)
"""
struct Reference{P, S, M, F, X, T}
    "Phase-domain line parameters, when available."
    phase::P
    "Sequence-domain line parameters, when available."
    sequence::S
    "Datasource sequence transformation matrix, when available."
    sequence_transform::X
    "Modal and propagation evidence, when available."
    modes::M
    "Imported vector fit, when available."
    fit::F
    "Terminal responses, when available."
    terminal::T
    "Explicit source-to-native terminal ordering."
    ports::Vector{Port}
end

function Reference(; phase = nothing, sequence = nothing, sequence_transform = nothing,
        modes = nothing, fit = nothing, terminal = nothing, ports = Port[])
    Reference(
        phase,
        sequence,
        sequence_transform,
        modes,
        fit,
        terminal,
        collect(ports)
    )
end

function _checks(reference::Reference)
    values = Check[]
    reference.phase === nothing || append!(values, (MatrixCheck{:Z}(), MatrixCheck{:Y}()))
    reference.sequence === nothing || append!(
        values, (ModalCheck{:Z}(), ModalCheck{:Y}())
    )
    if reference.modes !== nothing
        reference.modes.transform === nothing || append!(values, (
            ModalCheck{:transform}(), ModalCheck{:reconstruction}()
        ))
        reference.modes.propagation === nothing || push!(
            values, ModalCheck{:propagation}()
        )
        reference.modes.characteristic === nothing || push!(
            values, ModalCheck{:characteristic}()
        )
        reference.modes.phase_propagation === nothing || push!(
            values, ModalCheck{:phase_propagation}()
        )
        reference.modes.modal_propagation === nothing || push!(
            values, ModalCheck{:modal_propagation}()
        )
    end
    reference.fit === nothing || append!(values, (FitCheck{:Yc}(), FitCheck{:H}()))
    reference.terminal === nothing || push!(values, PhysicalCheck{:terminal}())
    append!(values,
        (
            PhysicalCheck{:symmetry}(),
            PhysicalCheck{:reciprocity}(),
            PhysicalCheck{:positive_real}()
        ))
    return tuple(values...)
end

"""
$(TYPEDEF)

Pair one native model with independently harvested reference evidence.

$(TYPEDFIELDS)
"""
struct Case{P, F, R, G <: Family, D <: Fidelity, U <: Reduction, C <: Tuple}
    "SHA-derived case identifier."
    id::String
    "Physical case family."
    family::G
    "Scientific fidelity of the native reconstruction."
    fidelity::D
    "Datasource conductor-reduction state."
    reduction::U
    "Native `LineParametersProblem`, or `nothing` for a source rejection."
    problem::P
    "Native formulation, or `nothing` for a source rejection."
    formulation::F
    "Imported reference evidence."
    reference::R
    "Checks supported by the available evidence."
    checks::C
    "Deliberate reconstruction assumptions."
    assumptions::Vector{Assumption}
    "Datasource and artifact provenance."
    provenance::Provenance
end

"""
$(TYPEDEF)

Set relative and absolute thresholds for one comparison class.

$(TYPEDFIELDS)
"""
struct Tolerance{T <: Real}
    "Maximum scaled relative error \\[dimensionless\\]."
    rtol::T
    "Absolute error floor in the compared quantity's units."
    atol::T

    function Tolerance(rtol::T, atol::T) where {T <: Real}
        0 <= rtol <= T(1e-3) || throw(DomainError(
            rtol, "gating relative tolerance must lie in [0, 1e-3]"
        ))
        atol >= zero(T) ||
            throw(DomainError(atol, "absolute tolerance must be nonnegative"))
        return new{T}(rtol, atol)
    end
end

Tolerance(; rtol = 1e-6, atol = 0.0) = Tolerance(promote(rtol, atol)...)

"""
$(TYPEDEF)

Select an indexed dataset, checks, fidelity policy, and performance cases.

$(TYPEDFIELDS)
"""
struct Suite{C, K, P <: Policy}
    "Suite identifier."
    name::Symbol
    "Indexed source dataset."
    dataset::C
    "Selected case identifiers."
    ids::Vector{String}
    "Check override; `nothing` uses each case's available checks."
    checks::K
    "Pass/fail policy for exact and approximate cases."
    policy::P
    "Frozen per-check tolerance overrides."
    tolerances::Dict{String, Tolerance{Float64}}
    "Selected case identifiers for performance measurement."
    performance::Vector{String}
end

"""
$(TYPEDEF)

Retain complete numerical diagnostics for one comparison.

$(TYPEDFIELDS)
"""
struct Metrics{T <: Real}
    "Maximum absolute error in the compared quantity's units."
    max_abs::T
    "Maximum scaled relative error \\[dimensionless\\]."
    max_rel::T
    "Frobenius relative error at every frequency \\[dimensionless\\]."
    frobenius::Vector{T}
    "Frequency of the maximum scaled relative error \\[Hz\\]."
    worst_frequency::T
    "Matrix entry of the maximum scaled relative error."
    worst_entry::Tuple{Int, Int}
    "Maximum symmetry residual \\[dimensionless\\]."
    symmetry::T
    "Maximum reciprocity residual \\[dimensionless\\]."
    reciprocity::T
    "Maximum normalized negative Hermitian eigenvalue \\[dimensionless\\]."
    positive_real::T
    "Maximum reconstruction residual \\[dimensionless\\]."
    reconstruction::T
end

"""
$(TYPEDEF)

Associate one typed check with its verdict, tolerance, and diagnostics.

$(TYPEDFIELDS)
"""
struct Comparison{C <: Check, V <: Verdict, M}
    "Executed check."
    check::C
    "Comparison verdict."
    verdict::V
    "Applied numerical tolerance."
    tolerance::Tolerance{Float64}
    "Numerical diagnostics, or `nothing` when unavailable."
    metrics::M
    "Concise scientific explanation."
    detail::String
end

"""
$(TYPEDEF)

Record warmed solver performance on an identified runtime.

$(TYPEDFIELDS)
"""
struct Perf
    "Minimum warmed execution time \\[s\\]."
    minimum_seconds::Float64
    "Median warmed execution time \\[s\\]."
    median_seconds::Float64
    "Allocated bytes per evaluation."
    bytes::Int
    "Allocation count per evaluation."
    allocations::Int
    "Evaluated frequency count."
    frequencies::Int
    "Retained conductor count."
    conductors::Int
    "Median evaluated frequency points per second \\[s⁻¹\\]."
    points_per_second::Float64
    "Julia version."
    julia::String
    "Operating-system kernel."
    os::String
    "Machine architecture."
    arch::String
    "Julia thread count."
    threads::Int
    "BLAS configuration."
    blas::String
    "Minimum-time ratio to the baseline, when available."
    minimum_ratio::Union{Nothing, Float64}
    "Median-time ratio to the baseline, when available."
    median_ratio::Union{Nothing, Float64}
    "Allocated-byte ratio to the baseline, when available."
    bytes_ratio::Union{Nothing, Float64}
    "Allocation-count ratio to the baseline, when available."
    allocations_ratio::Union{Nothing, Float64}
end

"""
$(TYPEDEF)

Store one case execution, its comparisons, and optional performance evidence.

$(TYPEDFIELDS)
"""
struct Trial{C, A, P}
    "Executed case."
    case::C
    "Native result, or `nothing` when execution was unavailable."
    actual::A
    "Complete comparison records."
    comparisons::Vector{Comparison}
    "Warmed performance record, when requested."
    performance::P
end

"""Cache one native modal decomposition during a Gauntlet trial."""
struct ModalState{P, M}
    parameters::P
    modes::M
end

"""
$(TYPEDEF)

Collect Gauntlet trials and their execution interval under one suite.

$(TYPEDFIELDS)
"""
struct Report{S, T}
    "Executed suite."
    suite::S
    "Case trials in declared suite order."
    trials::Vector{T}
    "Suite start time."
    started_at::DateTime
    "Suite completion time."
    finished_at::DateTime
end

"""
    reference(case::Case)

Return the immutable external reference associated with `case`.
"""
reference(case::Case) = case.reference

"""
    assumptions(case::Case)

Return a copy of the declared scientific approximations for `case`.
"""
assumptions(case::Case) = copy(case.assumptions)

"""
    provenance(case::Case)

Return the datasource campaign and artifact provenance associated with `case`.
"""
provenance(case::Case) = case.provenance

function show(io::IO, port::Port)
    print(io, "Port(\"", port.id, "\", cable=", port.cable,
        ", component=\"", port.component, "\", phase=", port.phase, ")")
end

function show(io::IO, case::Case)
    print(io, "Case(\"", case.id, "\", ", nameof(typeof(case.family)), ", ",
        nameof(typeof(case.fidelity)), ", ports=", length(case.reference.ports), ")")
end

function show(io::IO, reference_value::Reference)
    channels = Symbol[]
    reference_value.phase === nothing || push!(channels, :phase)
    reference_value.sequence === nothing || push!(channels, :sequence)
    reference_value.modes === nothing || push!(channels, :modes)
    reference_value.fit === nothing || push!(channels, :fit)
    reference_value.terminal === nothing || push!(channels, :terminal)
    print(io, "Reference(ports=", length(reference_value.ports), ", channels=",
        isempty(channels) ? "none" : join(channels, ", "), ")")
end

function show(io::IO, comparison::Comparison)
    print(io, "Comparison(", _check_name(comparison.check), ", ",
        lowercase(String(nameof(typeof(comparison.verdict)))))
    comparison.metrics === nothing ||
        @printf(io, ", max_rel=%.6g", comparison.metrics.max_rel)
    print(io, ")")
end

function show(io::IO, trial::Trial)
    print(io, "Trial(\"", trial.case.id, "\", ", length(trial.comparisons),
        " comparisons, performance=", trial.performance === nothing ? "no" : "yes", ")")
end

function show(io::IO, performance::Perf)
    @printf(io,
        "Perf(median=%.3f ms, bytes=%d, allocations=%d, frequencies=%d)",
        1e3 * performance.median_seconds,
        performance.bytes,
        performance.allocations,
        performance.frequencies)
end

function show(io::IO, report::Report)
    passed = count(
        trial -> all(
            comparison -> comparison.verdict isa Pass ||
                          comparison.verdict isa Diagnostic ||
                          comparison.verdict isa ReferenceRejected,
            trial.comparisons
        ),
        report.trials)
    print(io, "Report(", report.suite.name, ", ", passed, "/",
        length(report.trials), " acceptable trials)")
end
