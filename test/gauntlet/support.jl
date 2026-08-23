@testmodule GauntletSupport begin
    export ARTIFACT_ROOT, ARTIFACTS_TOML, GAUNTLET_VERSION, ISHEADLESS, WORK_ROOT,
           GauntletCase, artifact_name, backend_archive_name, backend_stage,
           benchmark_local,
           benchmark_metadata, case_digest, cleanup_work, comparison_passes,
           gauntlet_cleanup, gauntlet_force, gauntlet_instrumented, gauntlet_mode,
           finalize_artifacts,
           formulation_record,
           load_prior_snapshot, load_snapshot, performance_comparison, persist_snapshot,
           prepare_artifacts, publish_artifact, release_tag,
           run_case, run_snapshot, snapshot_path,
           validate_case, validate_structure, work_path

    include(joinpath(@__DIR__, "artifacts.jl"))
    using .GauntletArtifacts
    using BenchmarkTools
    using Dates
    using JLD2
    using LinearAlgebra: BLAS
    using Pkg.Artifacts
    using SHA
    using Statistics
    using LineCableModels: PhaseDomain, description, domain
    using LineCableModels.Engine

    const GAUNTLET_ROOT = @__DIR__
    const WORK_ROOT = joinpath(GAUNTLET_ROOT, "cases", ".work")
    const ISHEADLESS = haskey(ENV, "CI") ||
                       get(ENV, "LINECABLEMODELS_HEADLESS", "false") == "true"

    struct GauntletCase{RP, RF, P, F, T}
        name::Symbol
        backend::Symbol
        source_file::String
        reference_problem::RP
        reference_formulation::RF
        problem::P
        formulation::F
        port_order::Vector{String}
        expected_size::NTuple{3, Int}
        tolerances::T

        function GauntletCase(
                name::Symbol,
                backend::Symbol,
                source_file::String,
                reference_problem::RP,
                reference_formulation::RF,
                problem::P,
                formulation::F,
                port_order::Vector{String},
                expected_size::NTuple{3, Int},
                tolerances::T
        ) where {RP, RF, P, F, T}
            case = new{RP, RF, P, F, T}(
                name,
                backend,
                source_file,
                reference_problem,
                reference_formulation,
                problem,
                formulation,
                port_order,
                expected_size,
                tolerances
            )
            return validate_case(case)
        end
    end

    function GauntletCase(
            name::Symbol,
            backend::Symbol,
            source_file::AbstractString,
            reference_problem::RP,
            reference_formulation::RF,
            problem::P,
            formulation::F,
            port_order::AbstractVector{<:AbstractString},
            expected_size::NTuple{3, Int},
            tolerances::T
    ) where {RP, RF, P, F, T}
        return GauntletCase(
            name,
            backend,
            String(source_file),
            reference_problem,
            reference_formulation,
            problem,
            formulation,
            String.(port_order),
            expected_size,
            tolerances
        )
    end

    function _assignments(problem)
        return collect(Iterators.flatten(
            position.conn for position in problem.system.cables
        ))
    end

    function validate_case(case::GauntletCase)
        occursin(r"^[a-z][a-z0-9_]*$", string(case.backend)) ||
            throw(ArgumentError(
                "gauntlet backend must be a lowercase identifier; got $(repr(case.backend))",
            ))
        string(case.name) == case.reference_problem.system.system_id ||
            throw(ArgumentError(
                "gauntlet case name must match the reference system identifier",
            ))
        case.problem.system.system_id == case.reference_problem.system.system_id ||
            throw(ArgumentError(
                "gauntlet reference and candidate system identifiers differ",
            ))
        assignments = _assignments(case.problem)
        reference_assignments = _assignments(case.reference_problem)
        isempty(assignments) && throw(ArgumentError(
            "gauntlet case $(case.name) requires at least one explicit terminal",
        ))
        any(iszero, assignments) && throw(ArgumentError(
            "gauntlet case $(case.name) may not map a terminal to phase 0",
        ))
        all(>(0), assignments) || throw(ArgumentError(
            "gauntlet case $(case.name) phase assignments must be positive",
        ))
        length(unique(assignments)) == length(assignments) || throw(ArgumentError(
            "gauntlet case $(case.name) may not bundle terminals under one phase",
        ))
        sort(assignments) == collect(1:length(assignments)) || throw(ArgumentError(
            "gauntlet case $(case.name) phase assignments must be contiguous from 1",
        ))
        reference_assignments == assignments || throw(ArgumentError(
            "gauntlet case $(case.name) reference and candidate terminal mappings differ",
        ))
        length(case.port_order) == length(assignments) || throw(DimensionMismatch(
            "gauntlet case $(case.name) port order must name every explicit terminal",
        ))
        case.expected_size == (
            length(assignments), length(assignments), length(case.problem.frequencies)
        ) || throw(DimensionMismatch(
            "gauntlet case $(case.name) expected size must match its ports and frequencies",
        ))
        case.formulation.options.kron_reduction && throw(ArgumentError(
            "gauntlet case $(case.name) must disable Kron reduction",
        ))
        case.formulation.options.reduce_bundle && throw(ArgumentError(
            "gauntlet case $(case.name) must disable bundle reduction",
        ))
        case.reference_problem.frequencies == case.problem.frequencies ||
            throw(ArgumentError(
                "gauntlet case $(case.name) reference and candidate frequencies differ",
            ))
        required = (:reference, :regression, :performance)
        all(key -> haskey(case.tolerances, key), required) || throw(ArgumentError(
            "gauntlet tolerances must define reference, regression, and performance",
        ))
        return case
    end

    _method_record(value) = value === nothing ? nothing :
                            (
        type = string(parentmodule(typeof(value)), ".", nameof(typeof(value))),
        description = description(value)
    )

    function formulation_record(formulation::EMTFormulation)
        options = formulation.options
        return (
            type = string(parentmodule(typeof(formulation)), ".", nameof(typeof(formulation))),
            internal_impedance = _method_record(formulation.internal_impedance),
            insulation_impedance = _method_record(formulation.insulation_impedance),
            earth_impedance = _method_record(formulation.earth_impedance),
            insulation_admittance = _method_record(formulation.insulation_admittance),
            earth_admittance = _method_record(formulation.earth_admittance),
            earth_properties = _method_record(formulation.earth_properties),
            modal_transform = _method_record(formulation.modal_transform),
            equivalent_earth = _method_record(formulation.equivalent_earth),
            options = (
                reduce_bundle = options.reduce_bundle,
                kron_reduction = options.kron_reduction,
                ideal_transposition = options.ideal_transposition,
                temperature_correction = options.temperature_correction
            )
        )
    end

    formulation_record(case::GauntletCase) = (
        reference = formulation_record(case.reference_formulation),
        candidate = formulation_record(case.formulation)
    )

    case_digest(case::GauntletCase) = bytes2hex(sha256(read(case.source_file)))

    work_path(case::GauntletCase; root::AbstractString = WORK_ROOT) = joinpath(
        root, string(case.backend), string(case.name))

    function comparison_passes(error::RMSError, tolerance)
        return all((error.absolute .<= tolerance.absolute) .|
                   (error.relative .<= tolerance.relative))
    end

    function validate_structure(
            case::GauntletCase,
            parameters::LineParameters;
            port_order = case.port_order
    )
        domain(parameters) === PhaseDomain || throw(ArgumentError(
            "modal Gauntlet comparison is not implemented; compare canonical phase-domain Z and Y",
        ))
        size(Z(parameters)) == case.expected_size || throw(DimensionMismatch(
            "$(case.name) expected Z/Y size $(case.expected_size), got $(size(Z(parameters)))",
        ))
        port_order == case.port_order || throw(ArgumentError(
            "$(case.name) retained-conductor order does not match the case definition",
        ))
        return parameters
    end

    include("benchmark.jl")
    include("snapshots.jl")

    function benchmark_metadata end

    include("pscad/PSCADBenchmarks.jl")

    function _assert_comparison(comparison, tolerance, label::AbstractString)
        for quantity in (:Z, :Y)
            error = getproperty(comparison, quantity)
            limit = getproperty(tolerance, quantity)
            comparison_passes(error, limit) && continue
            failures = findall(.!((error.absolute .<= limit.absolute) .|
            (error.relative .<= limit.relative)))
            throw(ArgumentError(
                "$label $quantity comparison exceeds tolerance at matrix terms " *
                join(string.(Tuple.(failures)), ", "),
            ))
        end
        return nothing
    end

    function _assert_recordable(case::GauntletCase, regression, performance)
        regression === nothing || _assert_comparison(
            regression,
            case.tolerances.regression,
            "accepted regression"
        )
        performance === nothing || performance.passes === nothing ||
            performance.passes ||
            throw(ArgumentError(
                "local solver performance exceeds the accepted artifact tolerances",
            ))
        return nothing
    end

    function _run_live(
            case::GauntletCase;
            record::Bool = false,
            options::ComputeOptions = ComputeOptions()
    )
        prior = record ? load_prior_snapshot(case) : nothing
        root = work_path(case)
        @info "Starting external benchmark" case=case.name mode=record ? :record : :live work_directory=root
        reference = compute!(
            case.reference_problem,
            case.reference_formulation;
            options
        )
        reference_execution = benchmark_metadata(
            case.reference_problem,
            case.reference_formulation
        )
        reference_execution.backend === case.backend || throw(ArgumentError(
            "reference execution backend $(reference_execution.backend) does not match " *
            "case backend $(case.backend)",
        ))
        @info "Computing LineCableModels result" case = case.name
        candidate = compute!(case.problem, case.formulation; options)
        validate_structure(case, reference)
        validate_structure(case, candidate)
        reference_comparison = compare(reference, candidate)
        @info "Benchmarking LineCableModels calculation" case = case.name
        local_timing = benchmark_local(case; options)
        regression = prior === nothing ? nothing : compare(prior.accepted, candidate)
        performance = prior === nothing ? nothing :
                      performance_comparison(
            prior.metadata["julia_benchmark"],
            local_timing,
            case.tolerances.performance
        )
        record && _assert_recordable(case, regression, performance)
        persisted = nothing
        if record
            @info "Recording benchmark artifact" case = case.name
            persisted = persist_snapshot(
                case,
                reference,
                candidate,
                reference_comparison,
                local_timing;
                reference_execution
            )
        end
        @info "External benchmark completed" case=case.name mode=record ? :record : :live
        return (
            mode = record ? :record : :live,
            reference,
            candidate,
            comparison = reference_comparison,
            regression,
            performance,
            metadata = nothing,
            reference_execution,
            artifact = persisted,
            timings = (
                reference = reference_execution.elapsed_seconds,
                julia = local_timing
            )
        )
    end

    function run_case(case::GauntletCase; options::ComputeOptions = ComputeOptions())
        validate_case(case)
        mode = gauntlet_mode()
        mode === :snapshot && return run_snapshot(case; options)
        if mode === :record
            get(ENV, "LINECABLEMODELS_GAUNTLET_PERSIST", "false") == "true" ||
                throw(ArgumentError(
                    "record mode requires LINECABLEMODELS_GAUNTLET_PERSIST=true",
                ))
            gauntlet_instrumented() && throw(ArgumentError(
                "record mode cannot establish a performance baseline while code " *
                "coverage or allocation tracking is enabled",
            ))
            get(ENV, "LINECABLEMODELS_GAUNTLET_RUNNER", "false") == "true" ||
                throw(ArgumentError(
                    "record mode must run through test/gauntlet/runtests.jl so backend " *
                    "collections are finalized only after every case passes",
                ))
        end
        return _run_live(
            case;
            record = mode === :record,
            options
        )
    end
end
