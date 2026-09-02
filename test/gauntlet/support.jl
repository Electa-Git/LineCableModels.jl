@testmodule GauntletSupport begin
    export ARTIFACT_ROOT, ARTIFACTS_TOML, ISHEADLESS, SNAPSHOT_SCHEMA_VERSION,
           WORK_ROOT,
           AbstractCaseVariation, CaseDefinition, CaseParameter, CompositeVariation,
           ExactOverrides, LoadedCase, NoVariation, ParameterGrids,
           RelativeStandardUncertainty,
           BenchmarkCalculation, LineParametersPolicy, MomentBenchmark, MomentResult,
           OwnedBenchmark, UQMomentPolicy,
           UQ_MONTE_CARLO_TRIALS,
           GauntletCase, artifact_name, bind_published_artifact,
           collection_archive_name, collection_release, collection_stage,
           benchmark_local, benchmark_stage,
           benchmark_calculation, benchmark_definition, calculation_record,
           benchmark_metadata, benchmark_digest, case_digest, cleanup_work,
           comparison_passes,
           case_definition, case_index, case_parameter, compose_variations,
           correlation_record,
           gauntlet_cleanup, gauntlet_instrumented, gauntlet_mode,
           gauntlet_stage_force, finalize_staging,
           formulation_record,
           extract_moments, moment_comparison_passes, moment_error_summary,
           parameter_manifest,
           numerical_input_sha256, implementation_record, repository_provenance,
           load_case, load_prior_snapshot, load_snapshot, performance_comparison,
           reference_case, reference_grid,
           persist_snapshot,
           package_collection, prepare_staging, release_tag,
           report, run_case, run_snapshot, snapshot_path,
           run_benchmark, uq_inner_formulation, uq_moment_tolerances,
           validate_case, validate_structure, variation_record, work_path

    include(joinpath(@__DIR__, "artifacts.jl"))
    using .GauntletArtifacts
    using BenchmarkTools
    using Dates
    using JLD2
    using LinearAlgebra: BLAS
    using Pkg.Artifacts
    using SHA
    using Statistics
    import TOML
    import LineCableModels
    using LineCableModels: AbstractGrid, Grid, Gridspace, PhaseDomain,
                           description, domain, observe, observables
    using LineCableModels.Engine
    import LineCableModels.Grammar: ComputationOptions, computation_options

    const GAUNTLET_ROOT = @__DIR__
    const WORK_ROOT = joinpath(GAUNTLET_ROOT, "benchmarks", ".work")
    const ISHEADLESS = haskey(ENV, "CI") ||
                       get(ENV, "LINECABLEMODELS_HEADLESS", "false") == "true"

    include(joinpath(@__DIR__, "case_loader.jl"))
    include(joinpath(@__DIR__, "provenance.jl"))
    include(joinpath(@__DIR__, "reference_grid.jl"))
    include(joinpath(@__DIR__, "comparisons", "uq_moments.jl"))
    include(joinpath(@__DIR__, "uq_benchmarks.jl"))
    include(joinpath(@__DIR__, "owned.jl"))

    struct GauntletCase{RP, RF, P, F, T, V, C}
        name::Symbol
        case_id::Symbol
        backend::Symbol
        source_file::String
        case_source_file::String
        case_sha256::String
        parameter_manifest::Vector{NamedTuple}
        applied_variation::V
        correlation::C
        reference_problem::RP
        reference_formulation::RF
        problem::P
        formulation::F
        port_order::Vector{String}
        expected_size::NTuple{3, Int}
        tolerances::T

        function GauntletCase(
                name::Symbol,
                case_id::Symbol,
                backend::Symbol,
                source_file::String,
                case_source_file::String,
                case_sha256::String,
                parameter_manifest::Vector{NamedTuple},
                applied_variation::V,
                correlation::C,
                reference_problem::RP,
                reference_formulation::RF,
                problem::P,
                formulation::F,
                port_order::Vector{String},
                expected_size::NTuple{3, Int},
                tolerances::T
        ) where {RP, RF, P, F, T, V, C}
            case = new{RP, RF, P, F, T, V, C}(
                name,
                case_id,
                backend,
                source_file,
                case_source_file,
                case_sha256,
                parameter_manifest,
                applied_variation,
                correlation,
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
            name,
            backend,
            String(source_file),
            String(source_file),
            bytes2hex(sha256(read(source_file))),
            NamedTuple[],
            (kind = :inline_legacy,),
            (rule = :inline_legacy,),
            reference_problem,
            reference_formulation,
            problem,
            formulation,
            String.(port_order),
            expected_size,
            tolerances
        )
    end

    function GauntletCase(
            name::Symbol,
            backend::Symbol,
            source_file::AbstractString,
            model::LoadedCase,
            reference_formulation::RF,
            formulation::F,
            tolerances::T
    ) where {RF, F, T}
        problem = model.nominal_problem
        reference_problem = LineParametersProblem(
            problem.system;
            temperature = problem.temperature,
            earth_props = problem.earth_props,
            frequencies = problem.frequencies
        )
        return GauntletCase(
            name,
            model.id,
            backend,
            String(source_file),
            model.source_file,
            model.source_sha256,
            parameter_manifest(model),
            variation_record(model.variation),
            correlation_record(model),
            reference_problem,
            reference_formulation,
            problem,
            formulation,
            copy(model.port_order),
            model.expected_size,
            tolerances
        )
    end

    function computation_options(
            ::Val{GauntletCase},
            options::NamedTuple
    )::ComputationOptions
        allowed = (:output_basis, :reference, :candidate, :benchmark)
        unknown = filter(key -> key ∉ allowed, keys(options))
        isempty(unknown) || throw(ArgumentError(
            "unknown Gauntlet computation options: $(sort!(collect(unknown)))",
        ))
        normalized = merge(
            (
                output_basis = :pul,
                reference = (;),
                candidate = (;),
                benchmark = (samples = 10, seconds = 10.0)
            ),
            options
        )
        normalized.output_basis in (:pul, :total) || throw(ArgumentError(
            "Gauntlet output_basis must be :pul or :total; got " *
            repr(normalized.output_basis),
        ))
        normalized.reference isa NamedTuple || throw(ArgumentError(
            "Gauntlet reference computation options must be a named tuple",
        ))
        normalized.candidate isa NamedTuple || throw(ArgumentError(
            "Gauntlet candidate computation options must be a named tuple",
        ))
        normalized.benchmark isa NamedTuple || throw(ArgumentError(
            "Gauntlet benchmark options must be a named tuple",
        ))
        benchmark_unknown = filter(
            key -> key ∉ (:samples, :seconds),
            keys(normalized.benchmark)
        )
        isempty(benchmark_unknown) || throw(ArgumentError(
            "unknown Gauntlet benchmark options: " *
            "$(sort!(collect(benchmark_unknown)))",
        ))
        benchmark = merge(
            (samples = 10, seconds = 10.0),
            normalized.benchmark
        )
        benchmark.samples isa Integer && !(benchmark.samples isa Bool) &&
        benchmark.samples > 0 || throw(ArgumentError(
            "Gauntlet benchmark samples must be a positive integer",
        ))
        benchmark.seconds isa Real && isfinite(benchmark.seconds) &&
        benchmark.seconds > 0 || throw(ArgumentError(
            "Gauntlet benchmark seconds must be a positive finite number",
        ))
        return (
            output_basis = normalized.output_basis,
            reference = normalized.reference,
            candidate = normalized.candidate,
            benchmark = (
                samples = Int(benchmark.samples),
                seconds = Float64(benchmark.seconds)
            )
        )
    end

    function _assignments(problem)
        return copy(problem.system.connection_order)
    end

    function validate_case(case::GauntletCase)
        occursin(r"^[a-z][a-z0-9_]*$", string(case.name)) ||
            throw(ArgumentError(
                "gauntlet benchmark ID must be a lowercase identifier; got $(repr(case.name))",
            ))
        occursin(r"^[a-z][a-z0-9_]*$", string(case.backend)) ||
            throw(ArgumentError(
                "gauntlet backend must be a lowercase identifier; got $(repr(case.backend))",
            ))
        isfile(case.source_file) || throw(ArgumentError(
            "gauntlet benchmark source is missing: $(case.source_file)",
        ))
        isfile(case.case_source_file) || throw(ArgumentError(
            "gauntlet case source is missing: $(case.case_source_file)",
        ))
        case.case_sha256 == case_digest(case) || throw(ArgumentError(
            "gauntlet case source changed after it was loaded: $(case.case_source_file)",
        ))
        string(case.case_id) == case.reference_problem.system.system_id ||
            throw(ArgumentError(
                "gauntlet case identifier must match the reference system identifier",
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

    function formulation_record(formulation::LineParametersFormulation)
        options = formulation.options
        methods = formulation.methods
        return (
            type = string(parentmodule(typeof(formulation)), ".", nameof(typeof(formulation))),
            internal_impedance = _method_record(methods.internal_impedance),
            insulation_impedance = _method_record(methods.insulation_impedance),
            earth_impedance = _method_record(methods.earth_impedance),
            insulation_admittance = _method_record(methods.insulation_admittance),
            semicon_admittance = _method_record(methods.semicon_admittance),
            earth_admittance = _method_record(methods.earth_admittance),
            earth_properties = _method_record(methods.earth_properties),
            equivalent_earth = _method_record(methods.equivalent_earth),
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

    case_digest(case::GauntletCase) = bytes2hex(sha256(read(case.case_source_file)))
    benchmark_digest(case::GauntletCase) = bytes2hex(sha256(read(case.source_file)))

    _record_value(record::NamedTuple, key::Symbol, default = nothing) = get(record, key, default)
    _record_value(record::AbstractDict, key::Symbol,
        default = nothing) = get(
        record, key, get(record, String(key), default))

    function _semantic_method_record(record)
        record === nothing && return nothing
        return (description = String(_record_value(record, :description)),)
    end

    function _semantic_line_parameters_options(record)
        output = _record_value(record, :output, :parameters)
        output isa Val && (output = only(typeof(output).parameters))
        return (
            reduce_bundle = _record_value(record, :reduce_bundle),
            kron_reduction = _record_value(record, :kron_reduction),
            ideal_transposition = _record_value(record, :ideal_transposition),
            temperature_correction = _record_value(record, :temperature_correction),
            output
        )
    end

    function _semantic_formulation_record(record)
        type_name = String(_record_value(record, :type, ""))
        if occursin("EMTFormulation", type_name) ||
           occursin("AnalyticalFormulation", type_name) ||
           occursin("LineParametersFormulation", type_name)
            return (
                backend = :line_cable_models,
                internal_impedance = _semantic_method_record(_record_value(record, :internal_impedance)),
                insulation_impedance = _semantic_method_record(_record_value(record, :insulation_impedance)),
                earth_impedance = _semantic_method_record(_record_value(record, :earth_impedance)),
                insulation_admittance = _semantic_method_record(_record_value(record, :insulation_admittance)),
                semicon_admittance = _semantic_method_record(_record_value(record, :semicon_admittance)),
                earth_admittance = _semantic_method_record(_record_value(record, :earth_admittance)),
                earth_properties = _semantic_method_record(_record_value(record, :earth_properties)),
                equivalent_earth = _semantic_method_record(_record_value(record, :equivalent_earth)),
                options = _semantic_line_parameters_options(_record_value(record, :options))
            )
        elseif occursin("PSCADFormulation", type_name)
            return (
                backend = :pscad,
                earth_impedance = (
                    description = String(_record_value(
                        _record_value(record, :earth_impedance), :description
                    )),
                    pscad_field = _record_value(_record_value(record, :earth_impedance), :pscad_field),
                    pscad_value = _record_value(_record_value(record, :earth_impedance), :pscad_value),
                    pscad_readback = _record_value(_record_value(record, :earth_impedance), :pscad_readback)
                ),
                earth_admittance = _semantic_method_record(_record_value(record, :earth_admittance)),
                insulation_admittance = _semantic_method_record(_record_value(record, :insulation_admittance)),
                options = (;)
            )
        end
        throw(ArgumentError("unsupported Gauntlet formulation record $type_name"))
    end

    function _semantic_formulation_pair(record)
        return (
            reference = _semantic_formulation_record(_record_value(record, :reference)),
            candidate = _semantic_formulation_record(_record_value(record, :candidate))
        )
    end

    function _same_problem_structure(left, right)
        typeof(left) === typeof(right) || return false
        if left isa AbstractDict
            keys(left) == keys(right) || return false
            return all(key -> _same_problem_structure(left[key], right[key]), keys(left))
        elseif left isa AbstractArray
            size(left) == size(right) || return false
            return all(_same_problem_structure.(left, right))
        elseif left isa Tuple
            length(left) == length(right) || return false
            return all(_same_problem_structure.(left, right))
        elseif left isa Union{Nothing, Missing, Bool, Number, Symbol, AbstractString, Char}
            return isequal(left, right)
        elseif isstructtype(typeof(left))
            return all(
                name -> _same_problem_structure(getfield(left, name), getfield(right, name)),
                fieldnames(typeof(left))
            )
        end
        return isequal(left, right)
    end

    work_path(case::GauntletCase;
        root::AbstractString = WORK_ROOT) = joinpath(
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
        impedance = observe(parameters, Z)
        size(impedance) == case.expected_size || throw(DimensionMismatch(
            "$(case.name) expected Z/Y size $(case.expected_size), got " *
            "$(size(impedance))",
        ))
        port_order == case.port_order || throw(ArgumentError(
            "$(case.name) retained-conductor order does not match the case definition",
        ))
        return parameters
    end

    include(joinpath(@__DIR__, "benchmark.jl"))
    include(joinpath(@__DIR__, "snapshots.jl"))

    function benchmark_metadata end

    include(joinpath(@__DIR__, "pscad", "PSCADBenchmarks.jl"))

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
            options::NamedTuple = (;)
    )
        run = computation_options(Val(GauntletCase), options)
        reference = run.reference
        if case.reference_formulation isa PSCADBenchmarks.PSCADFormulation &&
           !haskey(reference, :remote)
            reference = merge(
                reference,
                (remote = PSCADBenchmarks._load_config(),)
            )
        end
        reference_options = merge(
            reference,
            (output_basis = run.output_basis,)
        )
        candidate_options = merge(
            run.candidate,
            (output_basis = run.output_basis,)
        )
        prior = record ? load_prior_snapshot(case) : nothing
        root = work_path(case)
        @info "Starting external benchmark" case=case.name mode=record ? :record : :live work_directory=root
        reference = compute(
            case.reference_problem,
            case.reference_formulation;
            options = reference_options
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
        candidate = compute(
            case.problem,
            case.formulation;
            options = candidate_options
        )
        validate_structure(case, reference)
        validate_structure(case, candidate)
        reference_comparison = compare(reference, candidate)
        @info "Benchmarking LineCableModels calculation" case = case.name
        local_timing = benchmark_local(
            case;
            options = candidate_options,
            samples = run.benchmark.samples,
            seconds = run.benchmark.seconds
        )
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

    function run_case(case::GauntletCase; options::NamedTuple = (;))
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
