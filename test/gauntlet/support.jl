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
           report, run_case, run_snapshot, snapshot_path,
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
    using LineCableModels: PhaseDomain, description, domain, observables
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

    function formulation_record(formulation::AnalyticalFormulation)
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

    const _IMMUTABLE_V1_SOURCE_SHA256 = (
        benchmark_132kV_630mm2_flathor_pscad =
        "4227d975e4192b1d0fb35ac3bb98e13b705d27be1fd963a79f7a59da61fee09b",
        benchmark_18kV_1000mm2_trifoil_pscad =
        "a91da88945f95a0a71b7000c6ffd6101c12a4aa5b494a766d59979a9eaee952a",
        benchmark_380kV_2000mm2_flatver_pscad =
        "9cf7f9b9a368e772d640c7f788ef3c55e1bad887e914bb03433d0fc5726b401d",
        benchmark_525kV_1600mm2_bipole_pscad =
        "47ac14e85a75917e964253f422fe228851ca882c08b05c6b13c17f5f149b52fe",
        benchmark_640kV_2000mm2_bipole_pscad =
        "d9d564aa2d6eaf8f717d13838590912be1ed0d3548ec1049f964b9dec6b7cd07",
        benchmark_solid_1000mm2_single_pscad =
        "113871806e03339b82f6aadf4b3acafd63487e8478199e7432a019938386a246",
        benchmark_two_bare_wires_pscad =
        "870cd2ee11011375bbe197df89a91c8d59a8604f7307013d21e0e1a867c3c95b"
    )

    function _reference_source_path(case::GauntletCase)
        return joinpath(
            GAUNTLET_ROOT,
            "reference_sources",
            "gauntlet_pscad_v1_0_0",
            basename(case.source_file) * ".source"
        )
    end

    function _reference_source_digest(case::GauntletCase)
        expected = get(_IMMUTABLE_V1_SOURCE_SHA256, case.name, nothing)
        expected === nothing && return case_digest(case)
        path = _reference_source_path(case)
        isfile(path) || throw(ArgumentError(
            "immutable Gauntlet v1 source is missing: $path",
        ))
        observed = bytes2hex(sha256(read(path)))
        observed == expected || throw(ArgumentError(
            "immutable Gauntlet v1 source SHA-256 does not match: $path",
        ))
        return observed
    end

    _record_value(record::NamedTuple, key::Symbol, default = nothing) = get(record, key, default)
    _record_value(record::AbstractDict, key::Symbol, default = nothing) = get(
        record, key, get(record, String(key), default))

    function _semantic_method_record(record)
        record === nothing && return nothing
        return (description = String(_record_value(record, :description)),)
    end

    function _semantic_analytical_options(record)
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
           occursin("AnalyticalFormulation", type_name)
            return (
                backend = :analytical,
                internal_impedance = _semantic_method_record(_record_value(record, :internal_impedance)),
                insulation_impedance = _semantic_method_record(_record_value(record, :insulation_impedance)),
                earth_impedance = _semantic_method_record(_record_value(record, :earth_impedance)),
                insulation_admittance = _semantic_method_record(_record_value(record, :insulation_admittance)),
                earth_admittance = _semantic_method_record(_record_value(record, :earth_admittance)),
                earth_properties = _semantic_method_record(_record_value(record, :earth_properties)),
                modal_transform = _semantic_method_record(_record_value(record, :modal_transform)),
                equivalent_earth = _semantic_method_record(_record_value(record, :equivalent_earth)),
                options = _semantic_analytical_options(_record_value(record, :options))
            )
        elseif occursin("PSCADFormulation", type_name)
            options = _record_value(record, :options)
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
                options = (output_stem = String(_record_value(options, :output_stem)),)
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
        elseif left isa AbstractArray || left isa Tuple
            size(left) == size(right) || return false
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
        values = observables(parameters)
        size(values.series_impedance) == case.expected_size || throw(DimensionMismatch(
            "$(case.name) expected Z/Y size $(case.expected_size), got " *
            "$(size(values.series_impedance))",
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
            options::NamedTuple = (;)
    )
        prior = record ? load_prior_snapshot(case) : nothing
        root = work_path(case)
        @info "Starting external benchmark" case=case.name mode=record ? :record : :live work_directory=root
        reference = compute(
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
        candidate = compute(case.problem, case.formulation; options)
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
