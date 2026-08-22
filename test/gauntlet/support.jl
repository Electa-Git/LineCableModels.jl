@testmodule GauntletSupport begin
    export ARTIFACT_ROOT, ARTIFACTS_TOML, ISHEADLESS, WORK_ROOT,
           GauntletCase, artifact_name, benchmark_local,
           case_digest, comparison_passes, gauntlet_mode,
           load_prior_snapshot, load_snapshot, performance_comparison, persist_snapshot,
           publish_artifact, run_case, run_snapshot, snapshot_path,
           validate_case, validate_structure, work_path

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
    const ARTIFACT_ROOT = joinpath(GAUNTLET_ROOT, ".artifacts")
    const ARTIFACTS_TOML = joinpath(GAUNTLET_ROOT, "Artifacts.toml")
    const WORK_ROOT = joinpath(GAUNTLET_ROOT, "cases", ".work")
    const ISHEADLESS = haskey(ENV, "CI") ||
                       get(ENV, "LINECABLEMODELS_HEADLESS", "false") == "true"
    const SNAPSHOT_FORMAT_VERSION = 2
    const VALID_MODES = (:snapshot, :live, :record)

    struct GauntletCase{P, F, T}
        name::Symbol
        source_file::String
        problem::P
        formulation::F
        pscad_formulation::Symbol
        port_order::Vector{String}
        expected_size::NTuple{3, Int}
        tolerances::T

        function GauntletCase(
                name::Symbol,
                source_file::String,
                problem::P,
                formulation::F,
                pscad_formulation::Symbol,
                port_order::Vector{String},
                expected_size::NTuple{3, Int},
                tolerances::T
        ) where {P, F, T}
            case = new{P, F, T}(
                name,
                source_file,
                problem,
                formulation,
                pscad_formulation,
                port_order,
                expected_size,
                tolerances
            )
            return validate_case(case)
        end
    end

    function GauntletCase(
            name::Symbol,
            source_file::AbstractString,
            problem::P,
            formulation::F,
            pscad_formulation::Symbol,
            port_order::AbstractVector{<:AbstractString},
            expected_size::NTuple{3, Int},
            tolerances::T
    ) where {P, F, T}
        return GauntletCase(
            name,
            String(source_file),
            problem,
            formulation,
            pscad_formulation,
            String.(port_order),
            expected_size,
            tolerances
        )
    end

    function _assignments(case::GauntletCase)
        return collect(Iterators.flatten(
            position.conn for position in case.problem.system.cables
        ))
    end

    function validate_case(case::GauntletCase)
        assignments = _assignments(case)
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

    function formulation_record(case::GauntletCase)
        formulation = case.formulation
        options = formulation.options
        return (
            type = string(parentmodule(typeof(formulation)), ".", nameof(typeof(formulation))),
            pscad = string(case.pscad_formulation),
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

    function gauntlet_mode()
        value = Symbol(get(ENV, "LINECABLEMODELS_GAUNTLET_MODE", "snapshot"))
        value in VALID_MODES || throw(ArgumentError(
            "LINECABLEMODELS_GAUNTLET_MODE must be snapshot, live, or record; got $(repr(value))",
        ))
        haskey(ENV, "CI") && value !== :snapshot &&
            throw(ArgumentError(
                "CI permits only LINECABLEMODELS_GAUNTLET_MODE=snapshot",
            ))
        return value
    end

    case_digest(case::GauntletCase) = bytes2hex(sha256(read(case.source_file)))

    work_path(case::GauntletCase; root::AbstractString = WORK_ROOT) = joinpath(root, string(case.name))

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
            "modal PSCAD validation is not implemented because PSCAD does not provide modal Z and Y",
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

    function _load_live_support!()
        if !isdefined(@__MODULE__, :PSCADHarness)
            Base.include(@__MODULE__, joinpath(GAUNTLET_ROOT, "pscad", "PSCADHarness.jl"))
        end
        return Base.invokelatest(Base.getglobal, @__MODULE__, :PSCADHarness)
    end

    function run_case(case::GauntletCase)
        validate_case(case)
        mode = gauntlet_mode()
        mode === :snapshot && return run_snapshot(case)
        if mode === :record
            get(ENV, "LINECABLEMODELS_GAUNTLET_PERSIST", "false") == "true" ||
                throw(ArgumentError(
                    "record mode requires LINECABLEMODELS_GAUNTLET_PERSIST=true",
                ))
        end
        harness = _load_live_support!()
        run_live = Base.invokelatest(Base.getglobal, harness, :run_live)
        return Base.invokelatest(run_live, case; record = mode === :record)
    end
end
