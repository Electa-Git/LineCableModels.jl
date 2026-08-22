@testmodule GauntletSupport begin
    export ARTIFACT_ROOT, ISHEADLESS, WORK_ROOT,
           GauntletCase, case_digest, gauntlet_mode,
           load_snapshot, run_case, run_snapshot, snapshot_path,
           validate_structure, work_path, write_snapshot

    using Dates
    using JLD2
    using SHA
    using LineCableModels.Engine

    const GAUNTLET_ROOT = @__DIR__
    const ARTIFACT_ROOT = joinpath(GAUNTLET_ROOT, ".artifacts")
    const WORK_ROOT = joinpath(GAUNTLET_ROOT, "cases", ".work")
    const ISHEADLESS = haskey(ENV, "CI") ||
                       get(ENV, "LINECABLEMODELS_HEADLESS", "false") == "true"
    const SNAPSHOT_FORMAT_VERSION = 1
    const VALID_MODES = (:snapshot, :live, :record)

    struct GauntletCase{P, F, T}
        name::Symbol
        source_file::String
        problem::P
        formulation::F
        pscad_formulation::Symbol
        port_order::Vector{String}
        kron_reduced::Bool
        expected_size::NTuple{3, Int}
        tolerances::T
        compare_legacy_exporter::Bool
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

    snapshot_path(case::GauntletCase; root::AbstractString = ARTIFACT_ROOT) = joinpath(
        root, string(case.name) * ".jld2")

    work_path(case::GauntletCase; root::AbstractString = WORK_ROOT) = joinpath(root, string(case.name))

    function validate_structure(
            case::GauntletCase,
            parameters::LineParameters;
            port_order = case.port_order,
            kron_reduced = case.kron_reduced
    )
        size(Z(parameters)) == case.expected_size || throw(DimensionMismatch(
            "$(case.name) expected Z/Y size $(case.expected_size), got $(size(Z(parameters)))",
        ))
        port_order == case.port_order || throw(ArgumentError(
            "$(case.name) retained-conductor order does not match the case definition",
        ))
        kron_reduced == case.kron_reduced || throw(ArgumentError(
            "$(case.name) Kron-reduction state does not match the case definition",
        ))
        return parameters
    end

    include("snapshots.jl")

    function _load_live_support!()
        if !isdefined(@__MODULE__, :PSCADHarness)
            Base.include(@__MODULE__, joinpath(GAUNTLET_ROOT, "pscad", "PSCADHarness.jl"))
        end
        return Base.invokelatest(Base.getglobal, @__MODULE__, :PSCADHarness)
    end

    function run_case(case::GauntletCase)
        mode = gauntlet_mode()
        mode === :snapshot && return run_snapshot(case)
        harness = _load_live_support!()
        run_live = Base.invokelatest(Base.getglobal, harness, :run_live)
        return Base.invokelatest(run_live, case; record = mode === :record)
    end
end
