using TestItemRunner

include("artifacts.jl")
using .GauntletArtifacts

const GAUNTLET_FILTER = testitem -> :gauntlet in testitem.tags

cleanup = gauntlet_cleanup()
completed = false

try
    force = gauntlet_force()
    mode = gauntlet_mode()
    if mode === :record
        gauntlet_instrumented() && throw(ArgumentError(
            "record mode cannot establish a performance baseline while code " *
            "coverage or allocation tracking is enabled",
        ))
        ENV["LINECABLEMODELS_GAUNTLET_RUNNER"] = "true"
        prepare_artifacts(; force)
    end

    @run_package_tests(filter=GAUNTLET_FILTER, verbose=true)

    if mode === :record
        for collection in finalize_artifacts(; force)
            @info "Recorded Gauntlet backend collection" backend=collection.backend version=collection.gauntlet_version archive=collection.archive
        end
    end
    global completed = true
finally
    if cleanup && completed
        work_root = cleanup_work()
        @info "Cleaned Gauntlet working directory" path=work_root
    elseif cleanup
        @warn "Preserving Gauntlet working directory after failed run" path=GauntletArtifacts.WORK_ROOT
    end
end
