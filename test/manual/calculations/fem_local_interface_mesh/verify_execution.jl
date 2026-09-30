# Public implementation verification after the completed scientific qualification.
isdefined(@__MODULE__,:ROOT) || include("qualify.jl")
using Test, DelimitedFiles
const OUTPUT = mkpath(joinpath(ROOT,"feature-both-media-execution"))
wrapper = joinpath(ROOT,"getdp-live")
write(wrapper,raw"""#!/usr/bin/env bash
set -o pipefail
if [[ "$1" == "-info" ]]; then
    exec "$LCM_INTERFACE_SOLVER" "$@"
fi
stdbuf -oL -eL "$LCM_INTERFACE_SOLVER" "$@" 2>&1 | tee -a "$LCM_INTERFACE_LOG"
""")
chmod(wrapper,0o755)

function verify_execution(fraction, frequency=.1; reference="localized-terminal-shift")
    dir = mkpath(joinpath(OUTPUT,"f$frequency-gamma-$fraction"))
    marker = joinpath(dir,"complete.toml")
    isfile(marker) && (say("REUSE PUBLIC VERIFICATION ",fraction); return)
    problem, form = fixture(:mixed,frequency,fraction)
    options = (;MESH...,keep_run_directory=true,trace=true,timing=true,
        mesh_policy=:reuse,resume_run_directory=:latest,frequency_workers=1,solver_threads=1,
        getdp_executable=wrapper,gmsh_verbosity=0,getdp_verbosity=4)
    say("PUBLIC MANAGED BEGIN Gamma fraction=",fraction)
    measured = withenv("LCM_INTERFACE_SOLVER"=>GETDP,"LCM_INTERFACE_LOG"=>joinpath(ROOT,"live.log")) do
        @timed compute(problem,form;options)
    end
    result = measured.value
    run = details(result).data.fem.run.run_directory
    native = joinpath(dir,"detached")
    prepare_bundle(native,problem,form)
    cp(joinpath(run,"mesh/model.msh"),joinpath(native,"study.msh");force=true)
    solve!(native)
    qualification = joinpath(ROOT,"mixed-f$(frequency)-gamma$fraction",reference)
    largest_change = 0.
    @testset "Managed and detached Gamma=$fraction" begin
        for (quantity,actual) in (("Z",Z(result)[:,:,1]),("Y",Y(result)[:,:,1]))
            detached = native_matrix(native,quantity)
            expected = native_matrix(qualification,quantity)
            for component in (real,imag)
                a,b,c = component.(actual),component.(detached),component.(expected)
                @test all(isapprox.(a,b;rtol=1e-7,atol=0.))
                @test sign.(a)==sign.(c)
                change = maximum(abs.((a-c)./c))
                largest_change = max(largest_change,change)
                @test change <= .02
            end
            writedlm(joinpath(dir,"$quantity.tsv"),actual)
        end
        reused = withenv("LCM_INTERFACE_SOLVER"=>GETDP,"LCM_INTERFACE_LOG"=>joinpath(ROOT,"live.log")) do
            compute(problem,form;options=(;options...,resume_run_directory=run))
        end
        @test details(reused).data.fem.run.reused
        @test Z(reused)==Z(result) && Y(reused)==Y(result)
    end
    record(marker,Dict("managed_run"=>run,"detached_entry"=>joinpath(native,"study.pro"),
        "seconds"=>measured.time,"compile_seconds"=>measured.compile_time,
        "recompile_seconds"=>measured.recompile_time,"gc_seconds"=>measured.gctime,
        "timed_call_reused_solve"=>details(result).data.fem.run.reused,
        "maximum_change_from_qualified_remesh"=>largest_change,
        "comparison_reference"=>qualification,
        "same_mesh_native_parity"=>true,"resume_passed"=>true))
    say("PUBLIC VERIFICATION COMPLETE fraction=",fraction," elapsed=",measured.time,
        " compile=",measured.compile_time," maximum change from qualified mesh=",largest_change)
end

if abspath(PROGRAM_FILE) == @__FILE__
    for fraction in (0.,.99)
        verify_execution(fraction,.1)
    end
    # At this frequency, the original managed and detached geometry routes
    # already differed by 2.0655% in weak X21 before localization. The captured
    # pre-edit managed builder isolates the feature with a matched reference.
    verify_execution(.99,1e6;reference="managed-reference-terminal-shift")
    say("COMPLETE public managed/native execution verification")
end
