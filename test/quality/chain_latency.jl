# The first call of a space and the first call of the modal analysis of its results, the
# explicit two-step `compute(Gridspace{ModalAnalysisProblem}(result), modal)`, compile in
# seconds. The run of a space is called through `Base.invokelatest`, and the compiler
# compiles it from its values. A change that lets inference walk the run again with partly
# known types makes these calls take minutes or hours. The item runs the `chain` sequence
# of `test/tools/latency.jl` in a fresh process and kills it at 300 s. It fails when the
# two first calls together take more than 30 s or when the process was killed.
@testitem "Quality / chain latency / a space and the modal analysis of its results within 30 s" tags=[:quality, :slow] begin
    tool = joinpath(@__DIR__, "..", "tools", "latency.jl")
    command = `$(Base.julia_cmd()) --project=$(dirname(Base.active_project())) --startup-file=no $tool chain`
    output = IOBuffer()
    process = run(pipeline(ignorestatus(command); stdout = output, stderr = output); wait = false)
    killed = timedwait(() -> process_exited(process), 300.0; pollint = 0.5) === :timed_out
    killed && kill(process)
    wait(process)
    text = String(take!(output))
    times = [parse(Float64, found.captures[1])
             for found in eachmatch(r"(?m)^(?:space|two-step)[^:\n]*first call[^:\n]*: ([0-9.]+) s", text)]
    println("first calls of the space and of the two-step: ", join(string.(times, " s"), " and "),
        killed ? "; killed at 300 s" : "")
    passed = !killed && success(process) && length(times) == 2 && sum(times) <= 30
    passed || println(text)
    @test passed
end
