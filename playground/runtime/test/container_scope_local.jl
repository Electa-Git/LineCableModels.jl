using Test
using LineCableModelsRuntime
const ScopeRT = LineCableModelsRuntime
runner = CommandRunner()
try
    host = check_container_host(runner; requested=get(ENV, "LCM_TEST_CONTAINER_RUNTIME", "auto"))
    @testset "actual local container identity without resource allocation" begin
        before = ScopeRT.scoped_container_command(runner, host,
            ["container", "ls", "--all", "--no-trunc", "--format", "{{.ID}}"])
        @test before.exitcode == 0
        scope = container_scope(runner, host)
        @test occursin(r"^[a-f0-9]{64}$", scope)
        @test container_scope(runner, host) == scope
        mktempdir() do directory
            journal = ResourceJournal(joinpath(directory, "receipts"), "scope-audit")
            try
                @test recover_containers!(journal, runner, host) === nothing
                @test isempty(resource_receipts(journal))
            finally
                close(journal)
            end
        end
        after = ScopeRT.scoped_container_command(runner, host,
            ["container", "ls", "--all", "--no-trunc", "--format", "{{.ID}}"])
        @test after.exitcode == 0
        @test sort(split(strip(before.output), '\n')) == sort(split(strip(after.output), '\n'))
        @test isempty(runner.active)
        println("Inspected engine: ", host.engine.name, "; prerequisites: ",
            isempty(host.failures) ? "present (not a launch attestation)" : join(string.(host.failures), ", "))
        println("No containers created, started, stopped or removed.")
    end
finally
    close(runner)
end
