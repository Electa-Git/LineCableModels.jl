using LineCableModelsRuntime, Test
directory = abspath(only(ARGS))
@testset "actual CLI completes owned teardown before exit" begin
    @test isempty(readdir(joinpath(directory, "hosts")))
    store = RuntimeStore(joinpath(directory, "runtime.sqlite"))
    try
        runs = list_runs(store, Principal("browser-fixture"))
        @test !isempty(runs)
        @test all(run -> run.state in (:stopped, :failed), runs)
    finally
        close(store)
    end
end
