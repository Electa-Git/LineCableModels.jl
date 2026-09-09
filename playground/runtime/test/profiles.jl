@testset "passive approved profile contracts" begin
    fingerprint = repeat("a", 64)
    native = ProfileDefinition("parameters", "/operator/approved/project", fingerprint;
        operations=("line_parameters",), preparation="representative-cable")
    @test native.kind == :scientific
    @test native.isolation == :trusted_process
    @test native.budget.memory_bytes == 1024^3
    terminal = ProfileDefinition("repl", "registry.test/julia@sha256:" * fingerprint, fingerprint;
        kind=:terminal, isolation=:container)
    @test terminal.kind == :terminal
    @test isempty(terminal.operations)
    @test_throws ArgumentError ProfileDefinition("unsafe", "julia:latest", fingerprint;
        kind=:terminal, isolation=:container)
    @test_throws ArgumentError ProfileDefinition("unsafe", "/julia", fingerprint; kind=:terminal)
    @test_throws ArgumentError ProfileDefinition("empty", "/julia", fingerprint)
    @test_throws ArgumentError ProfileDefinition("legacy", "/julia", fingerprint;
        operations=("line_parameters",), protocol_version=1)
    @test_throws ArgumentError ProfileDefinition("hash", "/julia", "wrong"; operations=("line_parameters",))
    @test_throws ArgumentError ResourceBudget(cpus=NaN)
    @test_throws ArgumentError ResourceBudget(memory_bytes=0)
    @test_throws ArgumentError ResourceBudget(pids=true)
    @test_throws ArgumentError ResourceBudget(scratch_bytes=-1)
    @test_throws ArgumentError ResourceBudget(prepare_seconds=Inf)
    @test_throws ArgumentError ResourceBudget(job_seconds=0)
    profiles = ProfileRegistry()
    @test register!(profiles, native) === native
    @test register!(profiles, terminal) === terminal
    @test_throws ArgumentError register!(profiles, native)
    apps = ApplicationRegistry()
    definition = ApplicationDefinition("study", "Study", :workbench, "/study";
        requirements=(RuntimeRequirement("science", ("parameters",)),))
    register!(apps, LocalApplication(definition, dirname(@__DIR__), joinpath(@__DIR__, "ui_child.jl")))
    @test validate_requirements(apps, profiles) === nothing
    @test_throws ArgumentError validate_requirements(apps, ProfileRegistry())
end
