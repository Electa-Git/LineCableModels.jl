@testset "private Julia terminal component" begin
    LCM = LineCableModelsPlayground
    client = RuntimeClient(uuid4())
    terminal = JuliaTerminal(client, :terminal)
    @test terminal.role == "terminal"
    @test terminal.rows == 18
    @test terminal.title == "Julia REPL"
    @test terminal.client === client
    @test JuliaTerminal(client, :other; rows=6).rows == 6
    @test JuliaTerminal(client, :other; rows=60).rows == 60
    @test_throws ArgumentError JuliaTerminal(client, "../terminal")
    @test_throws ArgumentError JuliaTerminal(client, :terminal; rows=5)
    @test_throws ArgumentError JuliaTerminal(client, :terminal; rows=61)
    @test_throws ArgumentError JuliaTerminal(client, :terminal; rows=true)
    @test_throws ArgumentError JuliaTerminal(client, :terminal; title="")
    @test_throws ArgumentError JuliaTerminal(client, :terminal; title=repeat("x",121))
    @test_throws ArgumentError JuliaTerminal(client, :terminal; title="Terminal\n")
    @test fieldnames(JuliaTerminal) == (:client, :role, :title, :rows)
    inspection = ComponentXRay.inspection(terminal)
    @test isempty(inspection.bindings)
    @test length(inspection.actions) == 7
    @test inspection.source.file == "src/widgets/JuliaTerminal.jl"
    @test ".lc-runtime-terminal" in inspection.css_scopes
    @test !occursin(string(client.run_id), repr(inspection))
    session = Bonito.Session()
    try
        @test Bonito.jsrender(session, terminal) !== nothing
    finally
        close(session)
    end
    @test "/widgets/julia-terminal" in first.(LCM.WIDGET_ROUTES)
    css = last(LCM.TERMINAL_STYLES)
    @test !occursin(r"#[0-9a-fA-F]{3,8}\b", css)
    @test !occursin(r"--lc-[a-z-]+\s*:", css)
    @test all(isfile(joinpath(LCM.PLAYGROUND_ROOT,"assets",name)) for name in LCM.TERMINAL_ASSET_NAMES)
end
