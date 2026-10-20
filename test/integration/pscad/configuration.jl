@testitem "PSCAD / user-owned TOML and transport arguments" tags=[:integration, :pscad] begin
    using TOML, Base64
    const P = LineCableModels.PSCAD
    mktempdir() do directory
        filename = joinpath(directory, "chosen-station.toml")
        fields = Dict{String, Any}(
            "host" => "user@station", "local_root" => "not-created/exchange",
            "shared_root" => raw"Z:\shared space\exchange",
            "remote_root" => raw"C:\scratch space", "julia_executable" => raw"C:\Julia\julia.exe",
            "python_executable" => raw"C:\Python\python.exe")
        save() = open(io -> TOML.print(io, fields), filename, "w")
        save()
        station = P.RemoteConfig(filename)
        @test station.local_root == joinpath(directory, "not-created", "exchange")
        @test !ispath(station.local_root)
        @test station.shared_root == fields["shared_root"]
        @test station.remote_root == fields["remote_root"]
        @test station.transport === :ssh
        @test station.timeout == 1800
        @test station.pscad_port === nothing
        @test P.validate(nothing, P.RemoteConfig) === nothing
        for port in (1024, 50000, 65535)
            fields["pscad_port"] = port
            save()
            @test P.RemoteConfig(filename).pscad_port == port
            @test P.RemoteConfig("fixture", "shared", "scratch", "julia", "python";
                local_root = directory, pscad_port = port).pscad_port == port
        end
        for port in (true, false, 0, 1023, 65536, -1, 50000.0, "50000")
            fields["pscad_port"] = port
            save()
            @test_throws ArgumentError P.RemoteConfig(filename)
            @test_throws ArgumentError P.RemoteConfig("fixture", "shared", "scratch", "julia", "python";
                local_root = directory, pscad_port = port)
            caught = try
                P.RemoteConfig("fixture", "shared", "scratch", "julia", "python";
                    local_root = directory, pscad_port = port)
            catch error
                error
            end
            @test occursin(repr(port), sprint(showerror, caught))
        end
        delete!(fields, "pscad_port")
        save()
        @test isempty(station.command)
        relative = cd(dirname(directory)) do
            P.RemoteConfig(joinpath(basename(directory), basename(filename)))
        end
        @test relative.local_root == station.local_root
        powershell = "Write-Output 'spaces and quotes'"
        arguments = collect(P.remote_command(station, powershell))
        @test arguments[1:2] == ["ssh", fields["host"]]
        @test occursin(powershell, transcode(String, reinterpret(UInt16, base64decode(last(arguments)))))
        @test first(P.remote_command(Val(:local), station, powershell)) == "powershell.exe"
        fields["transport"] = "local"
        fields["command"] = String[]
        save()
        local_station = P.RemoteConfig(filename)
        @test local_station.transport === :local
        @test collect(P.remote_command(local_station, powershell)) == P._powershell_argv(powershell)
        @test !ispath(local_station.local_root)
        @test_throws ArgumentError validate(Dict("version"=>5), station)
        fields["transport"] = "command"
        fields["command"] = ["ts", "ssh", "{host}", "--direct", "-o", "ServerAliveInterval=15", "--"]
        save()
        station = P.RemoteConfig(filename)
        arguments = collect(P.remote_command(station, powershell))
        @test arguments[1:7] == ["ts", "ssh", fields["host"], "--direct", "-o", "ServerAliveInterval=15", "--"]
        @test arguments[8:end] == P._powershell_argv(powershell)
        fields["command"] = ["wrapper path", "{host}", "literal-{host}", "\$(literal)"]
        save()
        @test collect(P.remote_command(P.RemoteConfig(filename), powershell))[1:4] ==
            ["wrapper path", fields["host"], "literal-{host}", "\$(literal)"]
        for (key, value) in (("unknown", 1), ("timeout", true),
                ("host", 42), ("local_root", ""), ("transport", 1),
                ("command", "ts ssh"), ("command", [1]), ("command", String[]))
            saved = copy(fields)
            fields[key] = value
            save()
            @test_throws ArgumentError P.RemoteConfig(filename)
            empty!(fields); merge!(fields, saved)
        end
        fields["transport"] = "ssh"
        save()
        @test_throws ArgumentError P.RemoteConfig(filename)
        @test !ispath(joinpath(directory, "not-created"))
    end
end

@testitem "PSCAD / station arguments for both modes" tags=[:integration, :pscad] begin
    const P = LineCableModels.PSCAD
    for port in (nothing, 50000)
        config = P.RemoteConfig("fixture", "shared", "station-root", "julia", "python";
            local_root = @__DIR__, pscad_port = port)
        for (formulation, expected_port, mode) in ((Formulation(:pscad), 0, "constants"),
                (P.PSCADTransientFormulation(), something(port, 0), "transient"))
            command = P._supervisor_command(config, "shared-case", "scratch-case", "generated",
                formulation, [1.0, 1e6]; output_stem = "fixture")
            @test occursin("-PSCADPort $expected_port", command)
            @test occursin("-RemoteRoot 'station-root'", command)
            @test occursin("-Mode '$mode'", command)
        end
    end
    cancel = P._cancel_command("shared-toolkit", "scratch-case")
    @test occursin("shared-toolkit\\lifecycle.ps1", cancel)
    @test occursin("Stop-PSCADCase -LocalCase 'scratch-case'", cancel)
    @test all(haskey(P.PSCAD_REMOTE_SOURCES, name) for name in ("probe.py", "lifecycle.ps1"))
end

@testitem "PSCAD / supervisor error decoding" tags=[:integration, :pscad] begin
    const P = LineCableModels.PSCAD
    directory = mkpath(joinpath(@__DIR__, "..", "..", "..", ".linecablemodels", "test", "supervisor-errors"))
    mktempdir(directory) do root
        path = joinpath(root, "transport-stderr.txt")
        @test P._supervisor_errors(path) == ""
        write(path, "#< CLIXML\n<Objs><S S=\"Error\">Primary &lt;failure&gt;_x000D__x000A_</S>" *
            "<S S=\"Error\">Cleanup failure: &amp; details_x000D__x000A_</S></Objs>")
        @test P._supervisor_errors(path) == "Primary <failure>\nCleanup failure: & details"
        write(path, "#< CLIXML\n<Objs><S S=\"Progress\">ignored</S></Objs>")
        @test P._supervisor_errors(path) == ""
        write(path, "plain supervisor error\n")
        @test P._supervisor_errors(path) == "plain supervisor error"
    end
end
