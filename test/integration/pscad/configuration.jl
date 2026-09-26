@testitem "PSCAD / user-owned TOML and transport arguments" tags=[:integration] begin
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
        @test station.timeout_seconds == 1800
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
        @test_throws ArgumentError P._validate_solver_identity(Dict("version"=>5), station)
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
        for (key, value) in (("unknown", 1), ("timeout_seconds", true),
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
