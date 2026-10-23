@testitem "PSCAD / user-owned TOML configuration" tags=[:integration, :pscad] begin
    using TOML
    const P = LineCableModels.PSCAD
    mktempdir() do directory
        filename = joinpath(directory, "chosen-station.toml")
        fields = Dict{String, Any}(
            "host" => "station", "port" => 54330, "local_root" => "not-created/exchange",
            "shared_root" => raw"Z:\shared space\exchange",
            "pscad_executable" => raw"C:\Program Files (x86)\PSCAD\5.1.0\bin\win64\Pscad.exe")
        save() = open(io -> TOML.print(io, fields), filename, "w")
        save()
        station = P.RemoteConfig(filename)
        @test fieldnames(P.RemoteConfig) ==
              (:host, :ssh, :port, :pscad_executable, :local_root, :shared_root, :timeout)
        @test station.host == "station"
        @test station.ssh == ["ssh"]
        @test station.port == 54330
        @test station.pscad_executable == fields["pscad_executable"]
        @test station.local_root == joinpath(directory, "not-created", "exchange")
        @test !ispath(station.local_root)
        @test station.shared_root == fields["shared_root"]
        @test station.timeout == 1800
        relative = cd(dirname(directory)) do
            P.RemoteConfig(joinpath(basename(directory), basename(filename)))
        end
        @test relative.local_root == station.local_root
        fields["timeout"] = 60
        fields["ssh"] = ["ts", "ssh"]
        save()
        @test P.RemoteConfig(filename).timeout == 60
        @test P.RemoteConfig(filename).ssh == ["ts", "ssh"]
        # The solver identity that the automation server reports.
        identity = Dict("schema" => "2", "version" => "5.1.0", "automation_version" => "3.1.2",
            "common_version" => "3.3.1", "line_constants" => raw"$(HomeDir)\bin\tline.exe")
        @test validate(identity, station) === identity
        @test_throws ArgumentError validate(Dict("version"=>5), station)
        for (key, value) in (("schema", "1"), ("version", "5.0.2"), ("line_constants", ""))
            @test_throws ArgumentError validate(merge(identity, Dict(key => value)), station)
        end
        @test_throws ArgumentError validate(delete!(copy(identity), "line_constants"), station)
        for (key, value) in (("unknown", 1), ("timeout", true), ("timeout", 1.5),
                ("host", 42), ("host", " "), ("local_root", ""), ("shared_root", ""),
                ("port", "54330"), ("port", true), ("port", 80), ("port", 65536),
                ("pscad_executable", 5), ("pscad_executable", ""), ("ssh", "ssh"),
                ("ssh", String[]), ("ssh", ["ts", 1]), ("remote_root", raw"C:\scratch"),
                ("transport", "ssh"), ("command", ["ts"]))
            saved = copy(fields)
            fields[key] = value
            save()
            @test_throws ArgumentError P.RemoteConfig(filename)
            empty!(fields); merge!(fields, saved)
        end
        for key in ("port", "pscad_executable")
            saved = copy(fields)
            delete!(fields, key)
            save()
            @test_throws ArgumentError P.RemoteConfig(filename)
            empty!(fields); merge!(fields, saved)
        end
        @test !ispath(joinpath(directory, "not-created"))
    end
    station = P.RemoteConfig("station", 54330, raw"Z:\exchange", raw"C:\PSCAD\Pscad.exe";
        local_root = "exchange")
    @test station.local_root == abspath("exchange")
    @test (station.ssh, station.timeout) == (["ssh"], 1800)
    for (host, port, shared, executable, timeout, ssh) in (
            ("", 54330, "Z:", "Pscad.exe", 1, ["ssh"]), ("station", 1023, "Z:", "Pscad.exe", 1, ["ssh"]),
            ("station", 54330, " ", "Pscad.exe", 1, ["ssh"]), ("station", 54330, "Z:", " ", 1, ["ssh"]),
            ("station", 54330, "Z:", "Pscad.exe", 0, ["ssh"]), ("station", 54330, "Z:", "Pscad.exe", 1, String[]))
        @test_throws ArgumentError P.RemoteConfig(host, port, shared, executable;
            local_root = "exchange", timeout, ssh)
    end
    @test_throws ArgumentError P.RemoteConfig("station", 54330, "Z:", "Pscad.exe"; local_root = "")
end
