@testitem "PSCAD / station from a TOML file or a named tuple" tags=[:integration, :pscad] begin
    using TOML, PythonCall
    const P = Base.get_extension(LineCableModels, :LineCableModelsPSCADExt)
    # The normalizer of the computation options builds the station of the `remote` option.
    station(remote) =
        computation_options(P.PSCADFormulation, ComputationOptions((; remote))).data.remote
    mktempdir() do directory
        filename = joinpath(directory, "chosen-station.toml")
        fields = Dict{String, Any}(
            "host" => "station", "port" => 54330, "local_root" => "not-created/exchange",
            "shared_root" => raw"Z:\shared space\exchange",
            "pscad_executable" => raw"C:\Program Files (x86)\PSCAD\5.1.0\bin\win64\Pscad.exe")
        save() = open(io -> TOML.print(io, fields), filename, "w")
        save()
        configured = station(filename)
        @test fieldnames(typeof(configured)) ==
              (:host, :ssh, :port, :pscad_executable, :local_root, :shared_root, :timeout)
        @test configured.host == "station"
        @test configured.ssh == ["ssh"]
        @test configured.port == 54330
        @test configured.pscad_executable == fields["pscad_executable"]
        @test configured.local_root == joinpath(directory, "not-created", "exchange")
        @test !ispath(configured.local_root)
        @test configured.shared_root == fields["shared_root"]
        @test configured.timeout == 1800
        relative = cd(dirname(directory)) do
            station(joinpath(basename(directory), basename(filename)))
        end
        @test relative.local_root == configured.local_root
        fields["timeout"] = 60
        fields["ssh"] = ["ts", "ssh"]
        save()
        @test station(filename).timeout == 60
        @test station(filename).ssh == ["ts", "ssh"]
        # The solver identity that the automation server reports.
        identity = Dict("schema" => "2", "version" => "5.1.0", "automation_version" => "3.1.2",
            "common_version" => "3.3.1", "line_constants" => raw"$(HomeDir)\bin\tline.exe")
        @test validate(identity, configured) === identity
        @test_throws ArgumentError validate(Dict("version"=>5), configured)
        for (key, value) in (("schema", "1"), ("version", "5.0.2"), ("line_constants", ""))
            @test_throws ArgumentError validate(merge(identity, Dict(key => value)), configured)
        end
        @test_throws ArgumentError validate(delete!(copy(identity), "line_constants"), configured)
        for (key, value) in (("unknown", 1), ("timeout", true), ("timeout", 1.5),
                ("host", 42), ("host", " "), ("local_root", ""), ("shared_root", ""),
                ("port", "54330"), ("port", true), ("port", 80), ("port", 65536),
                ("pscad_executable", 5), ("pscad_executable", ""), ("ssh", "ssh"),
                ("ssh", String[]), ("ssh", ["ts", 1]), ("remote_root", raw"C:\scratch"),
                ("transport", "ssh"), ("command", ["ts"]))
            saved = copy(fields)
            fields[key] = value
            save()
            @test_throws ArgumentError station(filename)
            empty!(fields); merge!(fields, saved)
        end
        for key in ("port", "pscad_executable")
            saved = copy(fields)
            delete!(fields, key)
            save()
            @test_throws ArgumentError station(filename)
            empty!(fields); merge!(fields, saved)
        end
        @test !ispath(joinpath(directory, "not-created"))
    end
    # A named tuple of the same fields. A relative `local_root` resolves against the current
    # directory.
    fields = (host = "station", port = 54330, shared_root = raw"Z:\exchange",
        pscad_executable = raw"C:\PSCAD\Pscad.exe", local_root = "exchange")
    configured = station(fields)
    @test configured.local_root == abspath("exchange")
    @test (configured.ssh, configured.timeout) == (["ssh"], 1800)
    @test station(merge(fields, (ssh = ["ts", "ssh"], timeout = 60))).ssh == ["ts", "ssh"]
    for (host, port, shared, executable, timeout, ssh) in (
            ("", 54330, "Z:", "Pscad.exe", 1, ["ssh"]), ("station", 1023, "Z:", "Pscad.exe", 1, ["ssh"]),
            ("station", 54330, " ", "Pscad.exe", 1, ["ssh"]), ("station", 54330, "Z:", " ", 1, ["ssh"]),
            ("station", 54330, "Z:", "Pscad.exe", 0, ["ssh"]), ("station", 54330, "Z:", "Pscad.exe", 1, String[]))
        @test_throws ArgumentError station((; host, port, shared_root = shared,
            pscad_executable = executable, local_root = "exchange", timeout, ssh))
    end
    @test_throws ArgumentError station(merge(fields, (local_root = "",)))
    @test_throws ArgumentError station(merge(fields, (remote_root = raw"C:\scratch",)))
    @test_throws ArgumentError station((host = "station", shared_root = raw"Z:\exchange",
        pscad_executable = raw"C:\PSCAD\Pscad.exe", local_root = "exchange"))
    # The refusal of any other value states the two forms.
    for value in (42, ["station.toml"], Dict("host" => "station"))
        @test_throws "path of a station TOML file or a named tuple of its fields" station(value)
    end
end
