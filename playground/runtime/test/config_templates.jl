@testset "operator proxy template is inert and owner-scoped" begin
    mktempdir() do directory
        source = joinpath(@__DIR__, "..", "proxy.example.toml")
        path = joinpath(directory, "runtime.toml")
        cp(source, path)
        @test_throws ArgumentError read_config(path) # Missing private key is fatal.
        secrets = mkdir(joinpath(directory, "secrets"); mode=0o700)
        keyfile = joinpath(secrets, "proxy-key")
        key = repeat("a", 64) # Disposable test value, never an operator credential.
        write(keyfile, key)
        chmod(keyfile, 0o600)
        config = read_config(path)
        @test config.enabled
        @test config.listen_host == "127.0.0.1"
        @test config.identity isa ProxyIdentity
        @test config.control === nothing
        @test !ispath(joinpath(directory, "state"))
        headers = ["X-LCM-Principal"=>"operator", "X-LCM-Proxy-Key"=>key]
        @test authenticate(config.identity, headers, "127.0.0.1").administrator
        @test !authenticate(config.identity,
            ["X-LCM-Principal"=>"researcher", "X-LCM-Proxy-Key"=>key],
            "127.0.0.1").administrator
        @test_throws AccessDenied authenticate(config.identity, headers, "192.0.2.1")
        @test_throws AccessDenied authenticate(config.identity, Pair{String,String}[], "127.0.0.1")
    end
end
