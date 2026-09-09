# Private, process-owned Caddy rehearsal of the shipped deployment template.
# No system service, public listener, credential store or trust-store mutation.
import Base64, Sockets

function physical_proxy_port()
    listener = Sockets.listen(Sockets.IPv4("127.0.0.1"),0)
    try Int(last(Sockets.getsockname(listener))) finally close(listener) end
end

function start_physical_proxy(directory, upstream_port, port, proxy_key)
    executable = realpath(ENV["LCM_TEST_CADDY"])
    origin = "https://127.0.0.1:$port"
    root = joinpath(directory,"proxy"); mkpath(root); chmod(root,0o700)
    passwords = Dict(user=>string(uuid4()) * string(uuid4()) for user in ("operator","researcher"))
    hashes = Dict{String,String}()
    for (user,password) in passwords
        # Feed an ephemeral test password through stdin, never process arguments.
        hashes[user] = strip(read(pipeline(`timeout 20s $executable hash-password`;
            stdin=IOBuffer(password * "\n"),stderr=devnull),String))
        startswith(hashes[user],"\$2") || error("Caddy did not produce a bcrypt hash")
    end
    certs = joinpath(directory,"certs")
    template = read(joinpath(@__DIR__,"..","Caddyfile.example"),String)
    # Only ephemeral listener/upstream/certificates differ from the operator template.
    template = replace(template,"{\$LCM_PUBLIC_HOST} {"=>"{\$LCM_PUBLIC_HOST} {\n    bind 127.0.0.1\n    tls \"$(joinpath(certs,"server-cert.pem"))\" \"$(joinpath(certs,"server-key.pem"))\"",
        "reverse_proxy 127.0.0.1:8080"=>"reverse_proxy 127.0.0.1:$upstream_port")
    config = joinpath(root,"Caddyfile")
    write(config,"{\n    admin off\n    auto_https disable_redirects\n}\n" * template); chmod(config,0o600)
    environment = merge(Dict(ENV),Dict("LCM_PUBLIC_HOST"=>origin,"LCM_PROXY_KEY"=>proxy_key,
        "LCM_OPERATOR_PASSWORD_HASH"=>hashes["operator"],"LCM_RESEARCHER_PASSWORD_HASH"=>hashes["researcher"],
        "XDG_DATA_HOME"=>joinpath(root,"data"),"XDG_CONFIG_HOME"=>joinpath(root,"config")))
    runner = CommandRunner(;timeout_seconds=20)
    try
        result = run_owned_command!(runner,setenv(`$executable validate --config $config --adapter caddyfile`,environment))
        result.exitcode == 0 || error("Caddy template validation failed (configuration is private)")
    finally
        close(runner)
    end
    logfile = joinpath(root,"caddy.log")
    output = open(logfile,"w"); chmod(logfile,0o600)
    process = try
        run(pipeline(setenv(`$executable run --config $config --adapter caddyfile`,environment);
            stdin=devnull,stdout=output,stderr=output);wait=false)
    finally
        close(output)
    end
    client = HTTP.Client(transport=HTTP.Transport(tls_config=HTTP.TLS.Config(ca_file=joinpath(certs,"ca.pem")),proxy=nothing))
    headers(user) = ["Authorization"=>"Basic " * Base64.base64encode(user * ":" * passwords[user]),"Origin"=>origin]
    function stop()
        close(client)
        if !process_exited(process)
            kill(process,Base.SIGTERM)
            timedwait(()->process_exited(process),10) == :ok || kill(process,Base.SIGKILL)
        end
        wait(process)
    end
    try
        timedwait(20;pollint=0.1) do
            process_exited(process) && error("Private Caddy exited before readiness")
            try HTTP.get(origin * "/health";client,request_timeout=2,retry=false).status == 200 catch; false end
        end == :ok || error("Private Caddy did not become ready")
    catch
        stop(); rethrow()
    end
    return (;origin,client,headers,stop)
end
