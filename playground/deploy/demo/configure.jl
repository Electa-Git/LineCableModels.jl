"""
    DemoSetup

Generate private, single-operator demo configuration and user service files.
This operator script does not install packages, start services, approve workers,
or run scientific code. Existing output files are never replaced.
"""
module DemoSetup

using TOML, JSON3, Random
import LineCableModelsRuntime as RT

"""Fixed demo worker identities; each receives its own broker and storage grants."""
const WORKERS = ("demo-line", "demo-power")
"""Immutable infrastructure images, independent of the selected executor engine."""
const INFRASTRUCTURE = (
    nats="docker.io/library/nats:2.11-alpine@sha256:e4bf19f15fd3218814a4e3c9e0064e1334bd8aa20d5984b9f1a0afd084f8cc00",
    minio="quay.io/minio/minio:RELEASE.2025-09-07T16-13-09Z@sha256:14cea493d9a34af32f524e538b8346cf79f3321eff8e708c1e2960462bd8936e",
    mc="quay.io/minio/mc:RELEASE.2025-08-13T08-35-41Z@sha256:a7fe349ef4bd8521fb8497f55c6042871b2ae640607cf99d9bede5e9bdf11727",
)

"""
    save(path, value) -> String

Write a newly generated private file, refusing existing files and symlinks.
Parent directories are private. Dictionaries are serialized as TOML.
"""
function save(path, value)
    (ispath(path) || islink(path)) && error("Refusing existing output: $path")
    mkpath(dirname(path); mode=0o700)
    open(path, "w") do io
        chmod(path, 0o600)
        value isa AbstractDict ? TOML.print(io, value) : write(io, value)
    end
    return path
end

"""Return a validated absolute operator path suitable for the fixed service files."""
function operator_path(value)
    value isa String && startswith(value, "/") &&
        occursin(r"^/[A-Za-z0-9_./+-]+$", value) || error("Use an absolute path without spaces: $value")
    return normpath(value)
end

"""Construct the three passive container profiles without importing their engines."""
function profiles(images)
    specifications = (("line-parameters", "scientific", ["line.frequency_scan"]),
        ("power-flow", "scientific", ["powerflow.prepare", "impedance.evaluate"]),
        ("julia-terminal", "terminal", String[]))
    return map(collect(specifications)) do (id, kind, operations)
        image = images[id]
        occursin(r"^[a-zA-Z0-9./_-]+@sha256:[a-f0-9]{64}$", image) || error("Use an immutable RepoDigest for $id")
        Dict("id"=>id, "environment"=>image, "fingerprint"=>last(split(image, "@sha256:")),
            "kind"=>kind, "isolation"=>"container", "operations"=>operations,
            "budget"=>(kind == "terminal" ?
                Dict("cpus"=>0.5, "memory_bytes"=>512*1024^2, "pids"=>64, "scratch_bytes"=>8*1024^2) :
                Dict("cpus"=>1.0, "memory_bytes"=>4*1024^3, "pids"=>256,
                    "scratch_bytes"=>128*1024^2, "prepare_seconds"=>600, "job_seconds"=>120)))
    end
end

"""Return TLS broker and HTTPS artifact settings for one private credential directory."""
function connection(root, identity; broker_port=14222, artifact_port=14443)
    Dict("schema_version"=>1,
        "broker"=>Dict("url"=>"tls://127.0.0.1:$broker_port", "server_name"=>"localhost",
            "password_file"=>joinpath(root, "$identity.password"), "ca_file"=>joinpath(root,"ca.pem"),
            "certificate_file"=>joinpath(root,"client-cert.pem"), "key_file"=>joinpath(root,"client-key.pem")),
        "artifacts"=>Dict("backend"=>"s3", "endpoint"=>"https://127.0.0.1:$artifact_port",
            "bucket"=>"lcm-demo-private", "prefix"=>"runtime-v1",
            "credentials_file"=>joinpath(root,"artifact-$identity.toml"), "ca_file"=>joinpath(root,"ca.pem")))
end

"""Render a private user service belonging to the given target; no boot enablement is implied."""
function service(description, command; target, after="", stop=nothing, pre=nothing, post=nothing, environment="", type="exec")
    """
    [Unit]
    Description=$description
    PartOf=$target
    $(isempty(after) ? "" : "Requires=$after\nAfter=$after")
    [Service]
    Type=$type
    $(pre === nothing ? "" : "ExecStartPre=" * RT.systemd_command(String.(pre)))
    ExecStart=$(RT.systemd_command(String.(command)))
    $(stop === nothing ? "" : "ExecStop=" * RT.systemd_command(String.(stop)))
    $(post === nothing ? "" : "ExecStopPost=" * RT.systemd_command(String.(post)))
    $environment
    UMask=0077
    Restart=no
    TimeoutStartSec=180
    TimeoutStopSec=90
    KillMode=control-group
    StandardInput=null
    """
end

"""
    remote(config) -> String

Generate one new private demo instance on the worker computer. Random credentials,
30-day development certificates, persistent data directories, and service files
are confined to `config["root"]`. Return that directory. The two agents use the
runtime's unchanged managed-agent start and post-stop recovery commands.
"""
function remote(config)
    root, source, depot, caddy = operator_path.([config[k] for k in ("root","source","depot","caddy")])
    (ispath(root) || islink(root)) && error("Use a new private state directory: $root")
    engine = config["engine"]
    engine in ("podman", "docker") || error("Choose podman or docker explicitly")
    executable = something(Sys.which(engine), "")
    isempty(executable) && error("Engine is not installed")
    format = engine == "podman" ? "json" : "{{json .}}"
    info = JSON3.read(read(`$executable info --format $format`, String))
    rootless = engine == "podman" ? get(get(get(info,:host,Dict()),:security,Dict()),:rootless,false) :
        any(x -> occursin("rootless",x), get(info,:SecurityOptions,String[]))
    rootless === true || error("This private demo requires a rootless engine")
    for name in ("lcm-demo-nats","lcm-demo-artifacts","lcm-demo-artifact-init")
        success(pipeline(`$executable container inspect $name`; stdout=devnull, stderr=devnull)) &&
            error("Refusing existing container $name")
    end
    installed_profiles = profiles(config["images"])
    for reference in [getfield.(Ref(INFRASTRUCTURE), propertynames(INFRASTRUCTURE))...;
            [p["environment"] for p in installed_profiles]...]
        success(pipeline(`$executable image inspect $reference`; stdout=devnull, stderr=devnull)) ||
            error("Install this image explicitly first: $reference")
    end
    isfile(caddy) && isdir(depot) && isfile(joinpath(source,"playground/lcm")) || error("Missing tools or source")
    mkpath(root; mode=0o700)
    save(joinpath(root,"instance"), "lcm-demo-v1\n")
    certs = joinpath(root,"certs")
    run(`bash $(joinpath(source,"playground/deploy/remote/generate-dev-certs.sh")) $certs`)
    stanzas = String[]
    broker_env = ["NATS_JETSTREAM_KEY=" * bytes2hex(rand(RandomDevice(), UInt8, 32))]
    artifact_env = ["MINIO_ROOT_USER=lcm-demo-admin", "MINIO_ROOT_PASSWORD=" * bytes2hex(rand(RandomDevice(), UInt8, 32))]
    for (identity, variable) in (("coordinator","COORDINATOR"),(WORKERS[1],"LINE"),(WORKERS[2],"POWER"))
        password = bytes2hex(rand(RandomDevice(), UInt8, 32))
        secret = bytes2hex(rand(RandomDevice(), UInt8, 32))
        credentials = joinpath(root, identity)
        save(joinpath(credentials,"$identity.password"), password)
        save(joinpath(credentials,"artifact-$identity.toml"), Dict("access_key_id"=>"artifact-$identity", "secret_access_key"=>secret))
        for (src, dest) in (("ca.pem","ca.pem"),("worker-cert.pem","client-cert.pem"),("worker-key.pem","client-key.pem"))
            save(joinpath(credentials,dest), read(joinpath(certs,src)))
        end
        push!(broker_env, "LCM_$(variable)_PASSWORD=$password")
        push!(artifact_env, "LCM_$(variable)_ARTIFACT=$secret")
        role = identity == "coordinator" ? RT.CoordinatorIdentity() : RT.WorkerIdentity(identity)
        push!(stanzas, RT.broker_user_config(role; password_environment="LCM_$(variable)_PASSWORD", worker_ids=WORKERS))
        prefix = identity == "coordinator" ? "runtime-v1/*" : "runtime-v1/workers/$identity/*"
        actions = identity == "coordinator" ? ["s3:GetObject"] : ["s3:PutObject","s3:DeleteObject"]
        statements = Any[Dict("Effect"=>"Allow","Action"=>actions,"Resource"=>["arn:aws:s3:::lcm-demo-private/$prefix"])]
        identity == "coordinator" && push!(statements, Dict("Effect"=>"Allow","Action"=>["s3:ListBucket"],"Resource"=>["arn:aws:s3:::lcm-demo-private"]))
        save(joinpath(root,"policies","$identity.json"), JSON3.write(Dict("Version"=>"2012-10-17","Statement"=>statements)))
    end
    save(joinpath(root,"broker.env"), join(broker_env,"\n") * "\n")
    save(joinpath(root,"artifacts.env"), join(artifact_env,"\n") * "\n")
    save(joinpath(root,"nats.conf"), """
    server_name: lcm-private-demo
    port: 4222
    max_payload: 262144
    jetstream { store_dir: "/data/jetstream", max_memory_store: 256MB, max_file_store: 2GB, key: \$NATS_JETSTREAM_KEY }
    tls { cert_file: "/certs/server-cert.pem", key_file: "/certs/server-key.pem", ca_file: "/certs/ca.pem", verify: true, timeout: 2 }
    authorization { users: [ $(join(stanzas,"\n")) ] }
    """)
    save(joinpath(root,"Caddyfile"), """
    {
        admin off
        auto_https off
    }
    https://127.0.0.1:14443 {
        bind 127.0.0.1
        tls $certs/artifact-server-cert.pem $certs/artifact-server-key.pem
        reverse_proxy 127.0.0.1:19000
    }
    """)
    units = joinpath(root,"units")
    target = "lcm-demo-infrastructure.target"
    for directory in ("broker-data","artifact-data")
        mkpath(joinpath(root,directory); mode=0o700)
    end
    # Infrastructure runs as namespace root inside the explicitly rootless engine.
    # It has no capabilities, no host devices, and only its own persistent directory.
    common = [executable,"run","--rm","--pull=never","--cap-drop=ALL","--security-opt=no-new-privileges",
        "--label=lcm.demo=lcm-demo-v1"]
    nats = [common; "--name=lcm-demo-nats"; "--publish=127.0.0.1:14222:4222";
        "--env-file=$root/broker.env"; "--volume=$root/nats.conf:/etc/nats.conf:ro";
        "--volume=$certs:/certs:ro"; "--volume=$root/broker-data:/data"; INFRASTRUCTURE.nats; "-c"; "/etc/nats.conf"]
    minio = [common; "--name=lcm-demo-artifacts"; "--publish=127.0.0.1:19000:9000";
        "--env-file=$root/artifacts.env"; "--env=HOME=/tmp"; "--tmpfs=/tmp:rw,noexec,nosuid,size=16m";
        "--volume=$root/artifact-data:/data"; INFRASTRUCTURE.minio; "server"; "/data";
        "--address=:9000"; "--console-address=:9001"; "--certs-dir=/tmp/empty-certs"]
    for (unit, name, command) in (("broker","lcm-demo-nats",nats),("artifacts","lcm-demo-artifacts",minio))
        save(joinpath(units,"lcm-demo-$unit.service"), service("LCM demo $unit",command; target,
            stop=[executable,"stop","--time=20",name],
            post=["/usr/bin/bash",joinpath(@__DIR__,"cleanup-infrastructure"),executable,root,name]))
    end
    save(joinpath(units,"lcm-demo-storage-tls.service"), service("LCM demo private artifact HTTPS",
        [caddy,"run","--config",joinpath(root,"Caddyfile"),"--adapter","caddyfile"];
        target, after="lcm-demo-artifacts.service"))
    initializer = [common; "--name=lcm-demo-artifact-init"; "--network=container:lcm-demo-artifacts";
        "--env=MC_CONFIG_DIR=/tmp/mc"; "--tmpfs=/tmp:rw,noexec,nosuid,size=16m";
        "--env-file=$root/artifacts.env"; "--volume=$root/policies:/policies:ro";
        "--volume=$(joinpath(@__DIR__,"initialize-artifacts")):/initialize:ro";
        "--entrypoint=/bin/sh"; INFRASTRUCTURE.mc; "/initialize"]
    save(joinpath(units,"lcm-demo-artifact-init.service"), service("LCM demo private storage accounts", initializer;
        target, after="lcm-demo-artifacts.service", type="oneshot",
        pre=["/usr/bin/curl","--fail","--silent","--retry","30","--retry-connrefused","--retry-delay","1",
            "--max-time","2","--retry-max-time","60","--output","/dev/null","http://127.0.0.1:19000/minio/health/ready"],
        post=["/usr/bin/bash",joinpath(@__DIR__,"cleanup-infrastructure"),executable,root,"lcm-demo-artifact-init"]))
    save(joinpath(units,target), "[Unit]\nDescription=LCM private demo infrastructure\nRequires=lcm-demo-broker.service lcm-demo-storage-tls.service lcm-demo-artifact-init.service\nAfter=lcm-demo-broker.service lcm-demo-storage-tls.service lcm-demo-artifact-init.service\n")
    for (id, subset, capacity) in ((WORKERS[1],installed_profiles[[1,3]],4),(WORKERS[2],installed_profiles[[2]],2))
        data = connection(joinpath(root,id),id)
        data["profiles"] = subset
        data["agent"] = Dict("worker_id"=>id, "scratch_root"=>joinpath(root,id,"owned"),"capacity"=>capacity,"container_runtime"=>engine)
        file = save(joinpath(root,id,"agent.toml"), data)
        unit = RT.agent_service_unit(file)
        unit = replace(unit, "[Unit]\n"=>"[Unit]\nPartOf=lcm-demo-workers.target\nAfter=$target\nRequires=$target\n",
            "[Service]\n"=>"[Service]\nEnvironment=JULIA_DEPOT_PATH=$depot\nEnvironment=JULIA_LOAD_PATH=@:@stdlib\nEnvironment=OPENBLAS_NUM_THREADS=1\n")
        save(joinpath(units,RT.agent_unit_name(id)),unit)
    end
    save(joinpath(units,"lcm-demo-workers.target"), "[Unit]\nDescription=LCM demo workers\nRequires=lcm-agent-demo-line.service lcm-agent-demo-power.service\nAfter=lcm-agent-demo-line.service lcm-agent-demo-power.service\n")
    save(joinpath(root,"coordinator","images.toml"), config["images"])
    return root
end

"""
    client(config) -> String

Generate a loopback gateway on port 8081 and a private SSH forwarding service.
Only the copied `coordinator` bundle is required. The existing publisher on 8080
is not referenced. Return the state directory; do not start or enable services.
"""
function client(config)
    root, source, julia, depot, adapter = operator_path.([config[k] for k in ("root","source","julia","depot","ssh_adapter")])
    host = config["host"]
    occursin(r"^[a-zA-Z0-9][a-zA-Z0-9._-]*$",host) || error("Invalid SSH host alias")
    link = joinpath(root,"coordinator")
    data = connection(link,"coordinator"; broker_port=24222, artifact_port=24443)
    data["profiles"] = profiles(TOML.parsefile(joinpath(link,"images.toml")))
    data["workers"] = [Dict("id"=>WORKERS[1],"credential_ref"=>WORKERS[1],"profiles"=>["line-parameters","julia-terminal"],"capacity"=>4),
        Dict("id"=>WORKERS[2],"credential_ref"=>WORKERS[2],"profiles"=>["power-flow"],"capacity"=>2)]
    save(joinpath(root,"control.toml"),data)
    gateway = Dict("schema_version"=>1,"enabled"=>true,
        "gateway"=>Dict("listen_host"=>"127.0.0.1","port"=>8081,"public_origin"=>"http://127.0.0.1:8081"),
        "identity"=>Dict("mode"=>"local-development","principal"=>"developer","administrator"=>true),
        "storage"=>Dict("database"=>joinpath(root,"runtime.sqlite"),"scratch_root"=>joinpath(root,"runs")),
        "publisher"=>Dict("site_directory"=>joinpath(source,"playground/_site")),
        "limits"=>Dict("max_runs"=>8,"max_runs_per_owner"=>2,"startup_seconds"=>120,"shutdown_seconds"=>5,"disconnect_grace_seconds"=>60),
        "control"=>Dict("config_file"=>joinpath(root,"control.toml")))
    file = save(joinpath(root,"gateway.toml"),gateway)
    RT.read_config(file)
    units = joinpath(root,"units")
    tunnel = [adapter,"ssh","-N","-T","-o","BatchMode=yes","-o","ExitOnForwardFailure=yes",
        "-o","ServerAliveInterval=15","-o","ServerAliveCountMax=3",
        "-L","127.0.0.1:24222:127.0.0.1:14222","-L","127.0.0.1:24443:127.0.0.1:14443",host]
    # The tunnel outlives gateway shutdown and remote executor cleanup.
    save(joinpath(units,"lcm-demo-tunnel.service"),service("LCM demo private Kubuntu connection",tunnel; target=""))
    gateway_unit = service("LCM demo on localhost 8081",
        [joinpath(source,"playground/lcm"),"runtime","start","--config",file];
        target="lcm-demo-local.target", after="lcm-demo-tunnel.service",
        environment="Environment=LCM_JULIA=$julia\nEnvironment=JULIA_DEPOT_PATH=$depot\nEnvironment=JULIA_LOAD_PATH=@:@stdlib\nEnvironment=OPENBLAS_NUM_THREADS=1")
    # Signal the launcher first; its FIFO requests cooperative runtime shutdown.
    # The whole group is killed only if that bounded shutdown fails.
    save(joinpath(units,"lcm-demo-gateway.service"), replace(gateway_unit,"KillMode=control-group"=>"KillMode=mixed"))
    save(joinpath(units,"lcm-demo-local.target"), "[Unit]\nDescription=LCM private demo browser gateway\nRequires=lcm-demo-gateway.service\nAfter=lcm-demo-gateway.service\n")
    save(joinpath(root,"operator.env"), "LCM_DEMO_SSH=$adapter\nLCM_DEMO_HOST=$host\nLCM_DEMO_JULIA=$julia\nLCM_DEMO_DEPOT=$depot\n")
    return root
end

"""Generate files for the explicitly selected computer, without starting the demo."""
function main(arguments)
    length(arguments)==2 || error("Usage: configure.jl remote|client /absolute/operator.toml")
    action, path = arguments
    config = TOML.parsefile(path)
    root = action == "remote" ? remote(config) : action == "client" ? client(config) : error("Unknown setup side")
    println("Generated private demo files in $root; install the listed user units next.")
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    DemoSetup.main(ARGS)
end
