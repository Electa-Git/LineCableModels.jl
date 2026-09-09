using LineCableModelsRuntime
import TOML, JSON3
const RT = LineCableModelsRuntime
directory = abspath(only(ARGS))
ids = ("worker-a", "worker-b")
stanzas = [
    broker_user_config(CoordinatorIdentity(); password_environment="LCM_TEST_COORDINATOR_PASSWORD", worker_ids=ids),
    broker_user_config(WorkerIdentity("worker-a"); password_environment="LCM_TEST_A_PASSWORD"),
    broker_user_config(WorkerIdentity("worker-b"); password_environment="LCM_TEST_B_PASSWORD"),
]
base = read(joinpath(@__DIR__, "..", "..", "deploy", "remote", "nats-tls.conf"), String)
config = "max_payload: 262144\n" * replace(base, "users: [" => "users: [\n" * join(stanzas, "\n") * "\n"; count=1)
write(joinpath(directory, "nats.conf"), config)
for (id, password) in (("coordinator", "coordinator-fixture-password"),
        ("worker-a", "worker-a-fixture-password"), ("worker-b", "worker-b-fixture-password"))
    path = joinpath(directory, id * ".password")
    write(path, password)
    chmod(path, 0o600)
end

for identity in ("coordinator","worker-a","worker-b")
    path=joinpath(directory,"artifact-$identity.toml")
    open(path,"w") do io
        TOML.print(io,Dict("access_key_id"=>"artifact-$identity","secret_access_key"=>"fixture-artifact-secret-$identity"))
    end
    chmod(path,0o600)
    prefix=identity=="coordinator" ? "runtime-v1/*" : "runtime-v1/workers/$identity/*"
    actions=identity=="coordinator" ? ["s3:GetObject"] : ["s3:PutObject","s3:DeleteObject"]
    statements=Any[Dict("Effect"=>"Allow","Action"=>actions,"Resource"=>["arn:aws:s3:::lcm-runtime-private/$prefix"])]
    identity=="coordinator" && push!(statements,Dict("Effect"=>"Allow","Action"=>["s3:ListBucket"],
        "Resource"=>["arn:aws:s3:::lcm-runtime-private"]))
    write(joinpath(directory,"artifact-$identity-policy.json"),JSON3.write(Dict("Version"=>"2012-10-17","Statement"=>statements)))
end
