using Dates
using Downloads
using JSON3
using UUIDs
using LineCableModelsPlayground
using LineCableModelsPlaygroundProtocol

function await(predicate, timeout_seconds::Real, description::AbstractString)
    timedwait(predicate, timeout_seconds; pollint=0.05) == :ok ||
        error("Timed out waiting for $description")
end

function await_terminal(handle; timeout_seconds=60)
    await(
        () -> handle.state[] in LineCableModelsPlayground.JOB_TERMINAL_STATES,
        timeout_seconds,
        "terminal job state (job=$(handle.job_id[]), state=$(handle.state[]))"
    )
    return handle
end

function submit_echo(client, message)
    request = new_job_request(
        "system.echo",
        Dict("message" => message);
        session_id="remote-stack-smoke-$(uuid4())",
        timeout=Second(90)
    )
    handle = JobHandle()
    submit!(client, handle, request)
    return await_terminal(handle)
end

function submit_cable_constants(client)
    request = new_job_request(
        "cable.constants",
        Dict("frequency_hz" => 50.0);
        session_id="remote-stack-science-$(uuid4())",
        timeout=Second(120)
    )
    handle = JobHandle()
    submit!(client, handle, request)
    return await_terminal(handle; timeout_seconds=120)
end

function request(method::AbstractString, url::AbstractString; headers=Pair{String,String}[])
    output = IOBuffer()
    response = Downloads.request(url; method, headers, output)
    return (; status=response.status, headers=response.headers, body=take!(output))
end

function response_header(response, name::AbstractString)
    index = findfirst(pair -> lowercase(first(pair)) == lowercase(name), response.headers)
    isnothing(index) && error("response did not contain the $name header")
    return last(response.headers[index])
end

broker_url = get(ENV, "NATS_TEST_PUBLISHER_URL", "")
isempty(broker_url) && error("NATS_TEST_PUBLISHER_URL is required")
publisher_url = rstrip(get(ENV, "LCM_SMOKE_PUBLISHER_URL", ""), '/')
isempty(publisher_url) && error("LCM_SMOKE_PUBLISHER_URL is required")
artifact_backend = get(ENV, "LCM_SMOKE_ARTIFACT_BACKEND", "s3")

client = BrokerClient(url=broker_url, heartbeat_ttl=10)
try
    await(() -> client.connection_status[] == :online, 30, "broker connection")
    await(
        () -> any(
            operations -> "system.echo" in operations,
            values(client.available_capabilities[])
        ),
        90,
        "container worker capability heartbeat"
    )

    scientific = submit_cable_constants(client)
    scientific.state[] == :ready || error(
        "container scientific job ended in $(scientific.state[])"
    )
    scientific_result = scientific.last_successful_result[]
    isnothing(scientific_result) && error("scientific job returned no result")
    constants = scientific_result.inline_result["constants"]
    isempty(constants) && error("scientific job returned no cable constants")
    constants[1]["resistance_ohm_per_m"] > 0 || error(
        "scientific job returned a nonpositive resistance"
    )

    marker = "remote-container-artifact-$(uuid4())"
    message = repeat("$marker\n", 1_536)
    first = submit_echo(client, message)
    first.state[] == :ready || error("first job ended in $(first.state[])")
    first_result = first.last_successful_result[]
    isnothing(first_result) && error("first job did not return a result")
    isnothing(first_result.inline_result) || error("large result was unexpectedly inline")
    artifact = first_result.artifact
    isnothing(artifact) && error("large result did not produce an artifact")
    artifact.storage_backend == artifact_backend || error(
        "expected $artifact_backend artifact, got $(artifact.storage_backend)"
    )

    artifact_url = publisher_url * artifact.retrieval_reference
    full = request("GET", artifact_url)
    full.status == 200 || error("artifact GET returned $(full.status)")
    length(full.body) == artifact.size || error("artifact GET returned the wrong size")
    document = JSON3.read(String(copy(full.body)))
    document.echo.message == message || error("artifact payload changed in transit")

    head = request("HEAD", artifact_url)
    head.status == 200 || error("artifact HEAD returned $(head.status)")
    parse(Int, response_header(head, "Content-Length")) == artifact.size ||
        error("artifact HEAD returned the wrong content length")

    range = request("GET", artifact_url; headers=["Range" => "bytes=17-80"])
    range.status == 206 || error("artifact range GET returned $(range.status)")
    range.body == full.body[18:81] || error("artifact byte range is incorrect")

    second = submit_echo(client, message)
    second.state[] == :ready || error("cached job ended in $(second.state[])")
    second_result = second.last_successful_result[]
    second_result.cache_status == "hit" || error(
        "repeated request did not hit the result cache"
    )
    second_result.artifact.sha256 == artifact.sha256 || error(
        "cached request returned a different artifact"
    )

    println(
        "Container worker, JetStream, $artifact_backend artifact, gateway, " *
        "and cache smoke passed"
    )
finally
    close!(client)
end
