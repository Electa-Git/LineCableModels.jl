using AWSS3
using LineCableModelsPlayground

const Publisher = LineCableModelsPlayground
const HTTP = Publisher.Bonito.HTTP
const PAYLOAD = collect(codeunits(repeat("distributed-artifact-payload\n", 8_192)))

length(ARGS) == 1 || error("expected the artifact digest")
digest = only(ARGS)
gateway = Publisher.S3ArtifactGateway(
    ENV["LCM_TEST_S3_ENDPOINT"],
    ENV["LCM_TEST_S3_BUCKET"],
    ENV["LCM_TEST_S3_PUBLISHER_ACCESS_KEY"],
    ENV["LCM_TEST_S3_PUBLISHER_SECRET_KEY"];
    prefix="integration",
    allow_insecure=startswith(ENV["LCM_TEST_S3_ENDPOINT"], "http://")
)

function request(method, range="")
    headers = isempty(range) ? Pair{String,String}[] : ["Range" => range]
    return (request=HTTP.Request(method, "/artifacts/sha256/$digest", headers),)
end

full = Publisher.Bonito.HTTPServer.apply_handler(gateway, request("GET"))
full.status == 200 || error("publisher full download returned $(full.status)")
collect(full.body) == PAYLOAD || error("publisher full download changed the payload")

partial = Publisher.Bonito.HTTPServer.apply_handler(
    gateway,
    request("GET", "bytes=1024-4095")
)
partial.status == 206 || error("publisher range download returned $(partial.status)")
collect(partial.body) == PAYLOAD[1025:4096] || error("publisher range was incorrect")

head = Publisher.Bonito.HTTPServer.apply_handler(gateway, request("HEAD"))
head.status == 200 || error("publisher HEAD returned $(head.status)")
isempty(head.body) || error("publisher HEAD unexpectedly returned a body")

backend = gateway.backend
try
    AWSS3.s3_put(
        backend.config,
        backend.bucket,
        "integration/forbidden-write",
        UInt8[0x01],
        "application/octet-stream";
        parse_response=false
    )
    error("read-only publisher identity unexpectedly wrote an artifact")
catch exception
    exception isa Publisher.AWS.AWSException || rethrow()
end

println("S3 artifact round trip and role isolation passed")
