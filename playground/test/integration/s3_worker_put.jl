using AWSS3
using LineCableModelsWorker

const Worker = LineCableModelsWorker
const PAYLOAD = collect(codeunits(repeat("distributed-artifact-payload\n", 8_192)))

store = Worker.S3ArtifactStore(
    ENV["LCM_TEST_S3_ENDPOINT"],
    ENV["LCM_TEST_S3_BUCKET"],
    ENV["LCM_TEST_S3_WORKER_ACCESS_KEY"],
    ENV["LCM_TEST_S3_WORKER_SECRET_KEY"];
    prefix="integration",
    allow_insecure=startswith(ENV["LCM_TEST_S3_ENDPOINT"], "http://"),
    max_inline_bytes=16
)

reference = Worker.store_artifact!(store, PAYLOAD; media_type="application/octet-stream")
reference.storage_backend == "s3" || error("worker did not return an S3 artifact")
reference.size == length(PAYLOAD) || error("worker returned the wrong artifact size")

try
    AWSS3.s3_get(
        store.config,
        store.bucket,
        Worker.artifact_object_key(store, "sha256", reference.sha256);
        raw=true,
        retry=false
    )
    error("write-only worker identity unexpectedly read an artifact")
catch exception
    exception isa LineCableModelsWorker.AWS.AWSException || rethrow()
end

println(reference.sha256)
