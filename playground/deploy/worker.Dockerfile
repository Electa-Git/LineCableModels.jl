FROM docker.io/library/julia:1.12-bookworm@sha256:709daad7eccdb0363b080df203601988cff8c6a9581162668b76804cca407888

WORKDIR /workspace

# Resolve and precompile the locked engine stack independently from worker
# orchestration source. The local LineCableModels package is copied only as
# far as dependency precompilation requires.
COPY Project.toml /workspace/Project.toml
COPY src /workspace/src
COPY ext /workspace/ext
COPY playground/worker/Project.toml playground/worker/Manifest.toml /workspace/playground/worker/
COPY playground/protocol /workspace/playground/protocol
COPY playground/vendor /workspace/playground/vendor
RUN julia --startup-file=no --project=playground/worker -e \
    'using Pkg; Pkg.instantiate(; allow_autoprecomp=false); Pkg.precompile(["AWS", "AWSS3", "JSON3", "LineCableModels", "LineCableModelsPlaygroundProtocol", "NATS", "PowerImpedance", "URIs"])'

COPY . .
RUN julia --startup-file=no --project=playground/worker -e 'using LineCableModelsWorker'

ENV PATH="/workspace/playground:${PATH}"
CMD ["lcm", "worker", "start"]
