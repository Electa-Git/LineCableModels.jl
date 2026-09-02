FROM docker.io/library/julia:1.12-bookworm@sha256:709daad7eccdb0363b080df203601988cff8c6a9581162668b76804cca407888

RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates curl tar util-linux \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace

# Keep network resolution, the pinned Quarto download, and dependency
# precompilation below the application-source layer. A widget or broker edit
# should not reinstall the publishing stack.
COPY playground/Project.toml playground/Manifest.toml playground/bootstrap.sh /workspace/playground/
COPY playground/protocol /workspace/playground/protocol
COPY playground/vendor /workspace/playground/vendor
RUN ./playground/bootstrap.sh --publisher-only --dependencies-only \
    && julia --startup-file=no --project=playground -e \
        'using Pkg; Pkg.precompile(["AWS", "AWSS3", "Bonito", "BonitoWidgets", "JSON3", "LineCableModelsPlaygroundProtocol", "NATS", "Quarto", "URIs"])'

COPY playground /workspace/playground
RUN ./playground/bootstrap.sh --publisher-only

ENV PATH="/root/.local/bin:${PATH}"
EXPOSE 8080
CMD ["lcm", "playground", "start", "--no-open", "--no-render"]
