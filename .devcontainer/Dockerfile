FROM julia:1.12.7-bookworm AS julia_runtime

FROM mcr.microsoft.com/devcontainers/base:1-bookworm

COPY --from=julia_runtime /usr/local/julia /usr/local/julia

ENV JULIA_PATH=/usr/local/julia
ENV PATH=/usr/local/julia/bin:${PATH}

RUN julia --version
