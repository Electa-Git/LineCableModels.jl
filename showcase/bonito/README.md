# Bonito cable-design showcase

This standalone prototype proves that Bonito, reveal.js, WGLMakie, and
LineCableModels can rebuild and render one live `CableDesign`. It is intentionally
independent of Pluto and does not add presentation dependencies to the package.

Instantiate the pinned application environment from the repository root:

```sh
julia --project=showcase/bonito -e 'using Pkg; Pkg.instantiate()'
```

Start the presentation:

```sh
julia --project=showcase/bonito showcase/bonito/serve.jl
```

Open <http://127.0.0.1:8080>. The launcher reads `HOST` (default
`127.0.0.1`), `PORT` (default `8080`), and `PROXY_URL` (default `.`). A relative
proxy URL requires the proxy to forward both HTTP and WebSocket traffic.

Generate the backend-neutral CairoMakie smoke render:

```sh
julia --project=showcase/bonito showcase/bonito/export.jl
```

Run the focused tests:

```sh
julia --project=showcase/bonito showcase/bonito/test/runtests.jl
```

The reveal.js 6.0.1 JavaScript, CSS, theme, and notes assets are pinned remote
assets for this prototype. The first app construction therefore requires
network access; vendoring the same files later only requires changing the
central asset declarations in `app.jl`.
