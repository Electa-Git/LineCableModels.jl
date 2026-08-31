# LineCableModels playground

This directory is the clean foundation for the LineCableModels publisher. Its
only page is a static Bonito home. NATS and Quarto are declared as architectural
dependencies but are intentionally not activated by the home page.

## Install the command

Julia 1.12 provides package applications. Install this checkout as a developed
application once:

```sh
julia playground/install.jl
```

Make sure `~/.julia/bin` is on `PATH`. The working-tree source remains live, so
editing the package does not require reinstalling the command.

## Start

```sh
linecablemodels playground start
```

The command starts the server and opens the default browser. Configuration can
come from `HOST`, `PORT`, and `PROXY_URL`, or from explicit options:

```sh
linecablemodels playground start --host 127.0.0.1 --port 8080 --proxy-url .
linecablemodels playground start --no-open
```

Use `Ctrl+C` to stop the server cleanly.
