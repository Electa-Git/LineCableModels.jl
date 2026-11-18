"""
    Playground

Render the Quarto project of this folder and serve its site: the landing page at `/` and each presentation at
`/presentations/<name>/`. Run it with `julia --project=playground -m Playground`, with the options `--host`
(default `127.0.0.1`), `--port` (default `8080`) and `--no-render`, which serves the site of the last render.
"""
module Playground

using HTTP: HTTP
using quarto_jll: quarto

function (@main)(arguments)
    host = "127.0.0.1"
    port = 8080
    render = true
    position = 1
    while position <= length(arguments)
        option = arguments[position]
        if option == "--host"
            host = arguments[position + 1]
            position += 2
        elseif option == "--port"
            port = parse(Int, arguments[position + 1])
            position += 2
        elseif option == "--no-render"
            render = false
            position += 1
        else
            throw(ArgumentError("unknown option $option; the options are --host, --port and --no-render"))
        end
    end
    root = pkgdir(@__MODULE__)
    render && run(`$(quarto()) render $root -M engine:markdown`)
    HTTP.serve(HTTP.fileserver(joinpath(root, "_site")), host, port)
    return 0
end

end
