include(joinpath(@__DIR__, "..", "..", "common", "PublishedAssets.jl"))
using .PublishedAssets: PublishedSite, published_asset

const STATIC_MIME = Dict(".html"=>"text/html; charset=utf-8", ".css"=>"text/css",
    ".js"=>"text/javascript", ".json"=>"application/json", ".svg"=>"image/svg+xml",
    ".png"=>"image/png", ".jpg"=>"image/jpeg", ".jpeg"=>"image/jpeg", ".gif"=>"image/gif",
    ".woff"=>"font/woff", ".woff2"=>"font/woff2", ".ttf"=>"font/ttf", ".pdf"=>"application/pdf")

function serve_published(stream, site::PublishedSite, path::String)
    bytes = published_asset(site, path)
    isnothing(bytes) && return false
    extension = endswith(path, "/") || !occursin('.', basename(path)) ? ".html" : lowercase(splitext(path)[2])
    mime = get(STATIC_MIME, extension, "application/octet-stream")
    gateway_response(stream, 200, bytes; content_type=mime)
    return true
end
