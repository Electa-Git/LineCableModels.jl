struct TracedWorkspace
    input::Int
    plan::Int
    buffers::Int
    trace::Int
end
struct PlainWorkspace
    input::Int
    plan::Int
    buffers::Int
end
initialize_buffers(selected, T, input, plan, buffers) = merge(buffers, (; kernel = plan))
function solve(buffers, n; cache_mesh = nothing)
    local quadrature = buffers.cache
    configure(; scratch = n, cache = :cache)
    return (; plan = n, storage_size = quadrature)
end
