struct ProbeWorkspace
    input::Any
    invariants::Any
    buffers::Any
end
mutable struct TracedWorkspace
    const input::Int
    const plan::Int
    const buffers::Int
    const trace::Int
end
struct FormulaBuffers end
struct KernelCache end
abstract type AbstractStorage end
struct Table
    cache::Int
end
sample_storage(n) = zeros(n)
function kernel_workspace(n)
    return n
end
trace_buffers(n) = n
initialize_buffers(selected, T, input, plan, buffers) = buffers
function solve(work, n; scratch = nothing)
    storage = zeros(n)
    for cache in 1:n
    end
    return (; work = n, cache_mesh = n)
end
