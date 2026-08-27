"""
$(SIGNATURES)

Build a detached plot recipe through `entitle`, `parse`, `resolve`, `fetch`,
and `finish`, in that order.
"""
@orchestrator AbstractPlotDefinition function make_render(
        ::Type{D},
        source;
        kwargs...
) where {D <: AbstractPlotDefinition}
    entitled = entitle(D, source)
    input = parse(D, entitled; kwargs...)
    resolved = resolve(D, entitled, input)
    pages = fetch(D, entitled, resolved)
    return finish(D, pages)
end
