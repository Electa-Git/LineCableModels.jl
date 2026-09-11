"""
    LineCableModelsWGLMakieExt

Activate WGLMakie for browser-based LineCableModels rendering.
"""
module LineCableModelsWGLMakieExt

using LineCableModels
using WGLMakie

activate!() = (WGLMakie.activate!(); :wgl)
make_screen(::AbstractString; kwargs...) = nothing

function __init__()
    activate!()
    return nothing
end

end # module LineCableModelsWGLMakieExt
