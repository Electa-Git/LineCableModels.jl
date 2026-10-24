"""
    LineCableModelsWGLMakieExt

Activate WGLMakie for browser-based LineCableModels rendering.
"""
module LineCableModelsWGLMakieExt

import WGLMakie

activate!() = (WGLMakie.activate!(); :wgl)
make_screen(::AbstractString; kwargs...) = nothing

end # module LineCableModelsWGLMakieExt
