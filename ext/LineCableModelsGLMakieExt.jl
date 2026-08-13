module LineCableModelsGLMakieExt

using GLMakie
using LineCableModels

activate!() = (GLMakie.activate!(); :gl)
function make_screen(title::AbstractString; kwargs...)
    GLMakie.Screen(; title = String(title), kwargs...)
end

function __init__()
    activate!()
    return nothing
end

end # module LineCableModelsGLMakieExt
