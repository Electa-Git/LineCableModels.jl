"""
    LineCableModelsGLMakieExt

Activate GLMakie and create interactive plot screens for LineCableModels.
"""
module LineCableModelsGLMakieExt

import GLMakie

activate!() = (GLMakie.activate!(); :gl)
function make_screen(
        title::AbstractString;
        minimum_size::Tuple{Int, Int} = (1, 1),
        kwargs...
)
    screen = GLMakie.Screen(; title = String(title), kwargs...)
    GLMakie.GLFW.SetWindowSizeLimits(
        screen.glscreen,
        minimum_size...,
        GLMakie.GLFW.DONT_CARE,
        GLMakie.GLFW.DONT_CARE
    )
    return screen
end

end # module LineCableModelsGLMakieExt
