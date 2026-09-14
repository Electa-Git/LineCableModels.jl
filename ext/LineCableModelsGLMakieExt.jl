"""
    LineCableModelsGLMakieExt

Activate GLMakie and create interactive plot screens for LineCableModels.
"""
module LineCableModelsGLMakieExt

using GLMakie
using LineCableModels

activate!() = (GLMakie.activate!(); :gl)
function make_screen(
        title::AbstractString;
        minimum_size::Tuple{Int, Int} = (1, 1),
        aspect_size::Tuple{Int, Int} = (4, 3),
        kwargs...
)
    screen = GLMakie.Screen(; title = String(title), kwargs...)
    GLMakie.GLFW.SetWindowSizeLimits(
        screen.glscreen,
        minimum_size...,
        GLMakie.GLFW.DONT_CARE,
        GLMakie.GLFW.DONT_CARE
    )
    divisor = gcd(aspect_size...)
    GLMakie.GLFW.SetWindowAspectRatio(screen.glscreen,
        aspect_size[1] ÷ divisor, aspect_size[2] ÷ divisor)
    return screen
end

function __init__()
    activate!()
    return nothing
end

end # module LineCableModelsGLMakieExt
