"""
    LineCableModelsCairoMakieExt

Activate CairoMakie for non-interactive LineCableModels rendering.
"""
module LineCableModelsCairoMakieExt

import CairoMakie

activate!() = (CairoMakie.activate!(); :cairo)
make_screen(::AbstractString; kwargs...) = nothing

function __init__()
    activate!()
    return nothing
end

end # module LineCableModelsCairoMakieExt
