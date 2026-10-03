"""
    LineCableModelsCairoMakieExt

Activate CairoMakie for non-interactive LineCableModels rendering.
"""
module LineCableModelsCairoMakieExt

import CairoMakie
public CairoMakie

activate!() = (CairoMakie.activate!(); :cairo)
make_screen(::AbstractString; kwargs...) = nothing

end # module LineCableModelsCairoMakieExt
