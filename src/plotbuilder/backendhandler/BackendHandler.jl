"""
    LineCableModels.PlotBuilder.BackendHandler

Coordinate optional Makie backends without loading them from the core package.
"""
module BackendHandler

export backend_available, current_backend_symbol, ensure_backend!, make_screen,
       next_fignum, renderfig, set_backend!, with_backend

const _BACKEND_EXTENSIONS = Dict(
    :cairo => :LineCableModelsCairoMakieExt,
    :gl => :LineCableModelsGLMakieExt,
    :wgl => :LineCableModelsWGLMakieExt
)

const FIG_NO = Base.Threads.Atomic{Int}(1)

_parent_package() = parentmodule(parentmodule(@__MODULE__))

function _makie_extension()
    return Base.get_extension(_parent_package(), :LineCableModelsMakieExt)
end

function _backend_extension(backend::Symbol)
    name = get(_BACKEND_EXTENSIONS, backend, nothing)
    name === nothing && throw(
        ArgumentError("Unknown backend :$(backend). Use :cairo, :gl, or :wgl."),
    )
    return Base.get_extension(_parent_package(), name)
end

"Return whether `backend` has been explicitly loaded."
backend_available(backend::Symbol) = _backend_extension(backend) !== nothing

"Return the active Makie backend as `:cairo`, `:gl`, `:wgl`, `:unknown`, or `:none`."
function current_backend_symbol()
    ext = _makie_extension()
    return ext === nothing ? :none : Base.invokelatest(ext.current_backend_symbol)
end

"""
    set_backend!(backend; force=false)

Activate an explicitly loaded Makie backend.

`force` is retained for source compatibility and has no effect.
"""
function set_backend!(backend::Symbol; force::Bool = false)
    ext = _backend_extension(backend)
    ext === nothing && throw(
        ArgumentError(
        "Backend :$(backend) is not loaded. Run `using $(_backend_package(backend))` first.",
    ),
    )
    return Base.invokelatest(ext.activate!)
end

function _backend_package(backend::Symbol)
    backend === :cairo && return "CairoMakie"
    backend === :gl && return "GLMakie"
    backend === :wgl && return "WGLMakie"
    throw(ArgumentError("Unknown backend :$(backend). Use :cairo, :gl, or :wgl."))
end

"Ensure that an explicitly loaded backend is active."
function ensure_backend!(backend::Union{Nothing, Symbol} = nothing)
    backend !== nothing && return set_backend!(backend)
    current = current_backend_symbol()
    current in keys(_BACKEND_EXTENSIONS) && return current
    throw(
        ArgumentError(
        "No Makie backend is active. Load CairoMakie, GLMakie, or WGLMakie first.",
    ),
    )
end

"Run `f` with an explicitly loaded backend and restore the previous backend."
function with_backend(f::Function, backend::Symbol; force::Bool = false)
    previous = current_backend_symbol()
    set_backend!(backend; force = force)
    try
        return f()
    finally
        if previous in keys(_BACKEND_EXTENSIONS) && previous != backend
            set_backend!(previous; force = force)
        end
    end
end

"Create a backend-specific screen when supported."
function make_screen(
        title::AbstractString;
        backend::Symbol = current_backend_symbol(),
        kwargs...
)
    ext = _backend_extension(backend)
    return ext === nothing ? nothing :
           Base.invokelatest(ext.make_screen, String(title); kwargs...)
end

function make_screen(backend::Symbol, title::AbstractString; kwargs...)
    make_screen(title; backend, kwargs...)
end

"Display a Makie figure through the loaded plotting extension."
function renderfig(fig)
    ext = _makie_extension()
    ext === nothing && throw(
        ArgumentError(
        "Makie is not loaded. Load CairoMakie, GLMakie, or WGLMakie first.",
    ),
    )
    return Base.invokelatest(ext.renderfig, fig)
end

next_fignum() = Base.Threads.atomic_add!(FIG_NO, 1)
reset_fignum!(n::Int = 1) = (FIG_NO[] = n)

end # module BackendHandler
