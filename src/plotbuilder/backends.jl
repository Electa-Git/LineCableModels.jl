function _backend_spec(::Val{:cairo})
    (
        extension = :LineCableModelsCairoMakieExt,
        package = "CairoMakie"
    )
end
_backend_spec(::Val{:gl}) = (
    extension = :LineCableModelsGLMakieExt,
    package = "GLMakie"
)
_backend_spec(::Val{:wgl}) = (
    extension = :LineCableModelsWGLMakieExt,
    package = "WGLMakie"
)
_backend_spec(::Val) = nothing

const FIG_NO = Base.Threads.Atomic{Int}(1)

function _makie_extension()
    return Base.get_extension(parentmodule(@__MODULE__), :LineCableModelsMakieExt)
end

function _backend_specification(backend::Val{B}) where {B}
    specification = _backend_spec(backend)
    specification === nothing && throw(ArgumentError(
        "Unknown backend :$B. Use :cairo, :gl, or :wgl.",
    ))
    return specification
end

"Return whether `backend` has been explicitly loaded."
backend_available(backend::Symbol) = backend_available(Val(backend))
function backend_available(backend::Val)
    specification = _backend_specification(backend)
    return Base.get_extension(parentmodule(@__MODULE__), specification.extension) !== nothing
end

"Return the active Makie backend as `:cairo`, `:gl`, `:wgl`, `:unknown`, or `:none`."
function current_backend_symbol()
    ext = _makie_extension()
    return ext === nothing ? :none : Base.invokelatest(ext.current_backend_symbol)
end

"""
    set_backend!(backend)

Activate an explicitly loaded Makie backend.

"""
set_backend!(backend::Symbol) = set_backend!(Val(backend))
function set_backend!(backend::Val{B}) where {B}
    specification = _backend_specification(backend)
    ext = Base.get_extension(parentmodule(@__MODULE__), specification.extension)
    ext === nothing && throw(
        ArgumentError(
        "Backend :$B is not loaded. Run `using $(specification.package)` first.",
    ),
    )
    return Base.invokelatest(ext.activate!)
end

"Activate an explicitly loaded backend when none is active."
function ensure_backend!()
    current = current_backend_symbol()
    _backend_spec(Val(current)) === nothing || return current
    throw(
        ArgumentError(
        "No Makie backend is active. Load CairoMakie, GLMakie, or WGLMakie first.",
    ),
    )
end

"Run `f` with an explicitly loaded backend and restore the previous backend."
with_backend(f::Function, backend::Symbol) = with_backend(f, Val(backend))
function with_backend(f::Function, backend::Val{B}) where {B}
    previous = current_backend_symbol()
    set_backend!(backend)
    try
        return f()
    finally
        if _backend_spec(Val(previous)) !== nothing && previous !== B
            set_backend!(Val(previous))
        end
    end
end

"Create a backend-specific screen when supported."
function make_screen(
        title::AbstractString;
        backend::Symbol = current_backend_symbol(),
        kwargs...
)
    return make_screen(Val(backend), title; kwargs...)
end

function make_screen(backend::Val, title::AbstractString; kwargs...)
    specification = _backend_specification(backend)
    ext = Base.get_extension(parentmodule(@__MODULE__), specification.extension)
    return ext === nothing ? nothing :
           Base.invokelatest(ext.make_screen, String(title); kwargs...)
end

function make_screen(backend::Symbol, title::AbstractString; kwargs...)
    return make_screen(Val(backend), title; kwargs...)
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
