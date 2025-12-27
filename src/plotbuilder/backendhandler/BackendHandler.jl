"""
Makie backend handler for LineCableModels.

Design goals:
- Precompile in any environment (never touch GL/WGL at load-time).
- Default to CairoMakie as a safe backend.
- Offer a single `set_backend!` API (no user `using` needed).
- In headless (e.g., Literate → Documenter), display PNG inline.

How to wire it (minimal integration):

1) Add this helper to your package module once:

   include(joinpath(@__DIR__, "..", "@INPROGRESS", "makie_backend_alt.jl"))
   using .BackendHandler

2) In Makie-based preview entrypoints, ensure a backend is active:

   # Use caller keyword `backend::Union{Nothing,Symbol}` if you keep it
   BackendHandler.ensure_backend!(backend === nothing ? :cairo : backend)

3) For GL interactive windows without directly referencing GLMakie:

   if BackendHandler.current_backend_symbol() == :gl
	   if (scr = BackendHandler.gl_screen("Title")) !== nothing
		   display(scr, fig)
	   else
		   display(fig)
	   end
   else
	   display(fig)
   end

4) For docs/headless builds (Literate/Documenter) PNG inline display:

   # in your final display branch
   BackendHandler.renderfig(fig)

5) Let users select backends interactively (no extra imports):

   BackendHandler.set_backend!(:gl)   # or :wgl, :cairo

Notes:
- No `@eval import` anywhere; backends are loaded via `Base.require` using PkgId.
- Calls into newly loaded modules go through `Base.invokelatest` to avoid
  world-age issues.
"""
module BackendHandler

using Makie
using UUIDs
using ...Utils: is_headless

export set_backend!,
	ensure_backend!, current_backend_symbol, renderfig, with_backend, make_screen

# ---------------------------------------------------------------------------
# Backend registry
# ---------------------------------------------------------------------------

const _BACKENDS = Dict{Symbol, Tuple{UUID, String}}(
	:cairo => (UUID("13f3f980-e62b-5c42-98c6-ff1f3baf88f0"), "CairoMakie"),
	:gl    => (UUID("e9467ef8-e4e7-5192-8a1a-b1aee30e663a"), "GLMakie"),
	:wgl   => (UUID("276b4fcb-3e11-5398-bf8b-a0c2d153d008"), "WGLMakie"),
)

_pkgid(sym::Symbol) = begin
	tup = get(_BACKENDS, sym, nothing)
	tup === nothing && throw(
		ArgumentError("Unknown backend: $(sym). Valid: $(collect(keys(_BACKENDS)))"),
	)
	Base.PkgId(tup[1], tup[2])
end

"""Return true if a backend package exists in the environment."""
backend_available(backend::Symbol) = Base.find_package(last(_BACKENDS[backend])) !== nothing

# Track the last activated backend symbol (separate from Makie.internal state)
const _active_backend = Base.RefValue{Symbol}(:none)

# ---------------------------------------------------------------------------
# Activation core (lazy, world-age safe)
# ---------------------------------------------------------------------------

function _activate_backend!(backend::Symbol; allow_interactive_in_headless::Bool = false)
	if is_headless() && backend != :cairo && !allow_interactive_in_headless
		@warn "Headless environment: forcing :cairo instead of $(backend)."
		return _activate_backend!(:cairo; allow_interactive_in_headless)
	end

	pid = _pkgid(backend)
	# Load the backend module into Julia's module world; idempotent if already loaded
	mod = Base.require(pid)
	# Call `activate!` safely with world-age correctness
	Base.invokelatest(getproperty(mod, :activate!))
	_active_backend[] = backend
	return backend
end

"""Ensure a backend is active. Defaults to :cairo the first time."""
function ensure_backend!(backend::Union{Nothing, Symbol} = nothing)
	if backend === nothing
		return _active_backend[] == :none ? _activate_backend!(:cairo) : _active_backend[]
	else
		return set_backend!(backend)
	end
end

"""Activate a specific backend (:cairo, :gl, :wgl).

In headless, :gl/:wgl requests fall back to :cairo unless `force=true`.
No `using` required by callers.
"""
function set_backend!(backend::Symbol; force::Bool = false)
	haskey(_BACKENDS, backend) || throw(
		ArgumentError("Unknown backend: $(backend). Valid: $(collect(keys(_BACKENDS)))"),
	)
	if backend != :cairo && !backend_available(backend)
		if !is_headless()
			@warn "Backend $(last(_BACKENDS[backend])) not in environment; using :cairo."
		end

		return _activate_backend!(:cairo)
	end
	return _activate_backend!(backend; allow_interactive_in_headless = force)
end

"""Symbol of the current Makie backend (:cairo, :gl, :wgl, :unknown, :none)."""
function current_backend_symbol()
	try
		nb = nameof(Makie.current_backend())
		nb === :CairoMakie && return :cairo
		nb === :GLMakie && return :gl
		nb === :WGLMakie && return :wgl
		return :unknown
	catch
		return :none
	end
end

"""
	with_backend(f, backend; force=false)

Temporarily activate `backend`, run `f()`, then restore the previous backend.

Intended for export workflows (switch to CairoMakie, render/save, restore).
Restoration is best-effort: if the previous backend can’t be determined, it’s skipped.
"""
function with_backend(f::Function, backend::Symbol; force::Bool = false)
	prev = current_backend_symbol()
	prev_ok =
		prev in (:cairo, :gl, :wgl) ? prev :
		(_active_backend[] in (:cairo, :gl, :wgl) ? _active_backend[] : :none)

	ensure_backend!(backend)

	try
		return f()
	finally
		if prev_ok != :none && prev_ok != current_backend_symbol()
			try
				set_backend!(prev_ok; force = force)
			catch e
				@warn "Failed to restore backend $(prev_ok)" exception=(
					e,
					catch_backtrace(),
				)
			end
		end
	end
end

"""
	make_screen(title; backend=current_backend_symbol(), kwargs...) -> screen | nothing

Create a backend-specific screen/window handle for interactive display.

Currently:
- `:gl` returns `GLMakie.Screen(; title=...)` (loaded lazily).
- other backends return `nothing`.
"""
function make_screen(title::AbstractString;
	backend::Symbol = current_backend_symbol(),
	kwargs...,
)
	backend == :gl || return nothing
	mod = Base.require(_pkgid(:gl))
	ctor = getproperty(mod, :Screen)
	return Base.invokelatest(ctor; title = String(title), kwargs...)
end

"""
	make_screen(backend, title; kwargs...) -> screen | nothing

Convenience overload.
"""
make_screen(backend::Symbol, title::AbstractString; kwargs...) =
	make_screen(title; backend = backend, kwargs...)



"""Display a figure appropriately in headless docs or interactive sessions.

- Headless: returns `DisplayAs.Text(DisplayAs.PNG(fig))` if `DisplayAs` exists;
			otherwise attempts to rasterize via CairoMakie and returns nothing.
- Interactive: calls `display(fig)` and returns its result.
"""
function renderfig(fig)
	if is_headless()
		try
			D = Base.require(
				Base.PkgId(UUID("0b91fe84-8a4c-11e9-3e1d-67c38462b6d6"), "DisplayAs"),
			)
			return D.Text(D.PNG(fig))
		catch
			try
				ensure_backend!(:cairo)
				cm = Base.require(_pkgid(:cairo))
				savef = getproperty(cm, :save)
				io = IOBuffer()
				Base.invokelatest(savef, io, fig)
				return nothing
			catch
				return nothing
			end
		end
	else
		return display(fig)
	end
end

const FIG_NO = Base.Threads.Atomic{Int}(1)
next_fignum() = Base.Threads.atomic_add!(FIG_NO, 1)
reset_fignum!(n::Int = 1) = (FIG_NO[] = n)

# ---------------------------------------------------------------------------
# Default backend at runtime load (safe: CairoMakie only)
# ---------------------------------------------------------------------------

function __init__()
	try
		# Only set a default if nothing looks active yet
		if current_backend_symbol() in (:none, :unknown)
			ensure_backend!(:cairo)
		end
	catch e
		@warn "Failed to initialize default CairoMakie" exception=(e, catch_backtrace())
	end
end

end # module

