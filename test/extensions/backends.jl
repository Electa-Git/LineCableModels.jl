@testitem "Extensions / backend dispatch / loaded extension activation" tags=[:extension] begin
    import LineCableModels

    backends = LineCableModels.PlotBuilder
    packages=(
        (:cairo, :CairoMakie, :LineCableModelsCairoMakieExt),
        (:gl, :GLMakie, :LineCableModelsGLMakieExt),
        (:wgl, :WGLMakie, :LineCableModelsWGLMakieExt)
    )
    installed=filter(entry->Base.find_package(String(entry[2])) !== nothing, packages)

    if isempty(installed)
        @test backends.current_backend_symbol() === :none
    else
        for (backend, package_name, extension_name) in installed
            Core.eval(@__MODULE__, Expr(:using, Expr(:., package_name)))
            @test Base.get_extension(LineCableModels, extension_name) !== nothing
            @test backends.backend_available(backend)
            @test backends.backend_available(Val(backend))
            @test backends.set_backend!(backend) === backend
            @test backends.set_backend!(Val(backend)) === backend
            @test backends.current_backend_symbol() === backend
            if backend !== :gl
                @test backends.make_screen(backend, "headless") === nothing
                @test backends.make_screen(Val(backend), "headless") === nothing
            end
            @test backends.with_backend(
                ()->backends.current_backend_symbol(),
                backend
            ) === backend
            @test backends.with_backend(
                ()->backends.current_backend_symbol(),
                Val(backend)
            ) === backend
            @test_throws MethodError backends.set_backend!(backend; force = true)
            @test_throws MethodError backends.with_backend(identity, backend; force = true)
        end
        if length(installed) >= 2
            first_backend=installed[1][1]
            second_backend=installed[2][1]
            backends.set_backend!(first_backend)
            @test backends.with_backend(
                ()->backends.current_backend_symbol(),
                second_backend
            ) === second_backend
            @test backends.current_backend_symbol() === first_backend
        end
        @test backends.ensure_backend!() in first.(installed)
    end
end
