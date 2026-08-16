@testitem "Extensions / backend dispatch / loaded extension activation" tags=[:extension] begin
    using Pkg
    import LineCableModels

    backends=LineCableModels.PlotBuilder.BackendHandler
    packages=(
        (:cairo, :CairoMakie, :LineCableModelsCairoMakieExt),
        (:gl, :GLMakie, :LineCableModelsGLMakieExt),
        (:wgl, :WGLMakie, :LineCableModelsWGLMakieExt)
    )
    active_dependencies=keys(Pkg.project().dependencies)
    installed=filter(entry->String(entry[2]) in active_dependencies, packages)

    if isempty(installed)
        @test backends.current_backend_symbol() === :none
    else
        for (backend, package_name, extension_name) in installed
            Core.eval(@__MODULE__, Expr(:using, Expr(:., package_name)))
            @test Base.get_extension(LineCableModels, extension_name) !== nothing
            @test backends.backend_available(backend)
            @test backends.set_backend!(backend) === backend
            @test backends.current_backend_symbol() === backend
            if backend !== :gl
                @test backends.make_screen(backend, "headless") === nothing
            end
            @test backends.with_backend(
                ()->backends.current_backend_symbol(),
                backend
            ) === backend
        end
        @test backends.ensure_backend!() in first.(installed)
    end
end
