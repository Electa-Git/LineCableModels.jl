# Compatibility helpers for manual voltage-path diagnostics.
using LineCableModels, Gmsh
quasi_full_path_points(args...; kwargs...) =
    Base.get_extension(LineCableModels, :LineCableModelsGmshExt)._voltage_path_points(args...; kwargs...)
function write_quasi_full_paths(path, mesh, plan, model, positions, radius; kwargs...)
    Base.get_extension(LineCableModels, :LineCableModelsGmshExt)._write_voltage_paths(
        path, mesh, plan, model; endpoints=[(x, y-radius) for (x,y) in positions], kwargs...)
end
