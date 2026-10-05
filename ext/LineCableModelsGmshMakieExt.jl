module LineCableModelsGmshMakieExt

import LineCableModels
import Gmsh
import Makie
import GeometryBasics
import InteractiveUtils
import LineCableModels.DataModel
import LineCableModels.ImportExport
import LineCableModels.PlotBuilder: plot
using DocStringExtensions: TYPEDSIGNATURES
using Printf: @sprintf
using Statistics: mean
using Colors: HSV, RGB, RGBA
using Makie: Auto, Axis, Button, DataAspect, Fixed, GridLayout, Label, Legend,
             Observable, colgap!, colsize!, hlines!, lines!, on, poly!, rowsize!,
             scatter!, translate!, with_theme
import Base: resize!

const FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
const Renderer = Base.get_extension(LineCableModels, :LineCableModelsMakieExt)
using .FEM: FEMMesh, FEMFieldMap
import .Renderer: _mesh_input, _spatial_layer!, _mesh_sidebar!, _mesh_controls!
using .Renderer: _addon_activate_backend, _addon_compose_guides!, _addon_detach!, _addon_display!,
                 _addon_finish!, _addon_panel!, _addon_release_frames!,
                 _addon_reset!, _addon_shell, _addon_theme, _addon_widget!,
                 _native_system_shapes

include("LineCableModelsGmshExt/plotting/fem.jl")

end
