# Current API gallery inputs shared with maintained rendering tests.
using LineCableModels
include(joinpath(@__DIR__,"../../../test/support/scenarios.jl"))
using .CurrentScenarios
function build_manual_plot_gallery(backend=:gl;display_plot=true,export_theme=:default)
    options=(;backend,display_plot,export_theme)
    parameters=two_conductor_results()
    gallery=Pair{String,UIPlot}[]
    for (title,requests) in (("RLCG",(R,L,G,C)),("Cartesian ZY",(Z,Y)),("Polar ZY",(abs,angle)))
        pages=LineCableModels.plot(parameters,requests;options...)
        append!(gallery,["$title $index"=>page for (index,page) in enumerate(pages)])
    end
    for (index,scale) in enumerate((1.0,1.2))
        push!(gallery,"Current cable $index"=>preview(coaxial_design(;scale);options...))
    end
    push!(gallery,"Current system"=>preview(three_phase_system();earth_model=homogeneous(rho=100.0),options...))
    push!(gallery,"Materials"=>show_material_scale(;options...))
    return gallery
end
