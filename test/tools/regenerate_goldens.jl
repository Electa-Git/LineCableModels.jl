# Validate selection before importing a renderer or constructing any scene.
get(ENV,"LINECABLEMODELS_UPDATE_PLOT_REFERENCES","false")=="true" ||
    error("set LINECABLEMODELS_UPDATE_PLOT_REFERENCES=true to generate provisional images")
requested=get(ENV,"LINECABLEMODELS_PLOT_REFERENCE","")
names=("line_rlcg","line_zy_cartesian","line_zy_polar","formulation_comparison",
    "uq_comparison","uncertainty_intervals","cable_preview","cable_preview_compact",
    "system_preview","material_scale","mc_hist","mc_pdf","mc_ecdf","mc_qq","custom_layout")
requested in names || throw(ArgumentError("select one current provisional view; unknown selection: $requested"))
using SHA, TOML
include(joinpath(@__DIR__,"..","support","golden_fixtures.jl"))
using .GoldenFixtures
output=mktempdir(;prefix="lcm-rendering-provisional-",cleanup=false)
println("Provisional output: ",output)
open(joinpath(output,"identity.toml"),"w") do io
    TOML.print(io,Dict("status"=>"provisional","view"=>requested,"julia"=>string(VERSION),
        "project"=>Base.active_project(),"threads"=>Threads.nthreads(),
        "recipe_sha256"=>bytes2hex(open(sha256,joinpath(@__DIR__,"..","support","golden_fixtures.jl"))),
        "source_commit"=>readchomp(`git -C $(dirname(dirname(@__DIR__))) rev-parse HEAD`)))
end
constructed=scene(requested)
handles=constructed isa AbstractVector ? constructed : [constructed]
for (index,handle) in enumerate(handles)
    name=length(handles)==1 ? requested : requested*"-page$index"
    save_pixels(joinpath(output,name*".png"),handle)
end
check_scene(requested,handles)
println("Saved ",length(handles)," current pages; all pages belong to the named view")
println("Saved provisional view: ",requested)
