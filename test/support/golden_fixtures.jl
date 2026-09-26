# Fresh scene recipes. Outputs are provisional repeatability records until a
# person accepts their stated validation scope. No saved image is an oracle.
module GoldenFixtures
using LineCableModels
using LineCableModels.ReportBuilder: BenchmarkTableDefinition
using CairoMakie
using Measurements: measurement
using Statistics: mean
include("scenarios.jl")
using .CurrentScenarios
export scene_names, scene, custom_layout_plot, save_pixels, pixel_error, check_scene, defect!
const scene_names=("line_rlcg","line_zy_cartesian","line_zy_polar",
    "formulation_comparison","uq_comparison","uncertainty_intervals",
    "cable_preview","cable_preview_compact","system_preview","material_scale",
    "mc_hist","mc_pdf","mc_ecdf","mc_qq","custom_layout")
const plot_options=(backend=:cairo,display_plot=false,controls=false,open_export=false,
    size=(1100,750))

function custom_layout_plot(;backend=:cairo,display_plot=false,controls=false)
    return LineCableModels.plotwindow(;backend,display_plot,controls,open_export=false,
            title="Asymmetric current layout",size=(1100,750),export_name="custom_layout") do canvas
        left=Axis(canvas[1,1];title="Quadratic",xlabel="position",ylabel="response")
        right=Axis(canvas[1,2];title="Discrete samples",xlabel="index",ylabel="response")
        x=[.25,.75,1.5,2.5]
        lines!(left,x,x.^2 .+ .3;label="x² + 0.3",color=:navy,linewidth=3)
        scatter!(right,[1,2,3],[.4,1.1,.7];label="three samples",color=:orange,markersize=18)
        axislegend(left;position=:lt)
        axislegend(right;position=:rb)
        colsize!(canvas,1,Relative(.65))
        nothing
    end
end

function actual_mc(;trials=128,seed=2027)
    design=coaxial_design()
    space=Gridspace{CableConstantsProblem}(t->CableConstantsProblem(design;temperature=t),
        (Grid(20.0,AbsoluteError(1.0)),))
    return compute(ParametricProblem(space),MonteCarlo(CableConstantsFormulation();
        trials,seed,distribution=:uniform,return_samples=true,return_histograms=true))
end

function comparison(;uncertain=false)
    if uncertain
        build_problem=scale->begin
            design=coaxial_design(;scale)
            system=build(LineCableSystem,[design,design],[(0.,-1.),(.2,-1.3)];
                connections=[Dict(:core=>1,:sheath=>0),Dict(:core=>2,:sheath=>0)])
            LineParametersProblem(system;temperature=20.,frequencies=[10.,100.,1000.],
                earth_props=homogeneous(rho=100.))
        end
        space=Gridspace{LineParametersProblem}(build_problem,(Grid(1.,AbsoluteError(.01)),))
        problem=ParametricProblem(space)
        selection=Formulation(options=(reduce_bundle=false,kron_reduction=true,ideal_transposition=false))
        run(seed)=compute(problem,MonteCarlo(selection;trials=128,seed,
            distribution=:uniform,return_samples=true,return_histograms=false))
        reference=run(2027)
        candidate=run(2039)
        metadata=(port_order=["west","east"],)
        artifact=report(BenchmarkTableDefinition(((statistics,R,mean),);bands=(:all,)),
            (reference=(result=reference,metadata=metadata),
                candidate=(result=candidate,metadata=metadata)))
        return LineCableModels.plot(artifact;ydata=(R,),plot_options...)
    end
    reference=two_conductor_results()
    omega=reshape(2pi.*reference.f,1,1,:)
    candidate=LineParameters(1.1R(reference).+im.*omega.*1.2L(reference),
        1.3G(reference).+im.*omega.*1.4C(reference),reference.f;
        details=ComputationDetails(;coordinates=["west","east"],))
    selections=[NamedTuple(Formulation()),NamedTuple(Formulation(earth_impedance=:pollaczek1926))]
    points=ParametricResult(nothing,[reference,candidate],
        (problems=[:current],formulations=selections), ComputationDetails((;)))
    baseline=(result=reference,metadata=(port_order=["west","east"],formulation=selections[1],axes=nothing))
    artifact=report(BenchmarkTableDefinition((R,);bands=(:all,)),(reference=baseline,candidate=points))
    return LineCableModels.plot(artifact;ydata=(R,),plot_options...)
end

# Dispatch is deliberately lazy: selecting one view never constructs another.
function scene(name::AbstractString)
    name in scene_names || throw(ArgumentError("unknown provisional view: $name"))
    CairoMakie.activate!()
    set_theme!(backgroundcolor=:white, fonts=(regular="DejaVu Sans",bold="DejaVu Sans"))
    name=="custom_layout" && return custom_layout_plot()
    name=="material_scale" && return show_material_scale(;plot_options...)
    name=="formulation_comparison" && return comparison()
    name=="uq_comparison" && return comparison(;uncertain=true)
    if startswith(name,"mc_")
        result=actual_mc()
        verb=name=="mc_hist" ? Makie.hist : name=="mc_pdf" ? Makie.stairs :
            name=="mc_ecdf" ? Makie.ecdfplot : Makie.qqplot
        return verb(result,R;plot_options...,length_unit=:base,quantity_units=:base)
    end
    if name in ("cable_preview","cable_preview_compact","system_preview")
        design=coaxial_design()
        if name=="system_preview"
            system=build(LineCableSystem,[design,coaxial_design(;scale=1.2,name="larger")],
                [(0.0,-.3),(.05,-.4)];connections=[Dict(:core=>1,:sheath=>0),Dict(:core=>2,:sheath=>0)])
            handle=preview(system;earth_model=homogeneous(rho=100.0),plot_options...)
            ylims!(only(handle.axes),-.45,.05)
            return handle
        end
        options=name=="cable_preview_compact" ? merge(plot_options,(size=(900,350),)) : plot_options
        return preview(design;options...)
    end
    parameters=two_conductor_results()
    if name=="uncertainty_intervals"
        variants=map((1.0,2.0,3.0)) do width
            f=measurement.(parameters.f,width*.01parameters.f)
            scale=[width*.01(i+j+k) for i in 1:2,j in 1:2,k in 1:3]
            z=complex.(measurement.(real.(Z(parameters)),scale.*real.(Z(parameters))),
                measurement.(imag.(Z(parameters)),scale.*imag.(Z(parameters))))
            y=complex.(measurement.(real.(Y(parameters)),scale.*real.(Y(parameters))),
                measurement.(imag.(Y(parameters)),scale.*imag.(Y(parameters))))
            LineParameters(z,y,f;details=ComputationDetails(;coordinates=["west","east"],))
        end
        return LineCableModels.plot(variants...;ydata=(R,),plot_options...,
            length_unit=:base,quantity_units=:base,freq_unit=:base,clip=false,
            series_labels=("one sigma","two sigma","three sigma"))
    end
    requests=name=="line_rlcg" ? (R,L,G,C) : name=="line_zy_polar" ? (abs,angle) : (Z,Y)
    return LineCableModels.plot(parameters,requests;plot_options...,
        length_unit=:base,quantity_units=:base,freq_unit=:base,clip=false)
end

function save_pixels(path,handle)
    handle.figure.scene.backgroundcolor[]=Makie.RGBAf(1,1,1,1)
    CairoMakie.save(path,handle.figure;px_per_unit=1)
    # Both sides use the same PNG decoder. No orientation or alignment search.
    pixels=CairoMakie.FileIO.load(path)
    return cat((round.(UInt8,255 .* channel.(pixels)) for channel in (Makie.red,Makie.green,Makie.blue))...;dims=3)
end
function pixel_error(a,b)
    size(a)==size(b) || return Inf
    return maximum(abs.(Int16.(a).-Int16.(b)))
end
include("rendering_checks.jl")
end
