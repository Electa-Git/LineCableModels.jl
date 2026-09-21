@testitem "Makie / retained units coordinates and ragged frequencies" tags=[:visual] begin
    using CairoMakie
    z=reshape(complex.(1.:12.,101.:112.),2,2,3)
    a=ObservedResult(LineParameters(z,2z,[1.,10.,100.]);length_unit=:base)
    b=ObservedResult(LineParameters(z,2z,[1.,10.,100.]),
        ((R,[2,1],[2,1],[1,3]),(X,[2,1],[2,1],[1,3]),
         (G,[2,1],[2,1],[2,3]),(B,[2,1],[2,1],[2,3]));
        length_unit=:kilo,frequency_unit=:kilo)
    options=(backend=:cairo,display_plot=false,controls=false,open_export=false)
    pages=LineCableModels.plot([a,b];ydata=(R,G),options...)
    @test length(pages)==2
    for (page,selector,scale,samples) in zip(pages,(R,G),(1.,2.),([1,3],[2,3]))
        @test length(page.axes)==4
        for (axis,(i,j)) in zip(page.axes,((1,1),(1,2),(2,1),(2,2)))
            curves=filter(p -> p isa Makie.Lines,axis.scene.plots)
            @test length(curves)==2
            expected1=Float32.(real.(z[i,j,:]).*scale)
            expected2=expected1[samples]
            @test getindex.(curves[1][1][],2)≈expected1
            @test getindex.(curves[2][1][],2)≈expected2
            @test getindex.(curves[1][1][],1)≈[1.,10.,100.]
            @test getindex.(curves[2][1][],1)≈[1.,10.,100.][samples]
        end
    end
    @test b.quantities[1].coordinates.rows==[2,1]
    @test b.quantities[1].coordinates.frequency_unit==LineCableModels.Units.units(:kilo,:hertz)
end
