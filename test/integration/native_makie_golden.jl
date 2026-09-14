@testitem "Makie / current rendering / asymmetric pixels and provisional export" tags=[:visual] begin
    using CairoMakie
    include(joinpath(pkgdir(LineCableModels),"test/support/golden_fixtures.jl"))
    using .GoldenFixtures
    @test_throws ArgumentError scene("unknown-scene")
    # A fresh four-corner graphic fixes the save/decode convention and rejects
    # reflection and blank frames. View acceptance needs separate calibration.
    figure=Figure(size=(120,80),backgroundcolor=:white)
    axis=Axis(figure[1,1];limits=(0,2,0,2))
    hidedecorations!(axis); hidespines!(axis)
    scatter!(axis,[.25,1.75,.25,1.75],[.25,.25,1.75,1.75];
        color=[:red,:green,:blue,:black],markersize=[7,11,15,19])
    mktempdir() do directory
        pixels=save_pixels(joinpath(directory,"first.png"),(;figure))
        repeat_pixels=save_pixels(joinpath(directory,"repeat.png"),(;figure))
        # Cairo's public colorbuffer already uses row/column image order. Check
        # that convention directly; never choose an orientation by comparison.
        buffer=Makie.colorbuffer(figure;backend=CairoMakie,px_per_unit=1)
        buffer_pixels=cat((round.(UInt8,255 .* channel.(buffer))
            for channel in (Makie.red,Makie.green,Makie.blue))...;dims=3)
        @test eltype(pixels) === UInt8
        @test size(pixels)==(80,120,3)
        @test size(buffer_pixels)==size(pixels)
        @test buffer_pixels==pixels
        @test pixel_error(pixels,repeat_pixels)==0
        @test pixel_error(pixels,reverse(pixels;dims=1))>0
        @test pixel_error(pixels,reverse(pixels;dims=2))>0
        @test pixel_error(pixels,fill(0xff,size(pixels)))>0
    end
end
