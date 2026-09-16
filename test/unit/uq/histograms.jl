@testitem "UQ / histogram bins preserve finite constant and narrow populations" tags=[:unit] begin
    using Statistics
    for values in ([3.0,3.0,3.0], [1.0,1.0,nextfloat(1.0)],
            [1.0,nextfloat(1.0),nextfloat(1.0,2)],
            [0.0,nextfloat(0.0)], [floatmax(Float64),floatmax(Float64)],
            [-floatmax(Float64),floatmax(Float64)],
            Float32[1,1,nextfloat(1.0f0)])
        original=copy(values)
        for bins in (nothing,1,2,32)
            density=HistogramDensity(values;bins)
            edges=density.edges
            @test all(isfinite,edges)
            @test all(>(0),diff(edges))
            @test first(edges)<=minimum(values)<=maximum(values)<=last(edges)
            @test all(isfinite,density.density)
            @test sum(density.density .* diff(edges))≈1
            # Count the original population in each documented half-open bin;
            # the final bin alone includes the right endpoint.
            expected=[count(x->edges[i]<=x &&
                (i==length(edges)-1 ? x<=edges[i+1] : x<edges[i+1]),values)
                for i in 1:length(edges)-1] ./ length(values)
            @test density.density .* diff(edges)≈expected
            @test values==original
            if minimum(values)!=maximum(values)
                @test first(edges)==minimum(values)
                @test last(edges)==maximum(values)
            else
                @test length(density.density)==1
            end
        end
    end
    constant=HistogramDensity([1.,1.,1.];bins=32)
    narrow=HistogramDensity([1.,1.,nextfloat(1.)];bins=32)
    @test constant.edges!=narrow.edges
    @test_throws ArgumentError HistogramDensity([NaN])
    @test_throws ArgumentError HistogramDensity([Inf])
    @test_throws ArgumentError HistogramDensity(Float64[])
    @test_throws ArgumentError HistogramDensity([1.];bins=0)
end
