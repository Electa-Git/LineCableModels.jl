@testitem "Engine / physical cutoffs shared by observations and comparisons" tags=[:unit] begin
    using LinearAlgebra: diag
    using LineCableModels.Engine: compare
    const GR=LineCableModels.Grammar
    f=[1.,1e3,1e7]
    z=fill(1.0+im,2,2,3)
    y=fill(1e-13+1e-18im,2,2,3)
    y[2,2,:].=1e-4 .+ 2π.*f.*1e-9im
    source=LineParameters(z,y,f)
    saved=deepcopy((z,y,f))
    selected=((G,1,1,:),(B,1,1,:))
    a=ObservedResult(source,selected;length_unit=:base)
    b=ObservedResult(source,selected;length_unit=:kilo)
    @test all(iszero,observe(a,G))
    @test observe(b,G)==1000observe(a,G)
    @test observe(ObservedResult(source,selected;clip=false,length_unit=:base),G)==real.(y[1,1,:])
    @test all(iszero,observe(a,B))
    polar=ObservedResult(source,((Y,abs,1,1,:),(Y,angle,1,1,:)))
    @test all(iszero,observe(polar,Y,abs))
    @test all(ismissing,observe(polar,Y,angle))
    @test observe(ObservedResult(source,selected;atol=(G=0.,B=0.),length_unit=:base),G)==real.(y[1,1,:])
    @test observe(ObservedResult(source),G)[1,1,:]==observe(b,G)
    diagonal=ObservedResult(source,((G,diag,:,:),(B,diag,:,:));length_unit=:base)
    @test observe(diagonal,G,diag)[1,:]==observe(a,G)
    @test observe(a,G,1,1,2:3)==observe(a,G)[2:3]
    product=GR.observation_product(a,G)
    @test product.thresholds.kind===:declared_floor
    @test count(product.engineering_zero)==3
    @test_throws ArgumentError ObservedResult(source,selected;atol=1e-12)
    @test_throws ArgumentError ObservedResult(source,selected;atol=(oops=0.,))
    for q in (R,X,L,G,B,C,Z,Y)
        comparison=compare(source,source,q)
        @test comparison.details.data.resolution.unit==LineCableModels.Units.native_unit(q,:pul)
        @test !haskey(comparison.details.data.resolution,:revision)
    end
    @test compare(source,source,B).details.data.atol≈2π.*f.*compare(source,source,C).details.data.atol
    @test compare(source,source,X).details.data.atol≈2π.*f.*compare(source,source,L).details.data.atol
    @test (z,y,f)==saved
end

@testitem "Engine / cutoff boundaries and unassessed components" tags=[:unit] begin
    using LineCableModels.Engine: compare
    const GR=LineCableModels.Grammar
    for T in (Float32,Float64,BigFloat)
        cutoff=T(1e-12)
        f=T[1,1e3,1e7]
        g=reshape(T[cutoff/2,cutoff,2cutoff],1,1,:)
        z=fill(complex(one(T),one(T)),1,1,3)
        y=complex.(g,zero(T))
        a=LineParameters(z,y,f);b=LineParameters(z,2y,f)
        base=observe(ObservedResult(a;length_unit=:base,atol=(G=cutoff,)),G)
        @test eltype(base)===T
        @test vec(base)==T[0,0,2cutoff]
        @test observe(ObservedResult(a;atol=(G=cutoff,)),G)≈1000base
        for normalization in (:reference_rms,:pointwise), (left,right) in ((a,b),(b,a))
            error=compare(left,right,G;normalization,atol=cutoff)
            @test ismissing(only(error.relative))
            @test ismissing(only(error.absolute))
            @test Base.nonmissingtype(eltype(error.absolute))===T
            @test only(error.details.data.unresolved_samples).reference>0
            @test only(error.details.data.unresolved_samples).candidate>0
        end
        @test !ismissing(only(compare(a,b,G;band=:wide,atol=cutoff).relative))
        total=LineParameters(10z,10y,f;basis=:total)
        @test_throws ArgumentError ObservedResult(total)
        @test observe(ObservedResult(total;atol=(R=0.,X=0.,G=10cutoff,B=0.)),G)≈10base
        @test compare(a,b,G;atol=0).relative≈fill(one(T),1,1)
        small=fill(T(1e-30),1,1,3)
        @test only(compare(small,2small).absolute)>0
        @test only(compare(small,2small).relative)≈one(T)
    end
    f=[1.,10.,100.]
    z=fill(1.0+im,2,2,3);y=fill(2e-12+1e-6im,2,2,3)
    a=LineParameters(z,copy(y),f);b=LineParameters(z,copy(y),f)
    b.Y.values[2,2,:].=1e20
    @test observe(ObservedResult(a),G,1,1,:)==observe(ObservedResult(b),G,1,1,:)
    standalone=ShuntAdmittance(fill(1e-18im,2,2,3))
    raw=ObservedResult(standalone;length_unit=:base)
    @test observe(raw,B)==fill(1e-18,2,2,3)
    @test GR.observation_product(raw,B).thresholds.kind===:unassessed
    @test all(iszero,observe(ObservedResult(standalone;frequencies=f),B))
    @test all(iszero,observe(ObservedResult(standalone;atol=(B=1e-12,)),B))
    @test_throws ArgumentError ObservedResult(a;frequencies=[1.,2.,3.])
    @test_throws DimensionMismatch ObservedResult(standalone;frequencies=[1.])
    for invalid in (-1,Inf,NaN,true,"bad",(G=-1.,))
        @test_throws ArgumentError ObservedResult(a;atol=invalid)
    end
    a.Y.values[1,1,1]=NaN+im
    @test ismissing(observe(ObservedResult(a),G,1,1,1))
    @test ismissing(compare(a,b,G).absolute[1,1])
    @test !ismissing(compare(a,b,G;band=(10.,100.)).absolute[1,1])
    dc=ObservedResult(LineParameters(z,y,[0.,10.,100.]),(R,L,G,C))
    @test ismissing(observe(dc,C,1,1,1))
    @test all(isfinite,observe(dc,C,1,1,2:3))
    single=LineParameters(ones(ComplexF32,1,1,3),fill(ComplexF32(1e-6),1,1,3),f)
    double=LineParameters(ones(ComplexF64,1,1,3),fill(complex((Float64(Float32(1e-12))+1e-12)/2),1,1,3),f)
    @test all(iszero,observe(ObservedResult(double),G))
    @test ismissing(only(compare(single,double,G).relative))
    @test ismissing(only(compare(double,single,G).relative))
    @test all(iszero,observe(ObservedResult(single;atol=(G=1e100,)),G))
end

@testitem "UQ / recentering preserves every dependency and spread" tags=[:extension] begin
    using Measurements
    using LineCableModels.Grammar: detach,observation_resolution
    summary=LineCableModels.UQ.SampleSummary([1e-18,2e-18,3e-18])
    detached=detach(summary,1000.)
    @test detached.std==1000summary.std>0
    @test detached.min<=detached.q05<=detached.median<=detached.q95<=detached.max
    @test detached.mean==1000summary.mean
    @test detached.n==summary.n
    for T in (Float32,Float64,BigFloat)
        cutoff=T(1e-12)
        original=measurement.(T[cutoff/2,cutoff/2,2cutoff,2cutoff],T[cutoff/2,2cutoff,cutoff/2,2cutoff])
        y=reshape(complex.(original,zero.(original)),1,1,:)
        source=LineParameters(copy(y),y,T[1,10,100,1000])
        for length_unit in (:base,:kilo)
            factor=length_unit===:base ? T(1) : T(1000)
            observed=ObservedResult(source;length_unit,atol=(G=cutoff,))
            values=vec(observe(observed,G))
            @test eltype(values)===eltype(original)
            @test Measurements.value.(values)≈factor.*T[0,0,2cutoff,2cutoff]
            @test Measurements.uncertainty.(values)≈factor.*Measurements.uncertainty.(original)
            @test all(iszero,Measurements.uncertainty.(values.-factor.*original))
            raw=vec(observe(ObservedResult(source;length_unit,clip=false),G))
            @test Measurements.value.(raw)==factor.*Measurements.value.(original)
            @test all(iszero,Measurements.uncertainty.(raw.-factor.*original))
        end
    end
    value=complex(measurement(1e-18,2e-12),measurement(1e-18,3e-12))
    resolution=observation_resolution(value,Y;atol=(G=1e-12,B=1e-12))
    @test resolution.unresolved && resolution.available
    observed=ObservedResult(LineParameters(fill(value,1,1,1),fill(value,1,1,1),[1.]),
        ((Y,abs),(Y,angle));atol=(G=1e-12,B=1e-12),length_unit=:base)
    @test ismissing(only(observe(observed,Y,angle)))
    @test nominal(only(observe(observed,Y,abs)))==0
    @test uncertainty(only(observe(observed,Y,abs)))>0
    for value in (NaN,Inf,missing)
        r=observation_resolution(Union{Missing,Float64}[value],G)
        @test !only(r.unresolved) && !only(r.available)
    end
end
