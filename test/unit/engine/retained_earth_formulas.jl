# The inherited author-output table is retired. The matrix below establishes
# supported equation/placement and cross-integrator consistency only. The
# independent real-axis controls establish their explicitly named equations.
@testitem "Engine / author equations / current placement and integration consistency" tags=[:unit] begin
    const E=LineCableModels.Engine
    controls(method)=method===:trapz ? (;max_refinements=14) :
        method===:cim ? (;samples=512,maxevals=10^6) : (;maxevals=10^6)
    rho=[Inf,100.0]; epsilon=8.8541878128e-12.*[1.0,10.0]; mu=fill(4pi*1e-7,2)
    cases=((E.EarthImpedance,:overhead,(:Carson1926,:Wise1934,:Gary1976)),
        (E.EarthImpedance,:underground,(:Pollaczek1926,:WedepohlWilcox1973,:Saad1996,:Xue2018)),
        (E.EarthImpedance,:mixed,(:Ametani2009,:Lucca1994)),
        (E.EarthImpedance,:mixedreverse,(:Ametani2009,:Lucca1994)),
        (E.EarthAdmittance,:overhead,(:Wise1948,)),
        (E.EarthAdmittance,:underground,(:Pollaczek1926,:Xue2018)))
    for (owner,placement,authors) in cases, author in authors, f in (50.0,10000.0), self in (false,true)
        placement in (:mixed,:mixedreverse) && self && continue
        @testset "$author / $placement / $f Hz / self=$self" begin
        h=placement===:overhead ? (1.0,1.5) : placement===:underground ? (-1.0,-1.5) : placement===:mixed ? (1.0,-1.5) : (-1.0,1.5)
        layers=placement===:overhead ? (1,1) : placement===:underground ? (2,2) : placement===:mixed ? (1,2) : (2,1)
        self && (h=(h[1],h[1]))
        pair=E.EarthPair(1,self ? 1 : 2,h,self ? 0.0 : 0.4,layers;radius=self ? 0.005 : nothing)
        selected=owner.Formula(author)
        binding=validate(selected,pair)
        if haskey(binding.options,:integration)
            selected=owner.Formula(author;options=(integration=(method=:quad,options=controls(:quad)),))
        end
        expected=selected(rho,epsilon,mu,2pi*im*f,pair)()
        @test isfinite(expected)
        if haskey(binding.options,:integration)
            @testset "$method" for method in (:quad,:trapz,:cim)
                actual=owner.Formula(author;options=(integration=(method,options=controls(method)),))(
                    rho,epsilon,mu,2pi*im*f,pair)()
                for component in (real,imag)
                    @test component(actual) ≈ component(expected) rtol=1e-5 atol=0
                end
            end
        end
    end
    end
    for owner in (E.EarthImpedance,E.EarthAdmittance)
        pair=E.EarthPair(1,2,(1.0,1.5),0.4,(1,1))
        @test_throws ArgumentError owner.Formula(:Pollaczek1926)(rho,epsilon,mu,100pi*im,pair)
    end
end

@testitem "Engine / Wise1948 / independent unsplit real-axis equation" tags=[:unit] begin
    using QuadGK
    const E=LineCableModels.Engine
    controls(method)=method===:trapz ? (;max_refinements=14) :
        method===:cim ? (;samples=512,maxevals=10^6) : (;maxevals=10^6)
    epsilon=8.8541878128e-12.*[1.0,10.0];mu=fill(4pi*1e-7,2);rho=[Inf,100.0]
    for f in (50.0,10000.0),self in (true,false)
        heights=self ? (1.0,1.0) : (1.0,1.5)
        y=self ? 0.0 : .4;H=sum(heights)
        distance=self ? .005 : hypot(y,heights[1]-heights[2])
        references=[setprecision(BigFloat,bits) do
            s=2big(pi)*im*f
            e0=big"8.8541878128e-12";u0=4big(pi)*big"1e-7"
            g0=s^2*u0*e0;g1=s*u0*(big".01"+s*10*e0)
            ratio=g1/g0;scale=sqrt(g1-g0)
            cutoff=big(40)/H
            kernel=lambda->exp(-H*lambda)*cos(y*lambda)/(ratio*lambda+sqrt(lambda^2+g1-g0))
            integral,estimate=quadgk(kernel,big"0",abs(scale/ratio),abs(scale),inv(big(H)),cutoff;
                rtol=big"1e-12",maxevals=10^6)
            # Re(ratio)=10 and Re(sqrt(lambda²+g1-g0))>0, hence
            # |denominator|>=10lambda: this bounds the omitted real-axis tail.
            tail=exp(-H*cutoff)/(10H*cutoff)
            value=(log(hypot(big(y),big(H))/big(distance))+2integral)/(2big(pi)*e0)
            (value=value,bound=(estimate+tail)/(big(pi)*e0))
        end for bits in (128,256,512)]
        expected=last(references)
        uncertainty=expected.bound+abs(expected.value-references[2].value)
        pair=E.EarthPair(1,self ? 1 : 2,heights,y,(1,1);radius=self ? .005 : nothing)
        @testset "$f Hz / self=$self / $method" for method in (:quad,:trapz,:cim)
            actual=E.EarthAdmittance.Formula(:Wise1948;options=(integration=(method,options=controls(method)),))(
                rho,epsilon,mu,2pi*im*f,pair)()
            for component in (real,imag)
                budget=1e-6*abs(component(expected.value))
                @test uncertainty<=budget/4
                @test abs(component(actual-expected.value))+uncertainty<=budget
            end
        end
    end
end

@testitem "Engine / Carson / independent conductive real-axis integral" tags=[:unit] begin
    using QuadGK
    const E=LineCableModels.Engine
    controls(method)=method===:trapz ? (;max_refinements=14) :
        method===:cim ? (;samples=512,maxevals=10^6) : (;maxevals=10^6)
    rho=[Inf,100.0]; epsilon=8.8541878128e-12.*[1,10]; mu=fill(4pi*1e-7,2)
    for f in (50.0,10000.0), self in (true,false)
        s=2pi*im*f; heights=self ? (1.0,1.0) : (1.0,1.5)
        distance=self ? .005 : hypot(.4,.5)
        y=self ? 0.0 : .4; h=sum(heights); image_distance=hypot(h,y)
        refs=[setprecision(BigFloat,bits) do
            sb=2big(pi)*im*f;mub=4big(pi)*big"1e-7"
            k2=sb*mub/big"100";hb=BigFloat(h);yb=BigFloat(y)
            cutoff=big"40"/hb
            kernel(t)=exp(-hb*t)*cos(yb*t)/(t+sqrt(t^2+k2))
            value,estimate=quadgk(kernel,big"0",abs(sqrt(k2)),inv(hb),cutoff;
                rtol=bits==128 ? big"1e-8" : bits==256 ? big"1e-10" : big"1e-12",maxevals=10^6)
            # Re sqrt(t²+i a)>0, so |denominator|>=t.
            tail=exp(-hb*cutoff)/(hb*cutoff)
            expected=sb*mub/(2big(pi))*(log(hypot(hb,yb)/BigFloat(distance))+2value)
            (value=expected,bound=abs(sb*mub/big(pi))*(estimate+tail))
        end for bits in (128,256,512)]
        expected=last(refs).value
        uncertainty=last(refs).bound+abs(last(refs).value-refs[2].value)
        pair=E.EarthPair(1,self ? 1 : 2,heights,y,(1,1);radius=self ? .005 : nothing)
        @testset "$f Hz / self=$self / $method" for method in (:quad,:trapz,:cim)
            actual=E.EarthImpedance.Formula(:Carson1926;options=(integration=(method,options=controls(method)),))(
                rho,epsilon,mu,s,pair)()
            for component in (real,imag)
                budget=1e-6*abs(component(expected))
                @test uncertainty<=budget/4
                @test abs(component(actual)-component(expected))+uncertainty<=budget
            end
        end
    end
end
