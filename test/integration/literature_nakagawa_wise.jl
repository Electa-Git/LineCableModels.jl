@testitem "Engine / literature assimilation / Nakagawa final coefficients and Wise SI reduction" begin
    using LineCableModels,QuadGK
    E=LineCableModels.Engine;EI=E.EarthImpedance
    mu0=4π*1e-7;eps0=8.8541878128e-12
    for frequency in (1.,50.,1e5,1e7),permeable in (false,true),N in 1:3
        rho=[Inf;[100.,500.,20.][1:N]]
        epsilon=eps0*[1.;[10.,4.,30.][1:N]]
        mu=mu0*[1.;(permeable ? [1.,100.,3.] : ones(3))[1:N]]
        thickness=[Inf;[2.,7.][1:(N-1)];Inf]
        s=2π*frequency*im
        leaf=EI.Formula(:Nakagawa1973)(rho,epsilon,mu,s,nothing,nothing,thickness)
        negative=EI.Formula(:Nakagawa1973)(rho,epsilon,mu,-s,nothing,nothing,thickness)
        g=s.*mu.*(inv.(rho).+s.*epsilon)
        for pair in (E.EarthPair(1,1,(10.,10.),.02,(1,1)),
                E.EarthPair(1,2,(10.,15.),3.,(1,1)))
            self=pair.row==pair.column
            H=sum(pair.heights);x=self ? 0. : pair.separation
            d=self ? pair.separation : hypot(x,pair.heights[1]-pair.heights[2])
            mirror=hypot(x,H)
            correction=quadgk(0.,Inf;rtol=1e-10) do t
                u=t/H
                a=sqrt.(u^2 .+ g[2:end].-g[1])
                b=a./mu[2:end]
                if N==1
                    B=inv(u+mu0*b[1])
                elseif N==2
                    c1=b[1]+b[2]
                    c2=(b[1]-b[2])*exp(-2a[1]*thickness[2])
                    B=(c1+c2)/((u+mu0*b[1])*c1+(u-mu0*b[1])*c2)
                else
                    d1=thickness[2];d2=sum(thickness[2:3])
                    middle=exp(2a[2]*(d1-d2))
                    c1=(b[1]+b[2])*(b[2]+b[3])+
                        (b[1]-b[2])*(b[2]-b[3])*middle
                    c2=((b[1]-b[2])*(b[2]+b[3])+
                        (b[1]+b[2])*(b[2]-b[3])*middle)*exp(-2a[1]*d1)
                    B=(c1+c2)/((u+mu0*b[1])*c1+(u-mu0*b[1])*c2)
                end
                B*exp(-t)*cos(x*u)/H
            end[1]
            expected=s*mu0/(2π)*(log(mirror/d)+2correction)
            kind=self ? Val(:self) : Val(:mutual)
            value=leaf(kind,pair)
            @test value ≈ expected rtol=2e-7
            @test negative(kind,pair) ≈ conj(value) rtol=2e-7
            if N==1&&!permeable
                wise=EI.Formula(:Wise1934)(rho,epsilon,mu,s,nothing)
                @test wise(kind,pair) ≈ value rtol=2e-7
                # Independent difference-of-roots form, normalized by sqrt(contrast).
                contrast=g[2]-g[1]
                scale=abs(sqrt(contrast))
                delta=contrast/scale^2
                integral=quadgk(0.,Inf;rtol=1e-10) do v
                    (sqrt(v^2+delta)-v)/delta*
                        exp(-v*scale*H)*cos(v*scale*x)
                end[1]
                @test wise(kind,pair) ≈
                    s*mu0/(2π)*(log(mirror/d)+2integral) rtol=2e-7
            end
        end
    end
    for T in (Float32,Float64,BigFloat)
        rho=T[Inf,100,500,20];epsilon=T(eps0)*T[1,10,4,30]
        mu=T(mu0)*T[1,1,3,2];thickness=T[Inf,2,7,Inf]
        s=complex(zero(T),100T(π))
        pair=E.EarthPair(1,2,(T(10),T(15)),T(3),(1,1))
        for id in (:Nakagawa1973,:Wise1934)
            leaf=id==:Nakagawa1973 ?
                EI.Formula(id)(rho,epsilon,mu,s,nothing,nothing,thickness) :
                EI.Formula(id)(rho[1:2],epsilon[1:2],mu[1:2],s,nothing)
            value=leaf(Val(:mutual),pair)
            @test value isa Complex{T}
            @test isfinite(value)
        end
    end
    @test_throws DimensionMismatch EI.Formula(:Nakagawa1973)(
        [Inf,100.,200.,300.,400.],fill(eps0,5),fill(mu0,5),100π*im,
        nothing,nothing,[Inf,1.,2.,3.,Inf])
end
