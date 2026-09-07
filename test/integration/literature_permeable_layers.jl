@testitem "Engine / literature assimilation / permeable layered source equations" begin
    using LineCableModels, QuadGK
    E=LineCableModels.Engine; EI=E.EarthImpedance
    μ0=4π*1e-7; ε0=8.8541878128e-12
    function literal_A(q,μ,depths)
        q2,q3,q4=q; m2,m3,m4=μ; t1,t2=depths
        upper=q2/sinh(t1*q2)
        lower=((m3/m4)*q4*cosh(t2*q3)+q3*sinh(t2*q3))/
            (q3*cosh(t2*q3)+(m3/m4)*q4*sinh(t2*q3))
        return upper*(cosh(t1*q2)-q2/
            (q2*cosh(t1*q2)+(m2/m3)*q3*sinh(t1*q2)*lower))
    end
    for frequency in (1.0,50.0,1e4,1e6), mur in ([1.0,1.0,1.0],[1.0,2.0,7.0])
        s=complex(0.0,2π*frequency); soilμ=μ0*mur
        rho=[100.0,500.0,50.0]; epsilon=ε0*[10.0,6.0,20.0]; depths=[2.0,3.0]
        f=EI.Formula(:Moghram1998)([Inf;rho],[ε0;epsilon],[μ0;soilμ],s,
            nothing,nothing,[Inf;depths;Inf])
        for pair in (E.EarthPair(1,1,(10.0,10.0),0.02,(1,1)),
                E.EarthPair(1,2,(10.0,15.0),3.0,(1,1)),
                E.EarthPair(1,2,(10.0,15.0),0.0,(1,1)))
            self=pair.row==pair.column; H=sum(pair.heights)
            x=self ? 0.0 : pair.separation
            ideal=self ? log(H/pair.separation) :
                log(hypot(H,x)/hypot(pair.heights[1]-pair.heights[2],x))
            correction=quadgk(0.0,5.0;rtol=1e-10) do λ
                q=sqrt.([λ^2+s*m/r for (m,r) in zip(soilμ,rho)])
                A=literal_A(q,soilμ,depths)
                exp(-H*λ)*cos(x*λ)/(λ+(μ0/soilμ[1])*A)
            end[1]
            expected=s*μ0/(2π)*(ideal+2correction)
            kind=self ? Val(:self) : Val(:mutual)
            @test f(kind,pair) ≈ expected rtol=3e-8
        end
        for displacement in (false,true)
            args=([Inf;rho[1:2]],[ε0;epsilon[1:2]],[μ0;soilμ[1:2]],s,
                nothing,nothing,[Inf,2.0,Inf])
            w=EI.Formula(:Wedepohl1966;displacement_current=displacement)(args...)
            pair=E.EarthPair(1,2,(10.0,15.0),3.0,(1,1))
            correction=quadgk(0.0,5.0;rtol=1e-10) do λ
                q=sqrt.([λ^2+s*m/r+(displacement ? s^2*(m*e-μ0*ε0) : 0) for
                    (m,r,e) in zip(soilμ[1:2],rho[1:2],epsilon[1:2])])
                q2,q3=q; m2,m3=soilμ[1:2]; d=2.0
                A=q2/sinh(d*q2)*(cosh(d*q2)-
                    m3*q2/(m3*q2*cosh(d*q2)+m2*q3*sinh(d*q2)))
                exp(-25λ)*cos(3λ)/(λ+(μ0/m2)*A)
            end[1]
            source=s*μ0/(2π)*(log(hypot(25.0,3.0)/hypot(5.0,3.0))+2correction)
            @test w(Val(:mutual),pair) ≈ source rtol=3e-8
            if displacement
                modern=EI.Formula(:Papadopoulos2009)(args...)
                @test w(Val(:mutual),pair) ≈ modern(Val(:mutual),pair) rtol=3e-8
            end
        end
    end
    for frequency in (1.0,50.0,1e4), mur in (1.0,5.0)
        s=complex(0.0,2π*frequency)
        f=EI.Formula(:Moghram1998)([Inf,100.0,100.0,100.0],ε0*[1.0,10.0,10.0,10.0],
            μ0*[1.0,mur,mur,mur],s,nothing,nothing,[Inf,2.0,3.0,Inf])
        homogeneous=EI.Formula(:Wise1931)([Inf,100.0],ε0*[1.0,10.0],μ0*[1.0,mur],s,nothing)
        pair=E.EarthPair(1,2,(10.0,15.0),3.0,(1,1))
        @test f(Val(:mutual),pair) ≈ homogeneous(Val(:mutual),pair) rtol=3e-8
    end
end
