@testitem "Engine / literature assimilation / Ametani snow potential equivalence" begin
    using LineCableModels, QuadGK
    E=LineCableModels.Engine; EA=E.EarthAdmittance
    μ0=4π*1e-7; ε0=8.8541878128e-12
    # Japanese coefficients (7), with the denominator product explicitly
    # printed in the English A2 witness. Multiply c3,c4 and their dependent
    # expressions by c0 before evaluating to avoid positive exponentials.
    function snow_kernel(λ,γ,μ,d)
        g0,g1,g2=γ; m0,m1,m2=μ
        a1=sqrt(λ^2+g1-g0); a2=sqrt(λ^2+g2-g0)
        b1=a1/m1; b2=a2/m2
        t1=g0/g1; t2=g1/g2
        e=exp(-2a1*d)
        c1=b1+b2; c2=(b1-b2)*e
        c5=a1*t1; c6=m2*a2*t2/m1
        C3=c6+a1; C4=(c6-a1)*e
        A1=(c1+c2)/((λ+m0*b1)*c1+(λ-m0*b1)*c2)
        A3=2m1*C4*c5+(C3-C4)*(m1*c5+m0*λ)
        A4=(λ/m0+b1)*c1+(λ/m0-b1)*c2
        A2=(4b1*e*c5*(1-t2)+(c1+c2)*(C3-C4)*(1-t1))/(A3*A4)
        return A1-λ*A2
    end
    for frequency in (10.0,1000.0,1e5,1e6), thickness in (0.0,2.0,10.0),
        mur in (1.0,3.0)
        s=complex(0.0,2π*frequency)
        rho=[Inf,1e7,100.0]; ε=ε0*[1.0,4.0,10.0]; μ=μ0*[1.0,1.0,mur]
        gamma=[s*m*(inv(r)+s*e) for (r,m,e) in zip(rho,μ,ε)]
        f=EA.Formula(:Papadopoulos2009)(rho,ε,μ,s,nothing,nothing,[Inf,thickness,Inf])
        for (hi,hj,x) in ((10.0,15.0,3.0),(2.0,4.0,0.0))
            pair=E.EarthPair(1,2,(hi,hj),x,(1,1))
            integral=quadgk(λ->snow_kernel(λ,gamma,μ,thickness)*
                exp(-(hi+hj)*λ)*cos(x*λ),0.0,Inf;rtol=1e-10)[1]
            source=(log(hypot(x,hi+hj)/hypot(x,hi-hj))+2integral)/(2π*ε0)
            @test f(Val(:mutual),pair) ≈ source rtol=3e-8
        end
    end
end
