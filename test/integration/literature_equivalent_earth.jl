@testitem "EarthProps / literature assimilation / equivalent earth source recursions" begin
    using LineCableModels,QuadGK
    E=LineCableModels.Engine;EP=LineCableModels.EarthProps;EH=EP.EHEM
    for T in (Float32,Float64,BigFloat),N in 1:6,frequency in (T(1),T(50),T(1e6))
        mu0=T(4)*T(π)/T(10)^7;eps0=T(88541878128)/T(10)^22
        rho=T[100,500,20,200,30,1000][1:N]
        epsilon=T[10,4,30,7,3,15][1:N]
        depths=T[2,7,1,3,4][1:N-1]
        s=complex(zero(T),2T(π)*frequency)
        pair=E.EarthPair(1,2,(T(10),T(15)),T(3),(1,1))
        for id in (:MartinsBritto2020,:Xue2021)
            magnetic=id==:MartinsBritto2020 ? T[1,3,2,1,4,1][1:N] : ones(T,N)
            layers=ntuple(N) do k
                EP.EarthLayer(rho[k],epsilon[k],magnetic[k],k==N ? T(Inf) : depths[k])
            end
            model=build(EP.EarthModel,layers)
            r=getfield.(model.layers,:rho);ep=getfield.(model.layers,:eps_r)
            mu=getfield.(model.layers,:mu_r)
            value=EH.Formula(id)(Val(:overhead),collect(r),collect(ep),collect(mu),
                model,pair,frequency)
            @test value isa EP.EarthMaterial{T}
            if id==:MartinsBritto2020
                expected=inv(rho[end])
                # Algebraically independent tanh form of (20)--(22).
                for k in N-1:-1:1
                    q=sqrt(inv(rho[k]));b=sqrt(expected)
                    t=tanh(depths[k]*sqrt(T(π)*frequency*mu0*magnetic[k]/rho[k]))
                    expected=(q*(b+q*t)/(q+b*t))^2
                end
                @test inv(value.rho) ≈ expected rtol=200eps(T)
                @test value.eps_r==epsilon[end]
                @test value.mu_r==magnetic[end]
            else
                air=s^2*mu0*eps0
                q=sqrt.(s*mu0.*(inv.(rho).+s*eps0.*epsilon).-air)
                expected=q[end]
                for k in N-1:-1:1
                    t=tanh(depths[k]*q[k])
                    expected=q[k]*(expected+q[k]*t)/(q[k]+expected*t)
                end
                recovered=s*mu0*(inv(value.rho)+s*eps0*value.eps_r)-air
                @test recovered ≈ expected^2 rtol=200eps(T)
            end
        end
    end
    # Published manuscript Table VI, model 1. Input conductivities and
    # thickness have two decimal places; output has four. Compare within
    # the interval implied by the printed input precision.
    model=build(EP.EarthModel,(EP.EarthLayer(1/.00268,1.,1.,2.69),
        EP.EarthLayer(1/.00688,1.,1.)))
    values=map(p->collect(getfield.(model.layers,p)),(:rho,:eps_r,:mu_r))
    pair=E.EarthPair(1,1,(10.,10.),.005,(1,1))
    result=EH.Formula(:MartinsBritto2020)(Val(:overhead),values...,model,pair,50.)
    bounds=[1000EH.equivalent_conductivity_step(Val(:MartinsBritto2020),
        (2.68+ds1)/1000,(6.88+ds2)/1000,2.69+dh,50.,4π*1e-7)
        for ds1 in (-.005,.005),ds2 in (-.005,.005),dh in (-.005,.005)]
    @test minimum(bounds)<=6.8580<=maximum(bounds)
    @test minimum(bounds)<=1000/result.rho<=maximum(bounds)
end
