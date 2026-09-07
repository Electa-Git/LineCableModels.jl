@testitem "Engine / literature assimilation / Tsiamitros independent layer boundary solve" begin
    using LineCableModels,LinearAlgebra,QuadGK
    E=LineCableModels.Engine;EI=E.EarthImpedance
    mu0=4π*1e-7;eps0=8.8541878128e-12
    root(z,s)=iszero(imag(z))&&real(z)<0 ?
        complex(0.,sign(imag(s))*sqrt(-real(z))) : sqrt(z)
    # Exact dynamic stiffness of each homogeneous interval. Source and
    # receiver are added as nodes; no spatial discretization is involved.
    function boundary_green(lambda,s,rho,epsilon,mu,thickness,zs,zt)
        boundaries=[0.;cumsum(thickness[2:end-1])]
        nodes=sort!(unique([boundaries;zs;zt]))
        n=length(nodes);matrix=zeros(ComplexF64,n,n)
        g=s.*mu.*(inv.(rho).+s.*epsilon)
        a=root.(lambda^2 .+ g,Ref(s))
        matrix[1,1]+=a[1]/mu[1]
        matrix[end,end]+=a[end]/mu[end]
        for k in 1:n-1
            mid=(nodes[k]+nodes[k+1])/2
            medium=searchsortedlast(boundaries,mid)+1
            q=a[medium];v=exp(-q*(nodes[k+1]-nodes[k]))
            diag=q/mu[medium]*(1+v^2)/(1-v^2)
            off=-2q/mu[medium]*v/(1-v^2)
            matrix[k,k]+=diag;matrix[k+1,k+1]+=diag
            matrix[k,k+1]+=off;matrix[k+1,k]+=off
        end
        rhs=zeros(ComplexF64,n);rhs[findfirst(==(zs),nodes)]=1
        return (matrix\rhs)[findfirst(==(zt),nodes)]
    end
    for N in 1:4,frequency in (50.,1e5)
        s=2π*frequency*im
        rho=[Inf;[100.,500.,20.,200.][1:N]]
        epsilon=eps0*[1.;[10.,4.,30.,7.][1:N]]
        mu=mu0*[1.;[1.,3.,.7,2.][1:N]]
        thickness=[Inf;[2.,3.,4.][1:N-1];Inf]
        tops=[0.;cumsum(thickness[2:end-1])]
        leaf=EI.Formula(:Tsiamitros2008)(rho,epsilon,mu,s,nothing,nothing,thickness)
        for m in 1:N,l in 1:m,lambda in (.01,.1,1.)
            zm=tops[m]+.8;zl=tops[l]+.3
            n=N+1;a=zeros(ComplexF64,n);num=similar(a);den=similar(a)
            EI.downward_coefficients!(Val(:Tsiamitros2008),num,den,a,leaf.state,lambda)
            value=mu[m+1]/2*EI.spectral_kernel(Val(:Tsiamitros2008),
                leaf.state,a,num,den,l,m,.8,.3)
            expected=boundary_green(lambda,s,rho,epsilon,mu,thickness,zm,zl)
            @test value ≈ expected rtol=5e-11
        end
    end
    rho=[Inf,100.,500.,20.];epsilon=eps0*[1.,10.,4.,30.]
    mu=mu0*[1.,1.,3.,.7];thickness=[Inf,2.,3.,Inf]
    s=100π*im
    leaf=EI.Formula(:Tsiamitros2008)(rho,epsilon,mu,s,nothing,nothing,thickness)
    negative=EI.Formula(:Tsiamitros2008)(rho,epsilon,mu,-s,nothing,nothing,thickness)
    for (zs,zt,ls,lt) in ((-10.,-12.,1,1),(-10.,.5,1,2),
            (-10.,5.5,1,4),(.5,1.5,2,2),(.5,5.5,2,4),(2.5,3.5,3,3))
        pair=E.EarthPair(1,2,(-zs,-zt),.7,(ls,lt))
        reversed=E.EarthPair(2,1,(-zt,-zs),.7,(lt,ls))
        # Break the air radiation cut explicitly in the independent integral.
        cut=abs(s)*sqrt(mu[1]*epsilon[1])
        integral=quadgk(0.,cut,Inf;rtol=2e-8) do lambda
            boundary_green(lambda,s,rho,epsilon,mu,thickness,zs,zt)*cos(.7lambda)
        end[1]
        value=leaf(Val(:mutual),pair)
        @test value ≈ s/π*integral rtol=5e-6
        @test leaf(Val(:mutual),reversed) ≈ value rtol=3e-9
        @test negative(Val(:mutual),pair) ≈ conj(value) rtol=3e-9
    end
end
