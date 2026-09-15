@testitem "Engine / internal shunt / geometry and local operator" tags=[:unit] begin
    using LinearAlgebra
    E = LineCableModels.Engine
    include(joinpath(pkgdir(LineCableModels),"test","support","internal_shunt.jl"))
    design = internal_shunt_test_design()
    domain,bp = internal_shunt_test_domain(design)
    @test domain.terminals == 1:3
    @test length(domain.wires) == 6
    @test length(domain.tapes) == 1
    @test domain.material.eps_r == 1
    @test only(domain.left).material.eps_r == 3.0
    @test only(domain.right).material.eps_r == 3.0
    renamed = internal_shunt_test_design(suffix="_renamed")
    other,_ = internal_shunt_test_domain(renamed)
    @test E._shunt_domain_equal(domain,other)
    values = E._shunt_values(domain,Formulation().methods)
    @test values[1:3] ≈ [domain.a, domain.b, domain.material.eps_r]
    host = only(filter(region->region.primitive isa LineCableModels.DataModel.DifferenceShape &&
        region.primitive.outer isa Annulus,design.geometry.regions))
    incomplete = merge((primitive=host.primitive,source=host.source),
        (primitive=LineCableModels.DataModel.DifferenceShape(host.primitive.outer,
            Base.front(host.primitive.holes)),))
    @test E._shunt_host_domain(design,bp,incomplete,1,0,Float64) === nothing
    duplicated = merge(incomplete,(primitive=LineCableModels.DataModel.DifferenceShape(
        host.primitive.outer,(first(host.primitive.holes),Base.front(host.primitive.holes)...)),))
    @test E._shunt_host_domain(design,bp,duplicated,1,0,Float64) === nothing
    g = @inferred E._shunt_data(values,domain)
    level = (wire=32,order=16,quadrature=128,modes=256)
    result = @inferred E._shunt_capacitance(g;level,audit=true)
    production = E._shunt_capacitance(g;level)
    @test production.C == result.C
    @test production.diagnostic.boundary_residual === nothing
    @test result.C ≈ transpose(result.C) rtol=1e-5
    @test minimum(eigvals(Symmetric(result.C))) > 0
    @test sum(result.C[1,:]) > 0
    @test result.diagnostic.boundary_residual < 0.06
    @test result.state === nothing
    @test result.diagnostic.matrix_bytes <= E.INTERNAL_SHUNT_MATRIX_BYTES
    # The explicit owner matrix budget is retained. The inherited warmed
    # allocation number had no current performance requirement and is retired.
    doubled = merge(g,(epsilon=2g.epsilon,
        left=[merge(l,(epsilon=2l.epsilon,)) for l in g.left],
        right=[merge(l,(epsilon=2l.epsilon,)) for l in g.right],
        Rleft=g.Rleft/2,Rright=g.Rright/2,Rtotal=g.Rtotal/2))
    @test E._shunt_capacitance(doubled;level).C ≈ 2result.C rtol=1e-8
    rotation = 0.317
    rotated = merge(g,(wires=[merge(w,(x=w.x*cos(rotation)-w.y*sin(rotation),
        y=w.x*sin(rotation)+w.y*cos(rotation))) for w in g.wires],
        tapes=[merge(t,(phi=t.phi+rotation,)) for t in g.tapes]))
    @test E._shunt_capacitance(rotated;level).C ≈ result.C rtol=2e-4
    @test_throws BoundarySolveError E._shunt_capacitance(g;
        level=(wire=100000,order=16,quadrature=128,modes=1))
    formulation=Formulation(shunt_model=formula(:boundary;options=(resolution=level,)))
    blueprints = only(E.flatten(LineCableModelsCoaxial(), [design,renamed], Float64, [formulation]))
    @test E.LocalCableData(blueprints).shunt_details.solves == 1
    @test blueprints[1].shunt[1].C === blueprints[2].shunt[1].C
    @test isempty(E.flatten(LineCableModelsCoaxial(), design,
        Formulation(insulation_admittance=:lossy)).shunt)
    cable = E.LocalCableData(first(blueprints))
    p,y = zeros(ComplexF64,3,3),zeros(ComplexF64,3,3)
    layers = zeros(ComplexF64,length(bp.dielectrics))
    coefficients,tails = zeros(ComplexF64,3),zeros(ComplexF64,3)
    local_shunt = cable.shunt
    E.cable_potential!(p,cable,Formulation().methods,50.0,20.0,100pi*im,
        layers,coefficients,tails)
    E.cable_admittance!(y,cable,Formulation().methods,50.0,20.0,100pi*im,layers)
    @test y ≈ 100pi*im*inv(p) rtol=1e-10
    @test @allocated(E._shunt_potential!(p,local_shunt)) == 0
    @test @allocated(E._shunt_admittance!(y,local_shunt,100pi*im)) == 0
end

@testitem "Engine / internal shunt / nested shields and reduction" tags=[:unit] begin
    using LinearAlgebra
    E = LineCableModels.Engine
    # Two local domains share a closed shield (terminal 3). The outer domain's
    # inner-anchor charge includes all terminals enclosed by that shield.
    a = (Diagonal([2.0,3.0])+[1.0,-1.0]*[1.0 -1.0]).*1e-9
    b = (Diagonal([5.0,7.0])+2*[1.0,-1.0]*[1.0 -1.0]).*1e-9
    blocks = [E.InternalShuntBlock(1:5,1:3,a,inv(a)),
        E.InternalShuntBlock(1:5,3:5,b,inv(b))]
    # Outermost jacket to the external reference remains a radial contribution.
    p = fill(inv(7e-9),5,5)
    y = zeros(ComplexF64,5,5)
    s = 100pi*im
    y[5,5] = s*7e-9
    E._shunt_potential!(p,blocks)
    E._shunt_admittance!(y,blocks,s)
    @test y ≈ s*inv(p) rtol=1e-12
    # Ordinary terminal bundling sums charges and imposes equal potentials.
    incidence = [1 0 0; 0 1 0; 0 1 0; 0 0 1; 0 0 1]
    bundled = transpose(incidence)*y*incidence
    @test bundled ≈ transpose(bundled)
    @test minimum(eigvals(Symmetric(imag.(bundled)))) > 0
end

@testitem "Engine / internal shunt / multiple independent open terminals" tags=[:unit] begin
    using LinearAlgebra
    E = LineCableModels.Engine
    include(joinpath(pkgdir(LineCableModels),"test","support","internal_shunt.jl"))
    original,_ = internal_shunt_test_domain(internal_shunt_test_design(tapes=false))
    # Independent test of the numerical terminal map, not permission to bypass
    # the existing public coaxial lowering rules for radial conductor ordering.
    wires = [merge(wire,(terminal=iseven(index) ? 2 : 3,))
        for (index,wire) in pairs(original.wires)]
    domain = E.InternalShuntDomain(original.design,1:4,1:4,original.a,original.b,
        original.material,original.left,original.right,wires,original.tapes)
    @test domain.terminals == 1:4
    @test sort(unique(getproperty.(domain.wires,:terminal))) == [2,3]
    g = E._shunt_data(E._shunt_values(domain,Formulation().methods),domain)
    result = E._shunt_capacitance(g;level=(wire=32,order=16,quadrature=128,modes=256))
    @test size(result.C) == (3,3)
    @test result.C ≈ transpose(result.C) rtol=1e-6
    @test result.C[2,3] < 0
    @test minimum(eigvals(Symmetric(result.C))) > 0
end
