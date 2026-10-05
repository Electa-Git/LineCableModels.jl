@testitem "FEM / earth sizing length, ceiling and domain option" tags=[:extension] setup=[NativeFEMFixtures] begin
    using Gmsh
    N=NativeFEMFixtures
    expressions=["length"=>"FEMEarthBaseLength", "ceiling"=>"FEMEarthSizingCeilingActive",
        "size"=>"FEMEarthSizingLength", "domain"=>"DomainHalfwidth", "bulk"=>"MeshBulk"]
    for (rho,f,eps_r) in ((100.,50.,1.),(1000.,1e8,12.),(Inf,50.,1.),(Inf,1e8,1.))
        p=N.problem(;rho,frequency=f,eps_r)
        N.parameters(p;expressions) do a,b
            @test isequal(a,b)
            w=2pi*f;mu=4pi*1e-7;eps=8.8541878128e-12*eps_r
            q=sqrt(complex(-w*w*mu*eps,w*mu/rho))
            expected=min(real(q)>0 ? inv(real(q)) : Inf,2pi/abs(q))
            cap=sqrt(2e5/(w*mu))
            @test a["length"]≈expected rtol=2e-14
            @test a["ceiling"]==Float64(expected>cap)
            @test a["size"]≈min(expected,cap) rtol=2e-14
            @test a["domain"]≈max(5.,2min(expected,cap)) rtol=2e-14
        end
    end
    E=N.FEM
    @test E.computation_options(E.LineCableModelsFEM,ComputationOptions(domain_size_factor=3.)).data.domain_size_factor==3.
    err=try E.computation_options(E.LineCableModelsFEM,ComputationOptions(domain_skin_depths=2.));catch e;e;end
    @test err isa ArgumentError
    @test sprint(showerror,err)=="ArgumentError: unknown LineCableModelsFEM computation options: (:domain_skin_depths,)"
    N.bundle(N.problem()) do directory,entry
        @test occursin("DomainSizeFactor",read(joinpath(directory,"model_data.pro"),String))
        @test !occursin("DomainSkinDepths",read(joinpath(directory,"model_data.pro"),String))
    end
end

@testitem "FEM / long buried measurement lines retain remote grading" tags=[:extension] setup=[NativeFEMFixtures] begin
    N=NativeFEMFixtures
    N.geometry(N.problem(;frequency=.1,rho=100.,eps_r=1.,positions=[(0.,-.1)]);
            mesh=true) do g,log
        @test only(g.parser.get_number("DomainHalfwidth"))>1e4
        lines=g.model.get_entities_for_physical_group(1,
            Int(only(g.parser.get_number("MEASUREMENT_LINE"))))
        @test !isempty(lines)
        counts=map(lines) do curve
            _,tags,_=g.model.mesh.get_elements(1,curve)
            sum(length,tags)
        end
        # The local conductor target must not seed the whole distant line.
        @test sum(counts)<10_000
        _,_,blocks=g.model.mesh.get_elements(2)
        nodes=Set(vcat(blocks...))
        for curve in lines
            tags,_,_=g.model.mesh.get_nodes(1,curve,true)
            @test isempty(intersect(nodes,Set(tags)))
        end
    end
end

@testitem "FEM / loss-aware PML intervals and native parser parity" tags=[:extension] setup=[NativeFEMFixtures] begin
    N=NativeFEMFixtures
    expressions=["a_air"=>"FEMRootA~{0}","b_air"=>"FEMRootB~{0}",
        "a_earth"=>"FEMRootA~{1}","b_earth"=>"FEMRootB~{1}","target"=>"FEMPmlTarget",
        "ppw_air"=>"FEMPmlPPW~{0}","ppw_earth"=>"FEMPmlPPW~{1}"]
    for (label,d) in (("side","Side"),("top","Top"),("bottom","Bottom"))
        append!(expressions,[label*"_eta"=>"Pml$(d)Eta",label*"_A"=>"Pml$(d)Strength",
            label*"_L"=>"Pml$(d)Thickness",label*"_N"=>"FEMPml$(d)Layers"])
    end
    for f in (50.,1e8,3e8)
        N.parameters(N.problem(;frequency=f);expressions) do a,b
            @test isequal(a,b)
            for m in ("air","earth")
                q=complex(a["a_"*m],a["b_"*m]);rate=real(q)
                factor=rate>0 ? clamp(sqrt(.1abs(q)/rate),1.,3.) : 3.
                @test a["ppw_"*m]≈10factor rtol=2e-14
            end
            for (d,media) in (("side",("air","earth")),("top",("air",)),("bottom",("earth",)))
                X=a[d*"_L"]*(1+a[d*"_A"]/4-im*a[d*"_eta"]*a[d*"_A"]/4)
                phases=map(media) do m
                    q=complex(a["a_"*m],a["b_"*m]);E=real(q*X)
                    a["ppw_"*m]*abs(q)*abs(X)*(E>0 ? min(1,a["target"]/E) : 1)
                end
                @test a[d*"_N"]==max(48,ceil(maximum(phases)/(2pi)))
            end
        end
    end
end

@testitem "FEM / earth interface layer stays clear of buried cable" tags=[:extension] setup=[NativeFEMFixtures] begin
    N=NativeFEMFixtures
    ex=["thickness"=>"FEMEarthLayerThickness","active"=>"FEMEarthLayerActive",
        "clipped"=>"FEMEarthLayerClippedOrOmitted","decay"=>"MeshDecayEarth",
        "wave"=>"MeshWaveEarth","remote"=>"MeshRemoteEarth","cable"=>"CableSize~{0}"]
    for y in (1.,-1.,-.1)
        N.parameters(N.problem(;frequency=1e8,rho=1000.,eps_r=12.,positions=[(0.,y)]);expressions=ex) do a,b
            @test isequal(a,b)
            candidate=y<0 ? min(a["decay"],.5*(abs(y)-.01-a["cable"])) : a["decay"]
            active=a["wave"]<a["remote"] && candidate>=2a["wave"]
            @test a["active"]==active
            @test a["thickness"]≈(active ? candidate : 0.)
            @test a["clipped"]==((a["wave"]<a["remote"]) && a["thickness"]<a["decay"])
            y<0 && @test a["thickness"]<=.5*(abs(y)-.01-a["cable"])
        end
    end
end

@testitem "FEM / box edge grading and independent measurement mesh" tags=[:extension] setup=[NativeFEMFixtures] begin
    using Gmsh
    N=NativeFEMFixtures
    problem=N.problem(;frequency=1e6,rho=.1,eps_r=1.,radius=.01)
    N.geometry(problem;options=(mesh_size_factor=2.5,),mesh=true) do g,log
        D=only(g.parser.get_number("DomainHalfwidth"));x=only(g.parser.get_number("Xcenter"))+D
        wave=only(g.parser.get_number("MeshWaveEarth"));bulk=only(g.parser.get_number("MeshBulk"))
        last=only(g.parser.get_number("MeshRemoteEarth"));first=min(wave,bulk,last)
        curves=filter(g.model.get_entities(1)) do entity
            box=g.model.get_bounding_box(entity...)
            abs(box[1]-x)<1e-6 && abs(box[4]-x)<1e-6 && abs(box[2]+D)<1e-6 && abs(box[5])<1e-6
        end
        @test length(curves)==1
        _,xyz,_=g.model.mesh.get_nodes(1,only(curves)[2],true)
        widths=diff(sort(reshape(xyz,3,:)[2,:]))
        @test first>0
        @test widths[end]<=first*(1+1e-8)
        @test widths[end]<widths[1]
        # Floating lines never share a node with a two-dimensional element.
        _,_,blocks=g.model.mesh.get_elements(2);nodes=Set(vcat(blocks...))
        lines=g.model.get_entities_for_physical_group(1,Int(only(g.parser.get_number("MEASUREMENT_LINE"))))
        @test !isempty(lines)
        form=Formulation(:LineCableModelsFEM;options=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
        model=N.FEM._resolved_fem_model(problem,form)
        ratios=N.FEM._inspect_loaded_mesh(model,"native measurement mesh")
        # Allow twice the 0.25 target because Delaunay edges can be smaller than their targets.
        @test all(r -> isfinite(r) && 0<r<=.5,ratios)
        for curve in lines
            tags,_,_=g.model.mesh.get_nodes(1,curve,true)
            @test isempty(intersect(nodes,Set(tags)))
        end
    end
end

@testitem "FEM / interface seeds merge near-duplicate cable abscissae" tags=[:extension] setup=[NativeFEMFixtures] begin
    N=NativeFEMFixtures
    p=N.problem(;positions=[(0.,1.),(4e-15,1.5)],frequency=1e6,rho=100.)
    N.geometry(p) do g,log
        x=g.parser.get_number("FEMInterfaceX")
        tol=only(g.parser.get_number("FEMInterfaceTolerance"))
        @test count(v->abs(v)<tol,x)==1
        @test minimum(diff(x))>=tol
        @test !any(l->occursin("closer than the geometrical tolerance",l),log)
        for curve in Int.(g.parser.get_number("FEMFiniteInterface"))
            endpoints=g.model.get_boundary([(1,curve)],false,false,false)
            xyz=[g.model.get_value(dim,tag,Float64[]) for (dim,tag) in endpoints]
            @test length(xyz)==2
            @test abs(xyz[2][1]-xyz[1][1])>=tol
        end
    end
end
