@testitem "DataModel / implicit power-cell and strand derivatives" tags=[:unit, :extension] begin
    using Measurements
    const DM = LineCableModels.DataModel
    derivative(v, z) = v isa Measurement ? Measurements.derivative(v, z) : 0.0
    function geometry(s, q, shift, angle, dx, count)
        boundary = [(-s,-s),(1.3s,-s),(1.3s,0.8s),(-s,0.8s)]
        sites = [(-0.55s+shift,-0.2s),(0.65s,-0.1s),(0s,0.5s),(0s,-0.65s)][1:count]
        a = DM.signed_polygon_area(boundary)
        targets = count == 2 ? [a*q,a*(1-q)] : count == 3 ?
            [a*q,0.2a,a*(0.8-q)] : [a*q,0.2a,0.15a,a*(0.65-q)]
        transform(p) = (dx+cos(angle)*p[1]-sin(angle)*p[2],sin(angle)*p[1]+cos(angle)*p[2])
        return transform.(boundary), transform.(sites), targets
    end
    for count in (2,3,4), axis in 1:5
        base = [1.0,0.4,0.0,0.2,0.1]
        z = measurement(base[axis],0.01)
        values = Any[base...]; values[axis] = z
        args = geometry(values...,count)
        cells, weights = DM.balance_power_cells(args...)
        @test nominal(weights[end]) == 0
        @test derivative(weights[end],z) == 0
        for h in (1e-5,3e-6)
            hi, lo = copy(base), copy(base)
            hi[axis] += h; lo[axis] -= h
            chi = first(DM.balance_power_cells(geometry(hi...,count)...))
            clo = first(DM.balance_power_cells(geometry(lo...,count)...))
            for (i,cell) in enumerate(cells)
                a = DM.signed_polygon_area(cell)
                @test nominal(a) ≈ nominal(args[3][i]) rtol=1e-9
                @test derivative(a,z) ≈ derivative(args[3][i],z) atol=1e-8
                for (v,p,m) in zip(DM.polygon_centroid(cell),
                        DM.polygon_centroid(chi[i]),DM.polygon_centroid(clo[i]))
                    @test derivative(v,z) ≈ (p-m)/(2h) atol=1e-7 rtol=1e-5
                end
            end
        end
    end
    function strand_geometry(s,width,a,angle,dx)
        ca, sa = cos(angle), sin(angle)
        cell = [(dx+ca*x-sa*y,sa*x+ca*y) for (x,y) in
            [(-0.2s,-0.5s),(width*s,-0.5s),(width*s,0.5s),(-0.2s,0.5s)]]
        return cell, a*s^2
    end
    for axis in 1:5, occupancy in (0.4,0.6)
        base = [1.0,0.4,occupancy,0.2,0.1]
        z = measurement(base[axis],0.001)
        values = Any[base...]; values[axis] = z
        cell, a = strand_geometry(values...)
        # Full occupancy is tested only along support-preserving directions.
        occupancy == 0.6 && axis in (2,3) && continue
        points = DM.area_preserving_strand(cell,a;angle=values[4])
        for p in points, i in eachindex(cell)
            left, right = cell[i],cell[mod1(i+1,length(cell))]
            @test nominal((right[1]-left[1])*(p[2]-left[2])-
                (right[2]-left[2])*(p[1]-left[1])) >= -1e-12
        end
        @test nominal(DM.signed_polygon_area(points)) ≈ nominal(a) atol=1e-10
        @test derivative(DM.signed_polygon_area(points),z) ≈ derivative(a,z) atol=1e-7
        for h in (1e-5,3e-6)
            hi, lo = copy(base), copy(base)
            hi[axis] += h; lo[axis] -= h
            phi = DM.area_preserving_strand(strand_geometry(hi...)...;angle=hi[4])
            plo = DM.area_preserving_strand(strand_geometry(lo...)...;angle=lo[4])
            for (v,p,m) in zip(DM.polygon_centroid(points),
                    DM.polygon_centroid(phi),DM.polygon_centroid(plo))
                @test derivative(v,z) ≈ (p-m)/(2h) atol=1e-7 rtol=1e-5
            end
        end
    end
    z = measurement(1.0,0.1)
    cell = [(-z,-z),(z,-z),(z,z),(-z,z)]
    cells, weights = DM.balance_power_cells(cell,[(0.,0.)],[4z^2])
    @test only(cells) == cell
    @test only(weights) == 0
    # Inventory transitions are discrete; no continuous derivative claims are
    # made across the boundary where the next complete course appears.
    @test DM.course_count(18-1e-5,1.) == 1
    @test DM.course_count(18+1e-5,1.) == 2
end

@testitem "DataModel / bounded construction retains scale pose and fillet derivatives" tags=[:unit, :extension] begin
    using Measurements
    const DM = LineCableModels.DataModel
    derivative(v,z) = v isa Measurement ? Measurements.derivative(v,z) : 0.0
    copper = Material(kind=:conductor,rho=1.72e-8)
    function design(kind,s,rotation,dx,fillet)
        part = if kind === :sector
            boundary = resolve(EmptyBoundary(),Sector(span=pi/3,
                r_base=fillet*0.006s,r_back=0.006s,fillet=fillet*0.006s))
            _, polygons, _ = DM.sector_courses(boundary,Disk(0.0007s))
            return [resolve(Pose2(dx,0.,rotation),p) for p in polygons]
        elseif kind === :rectangle
            stranded(copper;center=Disk(0.0002s),
                shape=Rectangle(0.0003s,0.0001s),boundary=Disk(0.0006s))
        else
            stranded(copper;shape=Disk(0.0002s),boundary=Disk(0.0009s),compact=true)
        end
        cable = build(CableDesign,"derivative",assembly(at(terminal(:core,part),
            Pose2(dx,0.,rotation))))
        return [r.primitive for r in cable.geometry.regions]
    end
    for kind in (:circle,:sector,:rectangle), axis in 1:4
        kind !== :sector && axis == 4 && continue
        base = [1.,0.2,0.001,0.04]
        z = measurement(base[axis],0.001)
        values = Any[base...]; values[axis] = z
        shapes = design(kind,values...)
        plain = design(kind,base...)
        @test length(shapes) == length(plain)
        for (p,p0) in zip(shapes,plain)
            @test nominal(area(p)) ≈ area(p0) rtol=1e-9
            @test all(isapprox.(nominal.(centroid(p)),centroid(p0);atol=1e-12))
        end
        for h in (1e-5,3e-6)
            hi, lo = copy(base), copy(base)
            hi[axis] += h; lo[axis] -= h
            plus, minus = design(kind,hi...), design(kind,lo...)
            @test length(plus) == length(minus) == length(shapes)
            for (p,pplus,pminus) in zip(shapes,plus,minus)
                @test derivative(area(p),z) ≈ (area(pplus)-area(pminus))/(2h) atol=1e-12 rtol=1e-5
                for (v,vp,vm) in zip(centroid(p),centroid(pplus),centroid(pminus))
                    @test derivative(v,z) ≈ (vp-vm)/(2h) atol=2e-8 rtol=1e-4
                end
            end
        end
    end
end
