# Integrate a retained first-order edge field from its three vertex values.
# This controls extraction of the solved field, not PDE assembly or physics.
module EdgeFieldControl
using LinearAlgebra
using Gmsh

function read_field(path)
    before=Set(gmsh.view.get_tags())
    gmsh.merge(path)
    views=setdiff(Set(gmsh.view.get_tags()),before)
    length(views)==1 || error("one retained edge-field view required")
    view=only(views)
    try
        kinds,counts,blocks=gmsh.view.get_list_data(view)
        fields=NamedTuple[]
        for (kind,count,block) in zip(kinds,counts,blocks)
            kind=="VT" || error("edge control requires vector triangles, received $kind")
            width=length(block)÷count
            width==27 || error("expected two complex field components in two time steps, width=$width")
            for element in 0:count-1
                row=block[element*width+1:(element+1)*width]
                vertices=[(row[i],row[i+3]) for i in 1:3]
                values=complex.(row[10:18],row[19:27])
                # bt=(ax-b*y, ay+b*x) on each first-order triangle.
                system=zeros(6,3);rhs=zeros(ComplexF64,6)
                for i in 1:3
                    x,y=vertices[i]
                    system[2i-1,:]=[1.,0.,-y]
                    system[2i,:]=[0.,1.,x]
                    rhs[2i-1]=values[3i-2];rhs[2i]=values[3i-1]
                end
                coefficients=system\rhs
                residual=maximum(abs,system*coefficients-rhs)
                scale=opnorm(system,Inf)*norm(coefficients,Inf)+norm(rhs,Inf)
                bound=residual+128eps(Float64)*scale
                residual<=1e-10scale || error("retained field is not a first-order edge field")
                push!(fields,(;vertices,coefficients,bound))
            end
        end
        return fields
    finally
        gmsh.view.remove(view)
    end
end

function integrate(fields,vertices)
    total=0.0im;uncertainty=0.0
    for (p,q) in zip(vertices[1:end-1],vertices[2:end]),field in fields
        triangle=field.vertices
        any(k->maximum(v[k] for v in triangle)<min(p[k],q[k]) ||
            minimum(v[k] for v in triangle)>max(p[k],q[k]),1:2) && continue
        origin=collect(triangle[1])
        transform=hcat(collect(triangle[2]).-origin,collect(triangle[3]).-origin)
        uv=transform\(collect(p).-origin)
        delta=transform\(collect(q).-collect(p))
        # Barycentric coordinates are affine along a straight segment. Solve
        # their three nonnegativity constraints; integrate the polynomial at
        # both endpoints, whose trapezoidal average is exact on this segment.
        alpha=[1-sum(uv),uv...];beta=[-sum(delta),delta...]
        lo,hi=0.,1.
        for (a,b) in zip(alpha,beta)
            if iszero(b)
                a<0 && (hi=-1.)
            elseif b>0
                lo=max(lo,-a/b)
            else
                hi=min(hi,-a/b)
            end
        end
        hi>lo || continue
        d=collect(q).-collect(p)
        start=collect(p).+lo.*d;stop=collect(p).+hi.*d
        ax,ay,b=field.coefficients
        value(x)=[ax-b*x[2],ay+b*x[1]]
        term=sum((value(start).+value(stop)).*(stop.-start))/2
        total+=term
        uncertainty+=field.bound*sum(abs,stop.-start)+128eps(Float64)*abs(term)
    end
    return (;value=total,bound=uncertainty)
end

function ray(endpoint,centre,inner,outer,segments)
    cx,cy=centre;x,y=endpoint;d=x-cx
    nodes=[(cx+(outer-(outer-inner)*t)*t*d/inner,
        cy-(outer-(outer-inner)*t)*sqrt(1-(t*d/inner)^2))
        for t in range(0,1;length=segments+1)]
    push!(nodes,(x,y))
    return nodes
end

function richardson(values,budget)
    length(values)==4 || return (resolved=false,uncertainty=Inf)
    d=abs.(diff(values))
    all(>(0),d) && d[3]<d[2]<d[1] || return (resolved=false,uncertainty=Inf)
    p1,p2=log2(d[1]/d[2]),log2(d[2]/d[3])
    min(p1,p2)>=.5 && abs(p1-p2)<=.5 || return (resolved=false,uncertainty=Inf)
    uncertainty=d[3]/(2^min(p1,p2)-1)
    return (resolved=uncertainty<=budget/4,uncertainty)
end
end
