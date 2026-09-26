# Independent boundary construction for the rounded disk/wedge intersection.
# Circle centres follow signed-distance constraints to the two wedge lines;
# Green integrals establish area and first moments without production contacts.
module SectorControl
using QuadGK
function boundary(span,base,back,fillet)
    pi=oftype(base,π)
    theta=span/2
    axis=(base+fillet,zero(base))
    along=-(base+fillet)*cos(theta)+sqrt((back-fillet)^2-((base+fillet)*sin(theta))^2)
    upper=axis.+along.*(cos(theta),sin(theta));lower=(upper[1],-upper[2])
    upnormal=(-sin(theta),cos(theta));downnormal=(-sin(theta),-cos(theta))
    upbase=axis.+fillet.*upnormal;downbase=axis.+fillet.*downnormal
    uplink=upper.+fillet.*upnormal;downlink=lower.+fillet.*downnormal
    beta=atan(upper[2],upper[1])
    arc(c,r,a,b)=t->begin
        phi=a+(b-a)*t
        (c[1]+r*cos(phi),c[2]+r*sin(phi),-r*sin(phi)*(b-a),r*cos(phi)*(b-a))
    end
    line(p,q)=t->(p[1]+t*(q[1]-p[1]),p[2]+t*(q[2]-p[2]),q[1]-p[1],q[2]-p[2])
    return (line(downbase,downlink),arc(lower,fillet,3pi/2-theta,2pi-beta),
        arc((zero(base),zero(base)),back,-beta,beta),arc(upper,fillet,beta,pi/2+theta),
        line(uplink,upbase),arc(axis,fillet,pi/2+theta,3pi/2-theta))
end
function moments(span,base,back,fillet)
    curves=boundary(span,base,back,fillet)
    sums=zeros(typeof(base),4);errors=zeros(typeof(base),4)
    goal=precision(base)<=128 ? 1e-8 : precision(base)<=256 ? 1e-10 : 1e-12
    for curve in curves,index in 1:4
        value,estimate=quadgk(zero(base),one(base);rtol=oftype(base,goal),maxevals=10^6) do t
            x,y,dx,dy=curve(t)
            ((x*dy-y*dx)/2,x^2*dy/2,-y^2*dx/2,hypot(dx,dy))[index]
        end
        sums[index]+=value;errors[index]+=estimate
    end
    return (area=sums[1],centroid=(sums[2]/sums[1],sums[3]/sums[1]),
        perimeter=sums[4],errors=(area=errors[1],perimeter=errors[4],
            centroid=((errors[2]+abs(sums[2]/sums[1])*errors[1])/(sums[1]-errors[1]),
                (errors[3]+abs(sums[3]/sums[1])*errors[1])/(sums[1]-errors[1]))))
end
end
