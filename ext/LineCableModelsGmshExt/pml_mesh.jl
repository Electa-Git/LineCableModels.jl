# Prescribed PML mesh construction. No discrete field solve, error acceptance,
# reference comparison, retry or refinement is performed here.
_pml_stretch(u, b, slope=1.) = 1 + (1-im*slope)*b*u^3
_pml_stretch_derivative(u, b, slope=1.) = 3(1-im*slope)*b*u^2
_pml_coordinate(u, b, slope=1.) = u + (1-im*slope)*b*u^4/4

function _pml_modal_response(Q, b, slope=1.)
    z = _pml_coordinate(1., b, slope)
    iszero(Q) && return inv(z)
    return Q*(1+exp(-2Q*z))/(-expm1(-2Q*z))
end

function _pml_modal_curvature(u, Q, b, slope=1.)
    s, ds = _pml_stretch(u,b,slope), _pml_stretch_derivative(u,b,slope)
    z, zend = _pml_coordinate(u,b,slope), _pml_coordinate(1.,b,slope)
    iszero(Q) && return -ds/zend
    denominator = -expm1(-2Q*zend)
    tail, outgoing = exp(-2Q*(zend-z)), exp(-Q*z)
    v = outgoing*(-expm1(-2Q*(zend-z)))/denominator
    dv = -Q*s*outgoing*(1+tail)/denominator
    return Q^2*s^2*v + ds/s*dv
end

function _pml_modes(air, earth, thickness, clearance, direction)
    media = direction == 1 ? (air,earth) : direction == 2 ? (air,) : (earth,)
    spectrum = [(Q=0.0im, weight=1.)]
    tangential = sort!(unique([0.; imag(air).*[.5,.9,.99,.999,1.,1.001,1.01,1.1,2.];
        exp.(range(log(.01/clearance),log(12/clearance);length=30))]))
    for gamma in media, k in tangential
        # Inputs already contain the prescribed-Gamma transverse wavenumbers.
        q = sqrt(complex(k^2+gamma^2))
        real(q) < 0 && (q = -q)
        push!(spectrum, (Q=q*thickness, weight=exp(-2real(q)*clearance)))
    end
    return spectrum
end

function _pml_density_integral(u, density)
    integral = zeros(length(u))
    for i in 2:length(u)
        integral[i] = integral[i-1] + (u[i]-u[i-1])*(density[i]+density[i-1])/2
    end
    return integral
end

function _pml_density_grid(wave_numbers, thickness, strength, clearance, direction, resolution; slope=1.)
    spectrum = _pml_modes(wave_numbers..., thickness, clearance, direction)
    transition = min(1., strength^(-1/3))
    # Fixed quadrature for the prescribed density, resolving b*u^3 ≈ 1.
    # These samples are neither field-mesh nodes nor trial FEM solutions.
    u = sort!(unique([0.; exp.(range(log(transition*1e-5),0.;length=2501));
        collect(range(0.,1.;length=2501))]))
    density = zeros(length(u))
    for mode in spectrum
        scale = abs(_pml_modal_response(mode.Q,strength,slope))
        for i in eachindex(u)
            weight = mode.weight*abs2(_pml_modal_curvature(u[i],mode.Q,strength,slope))/
                (abs(_pml_stretch(u[i],strength,slope))*scale)
            density[i] = max(density[i],cbrt(weight))
        end
    end
    integral = _pml_density_integral(u,density)
    density = max.(density.*(resolution.interpolation_cells/integral[end]),
        abs.(_pml_stretch_derivative.(u,strength,slope)./_pml_stretch.(u,strength,slope))./
        resolution.coefficient_change)
    return (;u, integral=_pml_density_integral(u,density))
end

function _pml_equal_density_nodes(grid, count)
    nodes = [0.]
    for target in range(0.,grid.integral[end];length=count+1)[2:end-1]
        j = searchsortedfirst(grid.integral,target)
        fraction = (target-grid.integral[j-1])/(grid.integral[j]-grid.integral[j-1])
        push!(nodes,grid.u[j-1]+fraction*(grid.u[j]-grid.u[j-1]))
    end
    push!(nodes,1.)
    return nodes
end

"""
$(TYPEDSIGNATURES)

Construct native geometric strips from a prescribed modal interpolation
density and logarithmic stretch variation. `wave_numbers` are the air/earth
propagation constants \\[1/m\\], `thickness` and `clearance` are distances
\\[m\\], and `strength` is the cubic stretch coefficient \\[dimensionless\\].
`direction` selects side, top or bottom with indices 1, 2 or 3.
`slope` is the negative imaginary-to-real ratio of the complex stretch increment.

Return normalized strip endpoints, interval counts and progression ratios.
The construction is fixed: sample the density, integrate, split at six equal
increments and the stretch transition, then prescribe geometric progressions.
It does not estimate or enforce the scientific error of the computed fields.
"""
function _physical_pml_strips(wave_numbers, thickness, strength, clearance, direction, resolution; slope=1.)
    grid = _pml_density_grid(wave_numbers,thickness,strength,clearance,direction,resolution;slope)
    boundaries = sort!(unique([_pml_equal_density_nodes(grid,6); min(1.,strength^(-1/3))]))
    function cumulative(u)
        i = clamp(searchsortedlast(grid.u,u),1,length(grid.u)-1)
        t = (u-grid.u[i])/(grid.u[i+1]-grid.u[i])
        return (1-t)*grid.integral[i]+t*grid.integral[i+1]
    end
    function position(target)
        i = clamp(searchsortedlast(grid.integral,target),1,length(grid.integral)-1)
        t = (target-grid.integral[i])/(grid.integral[i+1]-grid.integral[i])
        return (1-t)*grid.u[i]+t*grid.u[i+1]
    end
    strips = FEMPMLStrip{Float64}[]
    for (start,stop) in zip(boundaries[1:end-1],boundaries[2:end])
        lo, hi = cumulative(start), cumulative(stop)
        count = max(3,ceil(Int,hi-lo))
        count < typemax(Cint) || throw(ArgumentError("resolved PML interval count exceeds Gmsh's Cint range"))
        first_size = position(lo+(hi-lo)/count)-start
        last_size = stop-position(hi-(hi-lo)/count)
        ratio = (last_size/first_size)^(1/(count-1))
        push!(strips,FEMPMLStrip(start,stop,count,ratio))
    end
    return strips
end
