# Current local-domain inputs. The open screen permits testing wire and tape charge
# scattering. The independent annular operator controls are in internal_shunt.jl.
function internal_shunt_test_design(;
        radius = 0.0003, epsilon = 3.0, tapes = true, count = 6,
        suffix = "", core_radius = 0.002, wire_angle = 0.17)
    metal=Material(kind = :conductor, rho = 2e-8, mu_r = 1.0, T0 = 20.0, alpha = 0.004)
    insulation_material=Material(kind = :insulator, rho = Inf, eps_r = epsilon)
    host_material=Material(kind = :insulator, rho = Inf, eps_r = 1.0)
    a=core_radius+0.002
    screen_outer=a+2radius
    t=0.00015
    b=screen_outer+(tapes ? t : zero(radius))
    open_parts=AbstractCablePart[terminal(Symbol(:middle, suffix),
        Group(:round_screen, Region(:wire, Disk(radius), metal);
            pattern = Ring(count; r = a+radius, φ0 = wire_angle)))]
    if tapes
        push!(open_parts,
            terminal(Symbol(:middle, suffix),
                Group(:strip_screen, Region(:strip, Rectangle(0.0012, t), metal);
                    pattern = Ring(1; r = screen_outer+t/2))))
    end
    origin=Stack(
        terminal(Symbol(:inner, suffix), Region(:axial, Disk(core_radius), metal)),
        Region(:inner_dielectric, Shell(0.002), insulation_material),
        Enclosure(:screen_domain, Stack(open_parts); primitive = Annulus(a, b), fill = host_material),
        Region(:outer_dielectric, Shell(0.0008), insulation_material),
        terminal(Symbol(:reference, suffix), Region(:closed_screen, Shell(0.0002), metal)),
        Region(:cover, Shell(0.0007), insulation_material))
    return build(CableDesign, "open-screen-control"*suffix, origin)
end
# A single-core cable with a ring of screen wires and a ring of armour wires in two dielectric
# fills, separated by a polymer sheath, with no closed metal between them. The armour is the
# outermost conductor, so the jacket is the reference.
function screen_armour_design(; screen = 12, armour = 18, screen_angle = 0.0,
        armour_angle = 0.1, permittivities = (2.3, 2.0, 2.3, 3.0, 2.3), suffix = "")
    metal=Material(kind = :conductor, rho = 1.72e-8, mu_r = 1.0, T0 = 20.0, alpha = 0.00393)
    dielectric(eps_r)=Material(kind = :insulator, rho = Inf, eps_r = eps_r)
    insulation, filler, sheath, bedding, jacket = dielectric.(permittivities)
    core, screen_wire, armour_wire=0.006, 0.0008, 0.0015
    screen_inner=core+0.006
    armour_inner=screen_inner+2screen_wire+0.003
    origin=Stack(
        terminal(Symbol(:core, suffix), Region(:conductor, Disk(core), metal)),
        Region(:insulation, Shell(0.006), insulation),
        Enclosure(:screen_layer,
            Stack(terminal(Symbol(:screen, suffix),
                Group(:screen_wires, Region(:wire, Disk(screen_wire), metal);
                    pattern = Ring(screen; r = screen_inner+screen_wire, φ0 = screen_angle))));
            primitive = Annulus(screen_inner, screen_inner+2screen_wire), fill = filler),
        Region(:sheath, Shell(0.003), sheath),
        Enclosure(:armour_layer,
            Stack(terminal(Symbol(:armour, suffix),
                Group(:armour_wires, Region(:wire, Disk(armour_wire), metal);
                    pattern = Ring(armour; r = armour_inner+armour_wire, φ0 = armour_angle))));
            primitive = Annulus(armour_inner, armour_inner+2armour_wire), fill = bedding),
        Region(:jacket, Shell(0.004), jacket))
    return build(CableDesign, "screen-armour"*suffix, origin)
end

function internal_shunt_test_group(design)
    SM=LineCableModels.Engine.ShuntModel
    blueprint=LineCableModels.DataModel.flatten(design)
    decided=SM.dielectric_groups(SM.Formula(:boundary), design, blueprint,
        Formulation().formulas, 1)
    length(decided.groups)==1 || error("control requires exactly one resolved open-screen group")
    return only(decided.groups), blueprint
end

# A synthetic dielectric group of concentric `layers` (ri, ro, epsilon) [m], grounded at the
# innermost radius and referenced to grounded metal at the outermost radius. Each wire and tape
# carries its terminal and the index of its layer.
function layered_group(layers; wires = [], tapes = [])
    SM=LineCableModels.Engine.ShuntModel
    terminals=1:max(1, maximum(member -> member.terminal, [wires; tapes]; init = 1))
    return SM.DielectricGroup{Float64}(1, terminals, terminals,
        LineCableModels.Engine.ShuntJacket(Float64(last(layers).ro)),
        SM.ShuntLayer{Float64}[layers...], SM.ShuntWire{Float64}[wires...],
        SM.ShuntTape{Float64}[tapes...])
end

# A synthetic group whose open conductors lie in the layer of radii `a` < `b` [m] and relative
# permittivity `epsilon`, between the layers `left` and `right`.
function annular_group(a, b, epsilon; left = [], right = [], wires = [], tapes = [])
    layer=length(left)+1
    return layered_group([left...; (ri = a, ro = b, epsilon); right...];
        wires = [merge(wire, (; layer)) for wire in wires],
        tapes = [merge(tape, (; layer)) for tape in tapes])
end

# Independent reference of the layered Green function, by transfer matrices in BigFloat: the
# zero mode from the log-radius sums, and the modes 1 to `modes` from the solutions that
# vanish at the inner and the outer radius, carried through the layers with u and epsilon u'
# continuous and divided by their Wronskian. No logarithm is extracted. The series converges
# only for targets and sources on separated radii.
function reference_kernel(z, target_layer, source, source_layer, g; modes)
    setprecision(BigFloat, 512) do
        layers=g.layers
        t=[log(big(first(layers).ri)); [log(big(layer.ro)) for layer in layers]]
        epsilon=[big(layer.epsilon) for layer in layers]
        R=[big(0); cumsum([(t[k+1]-t[k])/epsilon[k] for k in eachindex(layers)])]
        position(r, layer)=R[layer]+(log(big(r))-t[layer])/epsilon[layer]
        rz, rs=position(abs(z), target_layer), position(abs(source), source_layer)
        result=min(rz, rs)*(R[end]-max(rz, rs))/R[end]
        lower, upper=rz <= rs ?
                     ((log(big(abs(z))), target_layer), (log(big(abs(source))), source_layer)) :
                     ((log(big(abs(source))), source_layer), (log(big(abs(z))), target_layer))
        advance(u, w, step, e, m)=(cosh(m*step)*u+sinh(m*step)*w/(e*m),
            e*m*sinh(m*step)*u+cosh(m*step)*w)
        for m in 1:modes
            # The solution that vanishes inside, carried to the lower point.
            u, w=big(0), big(1)
            for k in 1:(lower[2]-1)
                u, w=advance(u, w, t[k+1]-t[k], epsilon[k], m)
            end
            uL, wL=advance(u, w, lower[1]-t[lower[2]], epsilon[lower[2]], m)
            # The solution that vanishes outside, carried to the upper and then the lower point.
            u, w=big(0), big(1)
            for k in length(layers):-1:(upper[2]+1)
                u, w=advance(u, w, t[k]-t[k+1], epsilon[k], m)
            end
            u, w=advance(u, w, upper[1]-t[upper[2]+1], epsilon[upper[2]], m)
            uR=u
            for k in upper[2]:-1:(lower[2]+1)
                u, w=advance(u, w, t[k]-(k == upper[2] ? upper[1] : t[k+1]), epsilon[k], m)
            end
            start=lower[2] == upper[2] ? upper[1] : t[lower[2]+1]
            u, w=advance(u, w, lower[1]-start, epsilon[lower[2]], m)
            result+=2uL*uR/(wL*u-uL*w)*cos(m*(angle(z)-angle(source)))
        end
        Float64(result)
    end
end

# The production kernel matrix of the boundary shunt formula, with its own buffers.
function kernel_matrix(selected, targets, target_layers, sources, source_layers, g, k;
        split_images = false)
    SM=LineCableModels.Engine.ShuntModel
    plan=(rows = length(targets), columns = 0, ports = 1, order = 0, sources,
        modes = size(k.modal, 1), quadrature = 0, derivatives = 0, ring = 0, block_rows = 0,
        block_columns = 0)
    buffers=SM.initialize_buffers(selected, Float64, g, plan, (;))
    return SM._shunt_kernel_matrix!(Matrix{Float64}(undef, length(targets), length(sources)),
        targets, target_layers, sources, source_layers, k; split_images, buffers)
end

# The production tape columns of the boundary shunt formula, with its own buffers.
function tape_columns(selected, targets, target_layers, faces, g, k; integration)
    SM=LineCableModels.Engine.ShuntModel
    order=first(first(faces)).p
    nodes=maximum(tape -> maximum(face -> length(face.points), tape), faces)
    columns=sum(tape -> sum(face -> face.p+1, tape), faces)
    plan=(rows = length(targets), columns = 0, ports = 1, order, sources = ComplexF64[],
        modes = size(k.modal, 1), quadrature = nodes, derivatives = 0, ring = 0,
        block_rows = 0, block_columns = 0)
    buffers=SM.initialize_buffers(selected, Float64, g, plan, (;))
    return SM._shunt_tape_columns!(Matrix{Float64}(undef, length(targets), columns),
        targets, target_layers, faces, k; integration, buffers)
end

# Independent check of a boundary solve: the largest residual of its terminal excitations
# on a shifted, three times finer grid of wire targets and denser tape-face targets that
# approach the corners.
function boundary_residual(selected, g, solved)
    SM=LineCableModels.Engine.ShuntModel
    plan=solved.workspace.plan
    count=3plan.level.wire
    points, terminals, layers=ComplexF64[], Int[], Int[]
    for wire in g.wires, j in 0:(count - 1)
        point=complex(wire.x, wire.y)+wire.r*cis(2pi*(j+0.37)/count)
        all(other -> other === wire || other.terminal != wire.terminal ||
            abs(point-complex(other.x, other.y))-other.r >
            SM._shunt_tol(g.layers[wire.layer].ro), g.wires) || continue
        push!(points, point)
        push!(terminals, wire.terminal)
        push!(layers, wire.layer)
    end
    wire_rows=length(points)
    n=max(41, 12(plan.level.order+1))
    parameters=sort!(unique!([-cos.(pi .* ((1:n) .- 0.5) ./ n); collect(range(-1, 1; length = n+1));
        [-1+10.0^-k for k in 2:10]; [1-10.0^-k for k in 2:10]]))
    for tape in plan.faces, face in tape, t in parameters
        push!(points, SM._shunt_face_point(face.shape, t))
        push!(terminals, face.terminal)
        push!(layers, face.layer)
    end
    check=merge(plan, (; rows = length(points), columns = 0, ring = 0, block_rows = 0,
        block_columns = 0))
    buffers=SM.initialize_buffers(selected, Float64, g, check, (;))
    collocation=Matrix{Float64}(undef, length(points), plan.columns)
    SM._shunt_matrix!(collocation, points, layers, plan.sources, plan.source_layers,
        plan.faces, plan.k; plan.integration, buffers)
    residual=collocation*solved.coefficients-
             SM._shunt_rhs!(buffers.excitation, points, layers, terminals, plan.k)
    errors=vec(maximum(abs, residual; dims = 2))
    return (; boundary = maximum(errors), wire = maximum(errors[1:wire_rows]; init = 0.0))
end
