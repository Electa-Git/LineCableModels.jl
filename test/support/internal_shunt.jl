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
function internal_shunt_test_group(design)
    SM=LineCableModels.Engine.ShuntModel
    blueprint=LineCableModels.DataModel.flatten(design)
    decided=SM.dielectric_groups(SM.Formula(:boundary), design, blueprint,
        Formulation().formulas, 1)
    length(decided.groups)==1 || error("control requires exactly one resolved open-screen group")
    return only(decided.groups), blueprint
end

# A synthetic dielectric group of host radii `a` < `b` [m] and relative permittivity
# `epsilon`, referenced to grounded metal at its outermost radius.
function annular_group(a, b, epsilon; left = [], right = [], wires = [], tapes = [])
    SM=LineCableModels.Engine.ShuntModel
    outer=isempty(right) ? b : last(right).ro
    terminals=1:max(1, maximum(member -> member.terminal, [wires; tapes]; init = 1))
    return SM.DielectricGroup{Float64}(1, terminals, terminals,
        LineCableModels.Engine.ShuntJacket(Float64(outer)), a, b, epsilon,
        SM.ShuntLayer{Float64}[left...], SM.ShuntLayer{Float64}[right...],
        SM.ShuntWire{Float64}[wires...], SM.ShuntTape{Float64}[tapes...])
end

# Independent reference of the layered annular Green function: direct, nearest-image and
# zero-mode terms, plus the reflected Fourier modes summed one pair at a time. `regular`
# omits the direct logarithmic singularity.
function reference_kernel(z, source, g, k; regular = false, split_images = false)
    r, s = abs(z), abs(source)
    Rr = k.Rleft + log(r/g.a)/g.epsilon
    Rs = k.Rleft + log(s/g.a)/g.epsilon
    result = min(Rr, Rs)*(k.Rtotal-max(Rr, Rs))/k.Rtotal
    if split_images
        result += (log(max(r, s))+k.ra0*log(s)+k.rb0*(2log(g.b)-log(r)))/g.epsilon
    else
        result += (regular ? log(max(r, s)) : log(max(r, s)/abs(z-source)))/g.epsilon
        result -= k.ra0/g.epsilon*log(abs(1-g.a^2/(conj(z)*source)))
        result -= k.rb0/g.epsilon*log(abs(1-z*conj(source)/g.b^2))
    end
    qa, qb = g.a^2/(r*s), r*s/g.b^2
    qd = (g.a/g.b)^2
    qc, qe = qd*r/s, qd*s/r
    pa, pb, pc, pe = qa, qb, qc, qe
    cosine = real(z*conj(source))/(r*s)
    previous, current = 1.0, cosine
    for m in eachindex(k.A)
        result += current*(k.A[m]*pa + k.B[m]*pb + k.D[m]*(pc+pe))
        previous, current = current, 2cosine*current-previous
        pa *= qa
        pb *= qb
        pc *= qc
        pe *= qe
    end
    result
end

# The production kernel matrix of the boundary shunt formula, with its own buffers.
function kernel_matrix(selected, targets, sources, g, k; split_images = false)
    SM=LineCableModels.Engine.ShuntModel
    plan=(rows = length(targets), columns = 0, ports = 1, order = 0, sources,
        modes = length(k.A), quadrature = 0, derivatives = 0, ring = 0, block_rows = 0,
        block_columns = 0)
    buffers=SM.initialize_buffers(selected, Float64, g, plan, (;))
    return SM._shunt_kernel_matrix!(Matrix{Float64}(undef, length(targets), length(sources)),
        targets, sources, g, k; split_images, buffers)
end

# The production tape columns of the boundary shunt formula, with its own buffers.
function tape_columns(selected, targets, faces, g, k; integration)
    SM=LineCableModels.Engine.ShuntModel
    order=first(first(faces)).p
    nodes=maximum(tape -> maximum(face -> length(face.points), tape), faces)
    columns=sum(tape -> sum(face -> face.p+1, tape), faces)
    plan=(rows = length(targets), columns = 0, ports = 1, order, sources = ComplexF64[],
        modes = length(k.A), quadrature = nodes, derivatives = 0, ring = 0, block_rows = 0,
        block_columns = 0)
    buffers=SM.initialize_buffers(selected, Float64, g, plan, (;))
    return SM._shunt_tape_columns!(Matrix{Float64}(undef, length(targets), columns),
        targets, faces, g, k; integration, buffers)
end

# Independent check of a boundary solve: the largest residual of its terminal excitations
# on a shifted, three times finer grid of wire targets and denser tape-face targets that
# approach the corners.
function boundary_residual(selected, g, solved)
    SM=LineCableModels.Engine.ShuntModel
    plan=solved.workspace.plan
    count=3plan.level.wire
    points, terminals=ComplexF64[], Int[]
    for wire in g.wires, j in 0:(count - 1)
        point=complex(wire.x, wire.y)+wire.r*cis(2pi*(j+0.37)/count)
        all(other -> other === wire || other.terminal != wire.terminal ||
            abs(point-complex(other.x, other.y))-other.r > SM._shunt_tol(g.b), g.wires) || continue
        push!(points, point)
        push!(terminals, wire.terminal)
    end
    wire_rows=length(points)
    n=max(41, 12(plan.level.order+1))
    parameters=sort!(unique!([-cos.(pi .* ((1:n) .- 0.5) ./ n); collect(range(-1, 1; length = n+1));
        [-1+10.0^-k for k in 2:10]; [1-10.0^-k for k in 2:10]]))
    for tape in plan.faces, face in tape, t in parameters
        push!(points, SM._shunt_face_point(face.shape, t))
        push!(terminals, face.terminal)
    end
    check=merge(plan, (; rows = length(points), columns = 0, ring = 0, block_rows = 0,
        block_columns = 0))
    buffers=SM.initialize_buffers(selected, Float64, g, check, (;))
    collocation=Matrix{Float64}(undef, length(points), plan.columns)
    SM._shunt_matrix!(collocation, points, plan.sources, plan.faces, g, plan.k;
        plan.integration, buffers)
    residual=collocation*solved.coefficients-
             SM._shunt_rhs!(buffers.excitation, points, terminals, g, plan.k)
    errors=vec(maximum(abs, residual; dims = 2))
    return (; boundary = maximum(errors), wire = maximum(errors[1:wire_rows]; init = 0.0))
end
