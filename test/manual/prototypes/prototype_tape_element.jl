# Included by prototype_wire_screen.jl. Four exact tape faces, not a general BEM
# backend. Jacobi-weighted continuous charge measures replace fictitious tape
# sources. All kernels, core/foil functions and round-wire sources stay owned by
# the parent experiment. References for the representation and integration:
# https://idus.us.es/server/api/core/bitstreams/aae30e40-44f5-4243-9c30-35d64123f0e5/content
# https://ars.copernicus.org/articles/9/39/2011/
# https://math.dartmouth.edu/~fastdirect/notes/quadr.pdf

function ws_jacobi!(values,t,alpha,beta)
    p = length(values)-1
    values[1] = 1
    p == 0 && return values
    values[2] = ((alpha+beta+2)*t+alpha-beta)/2
    s = alpha+beta
    for n in 1:p-1
        A = (2n+s+1)*(2n+s+2)/(2(n+1)*(n+s+1))
        B = (alpha^2-beta^2)*(2n+s+1)/(2(n+1)*(n+s+1)*(2n+s))
        C = (n+alpha)*(n+beta)*(2n+s+2)/((n+1)*(n+s+1)*(2n+s))
        values[n+2] = (A*t+B)*values[n+1]-C*values[n]
    end
    values
end

ws_jacobi(t,alpha,beta,p) = ws_jacobi!(Vector{Float64}(undef,p+1),t,alpha,beta)

function ws_gauss_jacobi(n,alpha,beta)
    alpha > -1 && beta > -1 && n > 1 || error("Invalid Jacobi rule.")
    s = alpha+beta
    diagonal = [(beta-alpha)/(s+2);
        [(beta^2-alpha^2)/((2k+s)*(2k+s+2)) for k in 1:n-1]]
    off = [2/(s+2)*sqrt((1+alpha)*(1+beta)/(s+3));
        [2/(2k+s)*sqrt(k*(k+alpha)*(k+beta)*(k+s)/
            ((2k+s-1)*(2k+s+1))) for k in 2:n-1]]
    decomposition = eigen(SymTridiagonal(diagonal,off))
    # Normalized measure w(t)dt / integral(w); no arclength factor belongs here.
    decomposition.values,vec(decomposition.vectors[1,:].^2)
end

function ws_junction_exponent(epsilon_host,epsilon_outer)
    epsilon_host > 0 && epsilon_outer > 0 || error("Positive dielectric permittivities required.")
    # Host sector pi/2, outer dielectric sector pi. Determinant has no cotangent
    # poles. Exclude the trivial nu=0 root and bracket the first positive root.
    determinant(nu) = epsilon_host*cospi(nu/2)*sinpi(nu) +
        epsilon_outer*sinpi(nu/2)*cospi(nu)
    left,right = 1e-10,1.0
    @assert determinant(left) > 0 && determinant(right) < 0
    for _ in 1:60
        mid = (left+right)/2
        if determinant(mid) > 0
            left = mid
        else
            right = mid
        end
    end
    (left+right)/2
end

function ws_exposed_geometry_check(g)
    # Whole-face bases require whole exposed faces. Reject overlaps rather than
    # silently impose boundary conditions inside a metal union or move geometry.
    tolerance = 64eps(Float64)*g.b
    centres = [complex(w.at.x,w.at.y) for w in g.wires]
    gaps = [abs(centres[i]-centres[j])-g.wires[i].r-g.wires[j].r
        for i in eachindex(centres) for j in 1:i-1]
    wire_gap = isempty(gaps) ? Inf : minimum(gaps)
    wire_gap >= -tolerance || error("Overlapping screen wires: exposed-union boundaries required.")
    radial = [g.tape.ri-abs(c)-w.r for (c,w) in zip(centres,g.wires)]
    minimum(radial) >= -tolerance || error("Wire/tape radial overlap: whole-face tape prototype is not applicable.")
    @assert all(abs(c)-w.r >= g.a-tolerance for (c,w) in zip(centres,g.wires))
    @assert g.a < g.tape.ri < g.tape.ro <= g.b+tolerance
    @assert 0 < g.tape.span < 2pi
    (; wire_wire_gap_m = wire_gap, wire_tape_radial_gap_m = minimum(radial),
        tolerance_m = tolerance)
end

ws_face_point(face,t) = face.kind == :arc ?
    face.radius*cis(face.phi+face.span*t/2) :
    (face.mid+face.half*t)*cis(face.phi)

function ws_tape_faces(g,p,quadrature)
    p >= 0 && quadrature >= max(2,p+1) || error("Quadrature needs at least order+1 nodes.")
    s = g.tape
    at_interface = isapprox(s.ro,g.b; atol = 64eps(Float64)*g.b,rtol = 0)
    nu_inner = 2/3
    nu_outer = at_interface ? ws_junction_exponent(g.epsilon,first(g.right).epsilon) : 2/3
    descriptions = [
        (;kind = :arc,radius = s.ri,phi = s.at.φ,span = s.span,mid = 0.0,half = 0.0,
            alpha = nu_inner-1,beta = nu_inner-1),
        (;kind = :arc,radius = s.ro,phi = s.at.φ,span = s.span,mid = 0.0,half = 0.0,
            alpha = nu_outer-1,beta = nu_outer-1),
        [(;kind = :end,radius = 0.0,phi = s.at.φ+sign*s.span/2,span = 0.0,
            mid = (s.ri+s.ro)/2,half = (s.ro-s.ri)/2,
            alpha = nu_outer-1,beta = nu_inner-1) for sign in (-1,1)]...]
    map(descriptions) do face
        nodes,weights = ws_gauss_jacobi(quadrature,face.alpha,face.beta)
        polys = transpose(reduce(hcat,[ws_jacobi(t,face.alpha,face.beta,p) for t in nodes]))
        weighted = weights.*polys
        loggamma = LineCableModels.Engine.SpecialFunctions.loggamma
        beta_norm = exp(loggamma(face.alpha+1)+loggamma(face.beta+1)-
            loggamma(face.alpha+face.beta+2))
        merge(face,(;p,nodes,weights,weighted,beta_norm,
            points = ws_face_point.(Ref(face),nodes)))
    end
end

function ws_face_projection(z,face)
    if face.kind == :arc
        radius = abs(z)
        angle_difference = atan(sin(angle(z)-face.phi),cos(angle(z)-face.phi))
        u = 2angle_difference/face.span
        delta = abs(radius-face.radius)/(face.radius*face.span/2)
    else
        rotated = z*cis(-face.phi)
        u = (real(rotated)-face.mid)/face.half
        delta = abs(imag(rotated))/face.half
    end
    u,delta
end

function ws_log_moments(z,face; rtol = 1e-10)
    u,delta = ws_face_projection(z,face)
    if hypot(max(abs(u)-1,0),delta) > 0.2
        return transpose(face.weighted)*log.(abs.(face.points.-z)),0.0
    end
    # Weighted logarithmic product integration: compute the log moments of the
    # Jacobi modes, independently of the smooth Green remainder. t=cos(theta)
    # absorbs endpoint weights; split at the logarithmic singularity/near peak.
    theta0 = acos(clamp(u,-1,1))
    splits = [0.0,theta0,pi]
    if delta > 0
        width = min(sqrt(delta),delta/max(sin(theta0),eps(Float64)))
        for sign in (-1,1), factor in (1,4)
            push!(splits,clamp(theta0+sign*factor*width,0,pi))
        end
    end
    sort!(unique!(splits))
    function integrand!(output,theta)
        t = cos(theta)
        difference = abs(u) <= 1 ?
            -2sin((theta+theta0)/2)*sin((theta-theta0)/2) : t-u
        distance = if face.kind == :arc
            # Exact opposed-arc distance, also valid in the self limit.
            hypot(abs(z)-face.radius,
                2sqrt(abs(z)*face.radius)*sin(face.span*difference/4))
        else
            hypot(face.half*difference,imag(z*cis(-face.phi)))
        end
        # Roundoff may identify the endpoint only after its contribution is
        # below integration accuracy; never manufacture a finite log(0) floor.
        distance > 0 || error("Unresolved logarithmic integration point.")
        weight = sin(theta/2)^(2face.alpha+1)*cos(theta/2)^(2face.beta+1)/face.beta_norm
        ws_jacobi!(output,t,face.alpha,face.beta)
        output .*= log(distance)*weight
    end
    # Reuse quadrature work vectors; modal integration must not allocate a new
    # polynomial vector at each of its thousands of function evaluations.
    parentmodule(LineCableModels.Engine.quadgk).quadgk!(integrand!,zeros(face.p+1),splits...;
        rtol,atol = rtol*0.01,order = max(7,cld(face.p+1,2)),maxevals = 100000)
end

function ws_tape_columns(targets,faces,g,k; log_rtol = 1e-10)
    columns = Matrix{Float64}(undef,length(targets),sum(f.p+1 for f in faces))
    error_bound = 0.0
    offset = 0
    for face in faces
        block = ws_kernel_matrix(targets,face.points,g,k; split_images = true)*face.weighted
        for (i,z) in pairs(targets)
            inner,outer = g.a^2/conj(z),g.b^2/conj(z)
            direct_coefficient = 1.0
            images = Tuple{ComplexF64,Float64}[]
            for (point,coefficient) in ((inner,k.ra0),(outer,k.rb0))
                if abs(point-z) <= 8eps(Float64)*g.b
                    direct_coefficient += coefficient
                else
                    push!(images,(point,coefficient))
                end
            end
            push!(images,(z,direct_coefficient))
            for (point,coefficient) in images
                coefficient == 0 && continue
                moments,err = ws_log_moments(point,face; rtol = log_rtol)
                block[i,:] .-= (coefficient/g.epsilon).*moments
                error_bound = max(error_bound,abs(coefficient/g.epsilon)*err)
            end
        end
        columns[:,offset+1:offset+face.p+1] .= block
        offset += face.p+1
    end
    (;columns,log_moment_error = error_bound)
end

function ws_face_targets(faces,n; validation = false)
    parameters = -cos.(pi.*((1:n).-0.5)./n)
    if validation
        # Independent uniform targets plus progressively closer corner targets.
        parameters = sort!(unique!([parameters;collect(range(-1,1;length = n+1));
            [-1+10.0^-k for k in 2:10];[1-10.0^-k for k in 2:10]]))
    end
    [ws_face_point(face,t) for face in faces for t in parameters]
end

function ws_spectral_capacitance(g,level; source_fraction = 0.75,log_rtol = 1e-10,
        validation_factor = 3,penetration_target = 0.01)
    0 < source_fraction < 1 || error("Wire source fraction must be between zero and one.")
    level.wire >= 4 && level.modes >= 1 || error("Insufficient wire/kernel resolution.")
    log_rtol > 0 && penetration_target > 0 || error("Positive accuracy targets required.")
    geometry_check = ws_exposed_geometry_check(g)
    k = ws_kernel_coefficients(g,level.modes)
    faces = ws_tape_faces(g,level.order,level.quadrature)
    wires = ws_points(g,level.wire,source_fraction)
    targets = [ws_points(g,2level.wire,source_fraction).targets;
        ws_face_targets(faces,max(24,4(level.order+1)))]
    tape_matrix = ws_tape_columns(targets,faces,g,k;log_rtol)
    K = hcat(ws_kernel_matrix(targets,wires.sources,g,k),tape_matrix.columns)
    rhs = hcat(-ws_core_voltage.(targets,Ref(g)),ones(length(targets)))
    norms = vec(sqrt.(sum(abs2,K;dims = 1)))
    scaled = K./transpose(norms)
    factor = svd(scaled)
    cutoff = max(size(scaled)...)*eps(Float64)*first(factor.S)
    retained = factor.S .> cutoff
    coefficients = (factor.V[:,retained]*
        ((transpose(factor.U[:,retained])*rhs)./factor.S[retained]))./norms
    # Exact total-charge moments: higher Jacobi modes integrate to zero.
    mass = [ones(length(wires.sources));reduce(vcat,[[1.0;zeros(f.p)] for f in faces])]
    core_moment = ws_core_voltage.(wires.sources,Ref(g))
    for face in faces
        append!(core_moment,face.kind == :arc ?
            [ws_core_voltage(first(face.points),g);zeros(face.p)] :
            transpose(face.weighted)*ws_core_voltage.(face.points,Ref(g)))
    end
    epsilon0 = 8.8541878128e-12
    screen = 2pi*epsilon0*vec(transpose(mass)*coefficients)
    core = [2pi*epsilon0/g.Rtotal,0.0] -
        2pi*epsilon0*vec(transpose(core_moment)*coefficients)
    C = vcat(transpose(core),transpose(screen))
    check = [ws_points(g,validation_factor*level.wire,source_fraction;shift = 0.37).targets;
        ws_face_targets(faces,max(41,validation_factor*4(level.order+1));validation = true)]
    tape_check = ws_tape_columns(check,faces,g,k;log_rtol)
    errors = hcat(ws_kernel_matrix(check,wires.sources,g,k),tape_check.columns)*coefficients -
        hcat(-ws_core_voltage.(check,Ref(g)),ones(length(check)))
    wire_count = validation_factor*level.wire*length(g.wires)
    residual = maximum(abs,errors)
    wire_residual = maximum(abs,errors[1:wire_count,:])
    tape_residual = maximum(abs,errors[wire_count+1:end,:])
    penetration_residual = maximum(abs,vec(sum(errors;dims = 2)))
    ccs,cca,csa = -C[1,2],sum(C[1,:]),sum(C[2,:])
    ccs_scale = max(abs(C[1,2]),abs(C[2,1]))
    required_residual = penetration_target*abs(cca)/ccs_scale
    moment_error = max(tape_matrix.log_moment_error,tape_check.log_moment_error)
    log_integration_indicator_V = 6moment_error*
        maximum(vec(sum(abs.(coefficients[length(wires.sources)+1:end,:]);dims = 1)))
    (;C,q = coefficients,residual,wire_residual,tape_residual,penetration_residual,required_residual,
        sampled_penetration_relative_indicator = penetration_residual*ccs_scale/abs(cca),
        reciprocity = norm(C-transpose(C))/norm(C),
        mutual_reciprocity = abs(C[1,2]-C[2,1])/max(abs(C[1,2]),abs(C[2,1])),
        unknowns = length(mass),rank = count(retained),
        condition = first(factor.S)/last(factor.S),level,tape = true,model = :spectral,
        worst_point = check[argmax(abs.(errors))[1]],faces,geometry_check,
        log_moment_error = moment_error,log_integration_indicator_V)
end

function ws_tape_selfchecks()
    @testset "Four-face tape controls" begin
        @test ws_junction_exponent(1.0,1.0) ≈ 2/3
        @test ws_junction_exponent(1.0,32.3) ≈ 2/pi*acos(sqrt(32.3/(2*33.3)))
        @test 0.5 < ws_junction_exponent(1.0,32.3) < 2/3
        for (alpha,beta) in ((0.0,0.0),(-1/3,-1/3),(-0.4904397072754232,-1/3))
            nodes,weights = ws_gauss_jacobi(48,alpha,beta)
            values = reduce(hcat,[ws_jacobi(t,alpha,beta,8) for t in nodes])
            @test sum(weights) ≈ 1
            @test all(weights .> 0)
            @test dot(nodes,weights) ≈ (beta-alpha)/(alpha+beta+2) atol=1e-14
            @test values*weights ≈ [1.0;zeros(8)] atol=1e-13
        end
        nodes,weights = ws_gauss_jacobi(64,0.0,0.0)
        values = transpose(reduce(hcat,[ws_jacobi(t,0.0,0.0,3) for t in nodes]))
        face = (;kind = :end,radius = 0.0,phi = 0.0,span = 0.0,mid = 2.0,half = 0.1,
            alpha = 0.0,beta = 0.0,p = 3,nodes,weights,beta_norm = 1.0,
            weighted = weights.*values,points = complex.(2 .+ 0.1.*nodes))
        moments,_ = ws_log_moments(2.0+0im,face)
        @test moments ≈ [log(0.1)-1,0,1/3,0] atol=2e-10
        endpoint,_ = ws_log_moments(2.1+0im,face)
        @test endpoint[1] ≈ log(0.2)-1 atol=2e-10
        for distance in (0.1,1e-5)
            near,_ = ws_log_moments(2.0+im*distance,face)
            exact = log(hypot(0.1,distance))-1+distance/0.1*atan(0.1/distance)
            @test near[1] ≈ exact atol=2e-10
        end
        h = (;a = 1.0,b = 3.0,epsilon = 1.0,Rleft = 0.0,Rright = log(4/3)/4,
            Rtotal = log(3)+log(4/3)/4,left = NamedTuple[],
            right = [(ri = 3.0,ro = 4.0,epsilon = 4.0)])
        k = ws_kernel_coefficients(h,512)
        for (z,s) in ((1.3cis(0.2),2.1cis(1.1)),(3cis(0.2),2.5cis(-0.4)))
            split = ws_kernel(z,s,h,k;split_images = true)-
                (log(abs(z-s))+k.ra0*log(abs(h.a^2/conj(z)-s))+
                 k.rb0*log(abs(h.b^2/conj(z)-s)))/h.epsilon
            @test split ≈ ws_kernel(z,s,h,k) atol=2e-14
        end
        # Direct and coincident nearest image produce the interface coefficient.
        difference = ws_kernel(3.0+0im,3cis(1e-7),h,k)-
            ws_kernel(3.0+0im,3cis(2e-7),h,k)
        @test difference ≈ (2/5)*log(2) rtol=1e-5
        # The charge measure contains no extra face-length Jacobian. Scaling all
        # radii and target/source positions leaves the potential columns unchanged.
        targets = [1.4cis(0.1),2.5cis(0.3)]
        original = ws_tape_columns(targets,[face],h,k).columns
        scaled_h = merge(h,(a = 10h.a,b = 10h.b,
            right = [(ri = 30.0,ro = 40.0,epsilon = 4.0)]))
        scaled_face = merge(face,(mid = 10face.mid,half = 10face.half,points = 10face.points))
        scaled = ws_tape_columns(10targets,[scaled_face],scaled_h,
            ws_kernel_coefficients(scaled_h,512)).columns
        @test scaled ≈ original atol=1e-12
    end
end
