# Disposable 18 kV internal-shunt experiment.
# Run: include("test/manual/prototypes/prototype_wire_screen.jl")
# No engine extensions, FEM solves, campaign writes, or fitting to FEM data.
# A LOCAL layered-annulus Laplace Green function couples round-wire auxiliary
# sources to continuous, corner-weighted charge modes on all four tape faces. The core remains the
# analytical equivalent circle; the foil is the local voltage reference.
# This is a cross-section approximation, not a helical/full-wave field solver.
# Before: this experiment activated Gauntlet. Now preserve the active IDE project.
gauntlet_project = normpath(joinpath(@__DIR__, "..", "..", "..", "gauntlet"))
gauntlet_project in LOAD_PATH || push!(LOAD_PATH, gauntlet_project)
using LineCableModels, LinearAlgebra, DataFrames, JLD2, SHA, TOML, Test
isdefined(@__MODULE__, :Gauntlet) ||
    include(joinpath(gauntlet_project, "Gauntlet.jl"))

# ---- Manual controls ----
ws_campaign = joinpath(gauntlet_project, ".work", "all-references")
ws_case_id = :cable_18kv_1000mm2_trefoil
# Each named trial changes one numerical control from the baseline; last combines them.
ws_levels = (
    (label = "baseline", wire = 32, order = 16, quadrature = 128, modes = 512),
    (label = "tape order", wire = 32, order = 32, quadrature = 128, modes = 512),
    (label = "wire resolution", wire = 64, order = 16, quadrature = 128, modes = 512),
    (label = "quadrature", wire = 32, order = 16, quadrature = 256, modes = 512),
    (label = "kernel", wire = 32, order = 16, quadrature = 128, modes = 1024),
    (label = "combined", wire = 64, order = 32, quadrature = 256, modes = 1024))
ws_log_rtol = 1e-10
ws_penetration_target = 0.01  # 1% target on core-to-foil coupling, not the matrix norm.
ws_source_fraction = 0.75  # Numerical source location inside metal, NOT a geometry change.
ws_make_plots = true     # Set true in a graphical REPL; uses the owned plot recipe.

# All helper names are local to this disposable experiment, not engine methods.
function ws_read_result(path)
    digest = first(split(read(path * ".sha256", String)))
    bytes2hex(open(sha256, path)) == digest || error("Changed numerical payload: $path")
    jldopen(path, "r") do f
        f["basis"] == :pul || error("This experiment requires per-unit-length operands.")
        LineParameters(PhaseDomain, f["Z"], f["Y"], f["frequencies"]; basis = :pul)
    end
end

function ws_geometry(design)
    regions = design.geometry.regions
    bytag(tag) = only(filter(r -> r.source.tag == tag, regions))
    wires = [r.primitive for r in regions if r.source.tag == :sheath_wires]
    all(w -> w isa Disk, wires) || error("Prototype supports round screen wires only.")
    tape = bytag(:sheath_copper_tape).primitive
    filler = bytag(:screen_matrix_fill)
    host = filler.primitive.outer
    host isa Annulus || error("Screen filler must have a concentric annular envelope.")
    foil = bytag(:jacket_aluminum_tape).primitive
    inner = [bytag(tag)
             for tag in (:core_semicon_tape_inner, :core_semicon_inner,
        :core_insulation, :core_semicon_outer, :core_semicon_tape_outer)]
    outer = [bytag(:sheath_water_blocking)]
    all(r -> r.primitive isa Annulus, [inner; outer]) ||
        error("Concentric layers required.")
    all(r -> r.primitive.at.x == r.primitive.at.y == 0, [inner; outer]) ||
        error("Off-centre dielectric layers are outside this prototype.")
    layer(r) = (ri = r.primitive.ri, ro = r.primitive.ro, epsilon = r.source.material.eps_r)
    left, right = layer.(inner), layer.(outer)
    @assert last(left).ro ≈ host.ri
    @assert first(right).ri ≈ host.ro && last(right).ro ≈ foil.ri
    for layers in (left, right), i in 2:length(layers)

        @assert layers[i - 1].ro ≈ layers[i].ri
    end
    # Default owned dielectric laws are lossless with frequency-independent eps_r.
    # Do not silently extend this one-solve experiment to lossy/dispersive laws.
    epsilon = filler.source.material.eps_r
    Rleft = sum(log(l.ro/l.ri)/l.epsilon for l in left)
    Rright = sum(log(l.ro/l.ri)/l.epsilon for l in right)
    L = log(host.ro/host.ri)
    Rtotal = Rleft + L/epsilon + Rright
    (; wires, tape, left, right, a = host.ri, b = host.ro, epsilon,
        Rleft, Rright, Rtotal, core = first(left).ri, foil = foil.ri)
end

# Mode-m input normal admittance of a radial dielectric stack, terminated by
# a fixed-potential circle. In t=log(r), each layer solves v''-m^2*v=0.
# tanh recursion avoids overflow at large m or large dielectric contrasts.
function ws_load(layers, m; reverse_layers = false)
    load = Inf
    for layer in (reverse_layers ? reverse(layers) : layers)
        t = tanh(m * log(layer.ro/layer.ri))
        characteristic = layer.epsilon * m
        load = isinf(load) ? characteristic/t :
               characteristic * (load + characteristic*t)/(characteristic + load*t)
    end
    load
end

function ws_kernel_coefficients(g, modes)
    A, B, D = zeros(modes), zeros(modes), zeros(modes)
    ra0 = isempty(g.left) ? -1.0 :
          (g.epsilon-last(g.left).epsilon)/(g.epsilon+last(g.left).epsilon)
    rb0 = isempty(g.right) ? -1.0 :
          (g.epsilon-first(g.right).epsilon)/(g.epsilon+first(g.right).epsilon)
    for m in 1:modes
        l = ws_load(g.left, m)/g.epsilon
        r = ws_load(g.right, m; reverse_layers = true)/g.epsilon
        ra = isinf(l) ? -1.0 : (m-l)/(m+l)
        rb = isinf(r) ? -1.0 : (m-r)/(m+r)
        denominator = g.epsilon*m*(1-ra*rb*(g.a/g.b)^(2m))
        A[m] = ra/denominator-ra0/(g.epsilon*m)
        B[m] = rb/denominator-rb0/(g.epsilon*m)
        D[m] = ra*rb/denominator
    end
    (; A, B, D, ra0, rb0)
end

function ws_selfchecks()
    # Independent Dirichlet-annulus separated solution, not the reflection form.
    h = (; a = 1.0, b = 3.0, epsilon = 1.0, Rleft = 0.0, Rright = 0.0,
        Rtotal = log(3.0), left = NamedTuple[], right = NamedTuple[])
    k = ws_kernel_coefficients(h, 256)
    @testset "Local Green function controls" begin
        for (z, s) in ((1.3cis(0.2), 2.1cis(1.1)), (1.4cis(2.1), 2.5cis(-0.4)))
            x, y, L = log(abs(z)), log(abs(s)), log(3.0)
            exact = min(x, y)*(L-max(x, y))/L
            for m in 1:256
                exact += 2sinh(m*min(x, y))*sinh(m*(L-max(x, y))) /
                         (m*sinh(m*L))*cos(m*(angle(z)-angle(s)))
            end
            @test ws_kernel(z, s, h, k) ≈ exact atol=2e-14
            @test ws_kernel(z, s, h, k) ≈ ws_kernel(s, z, h, k) atol=2e-14
            @test abs(ws_kernel(cis(0.3), s, h, k)) < 2e-14
            @test abs(ws_kernel(3cis(0.3), s, h, k)) < 2e-14
        end
        h2 = merge(h, (epsilon = 2.0, Rtotal = log(3.0)/2))
        @test ws_kernel(1.4cis(0.2), 2.1cis(1.1), h2, ws_kernel_coefficients(h2, 256)) ≈
              ws_kernel(1.4cis(0.2), 2.1cis(1.1), h, k)/2
        targets, sources = [1.3cis(0.2), 1.4cis(2.1)], [2.1cis(1.1), 2.5cis(-0.4)]
        @test ws_kernel_matrix(targets, sources, h, k) ≈
              [ws_kernel(z, s, h, k) for z in targets, s in sources]
        regular = ws_kernel_matrix(targets,sources,h,k; regular = true)
        @test regular ≈ [ws_kernel(z,s,h,k; regular = true) for z in targets, s in sources]
        @test ws_kernel_matrix(targets,sources,h,k) ≈
            regular-[log(abs(z-s)) for z in targets, s in sources]
        @test isfinite(ws_kernel(2cis(0.3),2cis(0.3),h,k; regular = true))
        layers = [(ri = 1.0, ro = 1.2, epsilon = 3.0), (ri = 1.2, ro = 2.0, epsilon = 3.0)]
        @test ws_load(layers, 5) ≈ ws_load([(ri = 1.0, ro = 2.0, epsilon = 3.0)], 5)
        @test ws_load(layers, 1e-8) ≈ 3/log(2) rtol=1e-12

        # End-to-end charge extraction control: a complete metal annulus must
        # recover the two exact coaxial capacitances and perfect shielding.
        theta = 2pi .* (0:255) ./ 256
        targets = [1.8cis.(theta); 2.2cis.(theta)]
        sources = [1.98cis.(theta); 2.02cis.(theta)]
        kernel = ws_kernel_matrix(targets, sources, h, k)
        q = kernel \ hcat(-ws_core_voltage.(targets, Ref(h)), ones(length(targets)))
        screen = vec(sum(q; dims = 1))
        core = [1/h.Rtotal, 0.0] -
               vec(transpose(ws_core_voltage.(sources, Ref(h)))*q)
        capacitance = vcat(transpose(core), transpose(screen))
        left, right = 1/log(1.8), 1/log(3/2.2)
        @test capacitance ≈ [left -left; -left left+right] rtol=1e-9

    end
end

# Potential per auxiliary charge scaled by 2pi*epsilon0. The direct logarithm
# and the nearest dielectric images are exact; only the distant reflections are
# Fourier-truncated. This avoids slow convergence at the tape/dielectric contact.
# The zero mode is the exact radial series solution. Nonzero modes sum repeated
# reflections between the two layered Robin boundaries. No Sommerfeld kernel
# appears here: this domain ends at the CLOSED aluminium foil.
function ws_kernel(z, source, g, k; regular = false, split_images = false)
    r, s = abs(z), abs(source)
    Rr = g.Rleft + log(r/g.a)/g.epsilon
    Rs = g.Rleft + log(s/g.a)/g.epsilon
    result = min(Rr, Rs)*(g.Rtotal-max(Rr, Rs))/g.Rtotal
    if split_images
        result += (log(max(r,s))+k.ra0*log(s)+k.rb0*(2log(g.b)-log(r)))/g.epsilon
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
    @inbounds for m in eachindex(k.A)
        result += current*(k.A[m]*pa + k.B[m]*pb + k.D[m]*(pc+pe))
        previous, current = current, 2cosine*current-previous
        pa *= qa
        pb *= qb
        pc *= qc
        pe *= qe
    end
    result
end

function ws_kernel_matrix(targets, sources, g, k; regular = false, split_images = false)
    # Same Green function, batched into small Fourier blocks for BLAS. No dense
    # N-by-modes cache: working storage is only 4*64 columns per boundary set.
    zero_modes = merge(k, (; A = Float64[], B = Float64[], D = Float64[]))
    matrix = [ws_kernel(z, s, g, zero_modes; regular,split_images) for z in targets, s in sources]
    ta, tb = g.a ./ conj.(targets), targets ./ g.b
    sa, sb = g.a ./ conj.(sources), sources ./ g.b
    tpa, tpb = ones(ComplexF64, length(targets)), ones(ComplexF64, length(targets))
    spa, spb = ones(ComplexF64, length(sources)), ones(ComplexF64, length(sources))
    for first_mode in 1:64:length(k.A)
        modes = first_mode:min(first_mode + 63, length(k.A))
        U = Matrix{Float64}(undef, length(targets), 4length(modes))
        V = Matrix{Float64}(undef, length(sources), 4length(modes))
        for (column, m) in enumerate(modes)
            tpa .*= ta
            tpb .*= tb
            spa .*= sa
            spb .*= sb
            cross = k.D[m]*(g.a/g.b)^m
            @inbounds for i in eachindex(targets)
                U[i, 4column - 3]=real(tpa[i])
                U[i, 4column - 2]=imag(tpa[i])
                U[i, 4column - 1]=real(tpb[i])
                U[i, 4column]=imag(tpb[i])
            end
            @inbounds for j in eachindex(sources)
                a=k.A[m]*spa[j]+cross*spb[j]
                b=k.B[m]*spb[j]+cross*spa[j]
                V[j, 4column - 3]=real(a)
                V[j, 4column - 2]=imag(a)
                V[j, 4column - 1]=real(b)
                V[j, 4column]=imag(b)
            end
        end
        mul!(matrix, U, transpose(V), 1.0, 1.0)
    end
    matrix
end

ws_core_voltage(z, g) = 1-(g.Rleft+log(abs(z)/g.a)/g.epsilon)/g.Rtotal

# Round wires retain their existing auxiliary-source treatment. The tape has
# no source contour or inset: its physical faces are handled in the included file.
function ws_points(g, nw, fraction; shift = 0.0)
    targets, sources = ComplexF64[], ComplexF64[]
    for wire in g.wires
        centre = complex(wire.at.x, wire.at.y)
        for j in 0:(nw - 1)
            push!(targets, centre + wire.r*cis(2pi*(j+shift)/nw))
            push!(sources, centre + fraction*wire.r*cis(2pi*j/nw))
        end
    end
    (; targets, sources)
end

include(joinpath(@__DIR__, "prototype_tape_element.jl"))

function ws_replace_internal(baseline, C, expected_previous)
    f = frequencies(baseline)
    corrected = copy(observe(baseline, Y))
    # Core and screen voltages/charges relative to the foil; cancels both the
    # external jacket interval and the full existing earth contribution.
    H = [1.0 0.0; 0.0 1.0; -1.0 -1.0]
    replacement = inv(C)
    for (k, frequency) in pairs(f)
        s = 2pi*im*frequency
        P = s*inv(observe(baseline, Y)[:, :, k])
        updated = copy(P)
        for start in (1, 4, 7)
            indices = start:(start + 2)
            previous = transpose(H)*P[indices, indices]*H
            @test isapprox(inv(previous), expected_previous; rtol = 1e-7, atol = 1e-15)
            updated[start:(start + 1), start:(start + 1)] += replacement-previous
        end
        corrected[:, :, k] = s*inv(updated)
        foil_indices = [3, 6, 9]
        @assert updated[foil_indices, :] == P[foil_indices, :]
        @assert updated[:, foil_indices] == P[:, foil_indices]
        for left in (1:3, 4:6, 7:9), right in (1:3, 4:6, 7:9)

            left == right && continue
            @assert updated[left, right] == P[left, right]
        end
    end
    LineParameters(PhaseDomain, copy(observe(baseline, Z)), corrected, f; basis = :pul)
end

function ws_run(campaign, case_id, levels, source_fraction;
        log_rtol = ws_log_rtol, penetration_target = ws_penetration_target)
    ws_selfchecks()
    ws_tape_selfchecks()
    case_id == :cable_18kv_1000mm2_trefoil || error("This is a disposable 18 kV prototype.")
    folder = joinpath(campaign, "benchmark_18kv_1000mm2_trefoil_fem")
    state = TOML.parsefile(joinpath(folder, "state.toml"))
    state["state"] == "complete" || error("Select a completed saved benchmark.")
    attempt = joinpath(folder, state["current"])
    reference = ws_read_result(joinpath(attempt, "reference", "calculation.jld2"))
    baseline = ws_read_result(joinpath(attempt, "candidate", "points", "1", "calculation.jld2"))
    @assert frequencies(reference) == frequencies(baseline)
    original = deepcopy((observe(baseline,Z),observe(baseline,Y),observe(reference,Y)))
    loaded = Gauntlet.load_case(case_id)
    geometry = ws_geometry(first(loaded.problem.system.designs))
    geometry_check = ws_exposed_geometry_check(geometry)
    Cleft = 2pi*8.8541878128e-12/geometry.Rleft
    Cright = 2pi*8.8541878128e-12/geometry.Rright
    expected_previous = [Cleft -Cleft; -Cleft Cleft+Cright]
    println("\nFour-face tape prototype: ",length(geometry.wires)," wires + finite open tape.")
    println("Core remains an equivalent circle. Saved earth/Z retained. No FEM solve or campaign writes.")
    println("Exposed geometry check: ",geometry_check)
    println("Outer junction exponent: ",ws_junction_exponent(geometry.epsilon,first(geometry.right).epsilon))
    results = Any[]
    convergence_df = DataFrame()
    couplings(C) = (-C[1,2],sum(C[1,:]),sum(C[2,:]))
    for level in levels
        println("Solving ",level.label,": ",level)
        elapsed = @elapsed result = ws_spectral_capacitance(geometry,level;
            source_fraction,log_rtol,penetration_target)
        ccs,cca,csa = couplings(result.C)
        delta = isempty(results) ? (missing,missing,missing) :
            100 .* (([ccs,cca,csa]./collect(couplings(first(results).C))).-1)
        push!(results,result)
        push!(convergence_df,(refinement = level.label,wire_points = level.wire,
            tape_order = level.order,quadrature_nodes = level.quadrature,modes = level.modes,
            unknowns = result.unknowns,rank = result.rank,
            C_cs_nF_m = ccs*1e9,C_ca_nF_m = cca*1e9,C_sa_nF_m = csa*1e9,
            delta_cs_percent = delta[1],delta_ca_percent = delta[2],delta_sa_percent = delta[3],
            boundary_residual_V = result.residual,wire_residual_V = result.wire_residual,
            tape_residual_V = result.tape_residual,
            penetration_residual_V = result.penetration_residual,
            required_residual_V = result.required_residual,
            sampled_penetration_indicator_percent = 100result.sampled_penetration_relative_indicator,
            log_moment_error = result.log_moment_error,
            log_integration_indicator_V = result.log_integration_indicator_V,
            reciprocity = result.reciprocity,mutual_reciprocity = result.mutual_reciprocity,
            condition = result.condition,seconds = elapsed);promote = true)
    end
    println("\nChanges are against the baseline row; all three couplings are tracked individually.")
    show(stdout,MIME"text/plain"(),select(convergence_df,
        [:refinement,:C_cs_nF_m,:C_ca_nF_m,:C_sa_nF_m,:delta_ca_percent,
         :wire_residual_V,:tape_residual_V,:penetration_residual_V,:required_residual_V]);allcols = true)
    println()
    best = last(results)
    prototype = ws_replace_internal(baseline,best.C,expected_previous)
    @test eigmin(Symmetric((best.C+transpose(best.C))/2)) > 0
    @test best.mutual_reciprocity < 1e-3
    converged_boundary = best.penetration_residual <= best.required_residual
    if !converged_boundary
        @warn "Sampled boundary error still exceeds the core-to-foil coupling target." measured_V=best.penetration_residual required_V=best.required_residual
    end
    println("Sampled residuals are diagnostics, not a certified maximum; inspect the independent refinements.")
    couplings_df = DataFrame()
    for (label,value) in (("saved FEM",reference),("saved annular",baseline),("four-face tape",prototype)),
            k in (1,length(frequencies(value)))
        f = frequencies(value)[k]
        coupling = -imag.(observe(value,Y)[:,:,k])/(2pi*f)*1e9
        push!(couplings_df,(model = label,frequency_Hz = f,core_screen = coupling[1,2],
            core_foil = coupling[1,3],screen_foil = coupling[2,3]))
    end
    println("\nCouplings -B/omega [nF/m]; FEM values are apparent capacitances.")
    show(stdout,MIME"text/plain"(),couplings_df;allcols = true);println()
    summary_df = DataFrame()
    for (label,value) in (("saved annular",baseline),("four-face tape",prototype)),
            (term,i,j) in (("core-screen",1,2),("core-foil",1,3),("screen-foil",2,3))
        row = Dict{Symbol,Any}(:model=>label,:term=>term)
        for band in (:all,:dc,:harmonic,:narrow,:wide)
            comparison = Base.invokelatest(LineCableModels.Engine.compare,reference,value,B;band)
            metric = comparison.relative[i,j]
            row[band] = ismissing(metric) ? missing : 100metric
        end
        push!(summary_df,row;cols = :union)
    end
    select!(summary_df,[:model,:term,:all,:dc,:harmonic,:narrow,:wide])
    println("\nB relative RMS [%] against saved FEM. Missing: at least one series is near zero.")
    show(stdout,MIME"text/plain"(),summary_df;allcols = true);println()
    @test original == (observe(baseline,Z),observe(baseline,Y),observe(reference,Y))
    @test observe(prototype,Z) == observe(baseline,Z)
    @test all(isfinite,observe(prototype,Y))
    (;reference,baseline,prototype,geometry,geometry_check,results,convergence_df,
        couplings_df,summary_df,attempt,converged_boundary)
end

# One BLAS thread; restore the user's setting even if the prototype fails.
ws_previous_threads = BLAS.get_num_threads()
ws_experiment = try
    BLAS.set_num_threads(1)
    ws_run(ws_campaign, ws_case_id, ws_levels, ws_source_fraction)
finally
    BLAS.set_num_threads(ws_previous_threads)
end
ws_convergence_df = ws_experiment.convergence_df
ws_summary_df = ws_experiment.summary_df
ws_couplings_df = ws_experiment.couplings_df
ws_result = ws_experiment.prototype
ws_strip_result = ws_result  # The finite tape now uses the four-face charge expansion.
ws_reference = ws_experiment.reference
ws_baseline = ws_experiment.baseline
if ws_make_plots
    @eval using GLMakie
    # Before: a fabricated ParametricResult supplied plot labels. Now ordinary
    # completed-result collections delegate through observation construction.
    ws_candidates = [ws_baseline, ws_result]
    ws_plots = LineCableModels.plot(ws_candidates; reference = ws_reference,
        ydata = (G, B), layout = (3, 3), length_unit = :base,
        # Before: reference-first labels. Now candidates precede the separate reference.
        series_labels = ("Saved annular", "Four-face tape", "Saved FEM"))
end
nothing
