# Exploratory audit of the supplied September 2026 half-space manuscript.
# This does not register formulas or change production assembly.
# Run with --project=test/gauntlet; --fem additionally runs the existing FEM.
# --legacy-baseline separately compares the retained production formulas.
# --cim-atol=VALUE explicitly supplies a dimensionless integral absolute tolerance.
# --pec selects the library PEC material for FEM and exactly zero analytical Zc.
using LineCableModels
using LinearAlgebra
using Printf
using SHA
using Gmsh
using LineCableModels.Engine: SpectralIntegral, integrate, computation_options
using LineCableModels.Engine: special_besselk
const SF = LineCableModels.Engine.SpecialFunctions
const JSON3 = LineCableModels.ImportExport.JSON3

module GauntletSupport
include(joinpath(@__DIR__, "runtime.jl"))
end

const AUDIT_FREQUENCIES = [0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0, 1000000.0]

relative(a, b) = norm(a - b) / max(norm(b), floatmin(Float64))
encoded(a) = (shape = collect(size(a)), real = vec(real.(a)), imag = vec(imag.(a)))

# Manuscript eq:internal-impedance, evaluated directly with a scaled I0/I1 ratio.
# The material is read from the case; no registered impedance formula is called.
function manuscript_internal(f, radius, material, temperature; Gamma = 0.0im)
    @assert temperature == material.T0 # This case is at the reference temperature.
    s = complex(0.0, 2pi*f)
    epsilon0, mu0 = 8.8541878128e-12, 4pi*1e-7
    shc = inv(material.rho) + s*epsilon0*material.eps_r
    kappac = sqrt(s*mu0*material.mu_r*shc-Gamma^2)
    x = kappac*radius
    return kappac/(2pi*radius*shc) * SF.besselix(0,x)/SF.besselix(1,x)
end

function interface_audit()
    worst = 0.0
    checks = 0
    # Unequal permeability, conducting air, arbitrary nonzero longitudinal Γ.
    omega = 3.0
    mu = (0.7, 1.4)
    sh = (0.2 + 0.3im, 2.0 + 0.6im)
    gamma = 0.03 + 0.2im
    k2 = (-im * omega) .* mu .* sh
    kap2 = .-gamma^2 .- k2
    for source in 1:2, lambda in (0.0, 0.01, 0.3, 3.0, 30.0)
        target = 3 - source
        e = source == 1 ? 1.0 : -1.0
        as, at = sqrt(lambda^2 + kap2[source]), sqrt(lambda^2 + kap2[target])
        ks, kt = kap2[source], kap2[target]
        ss, st = sh[source], sh[target]
        ms, mt = mu[source], mu[target]
        # Independently solve the four unsimplified interface equations for R,RH,T,TH.
        bc = ComplexF64[
            ks 0 -kt 0;
            0 ks 0 -kt;
            -im*gamma*lambda im*omega*ms*e*as im*gamma*lambda im*omega*mt*e*at;
            -e*ss*as -im*gamma*lambda -e*st*at im*gamma*lambda
        ]
        rhs = ComplexF64[-ks, 0, im*gamma*lambda, -e*ss*as]
        solved = bc \ rhs
        ratio = ks / kt
        det = (ss*as + st*at*ratio) * (ms*as + mt*at*ratio) +
              im*gamma^2*lambda^2*(1-ratio)^2/omega
        R = ((ss*as-st*at*ratio)*(ms*as+mt*at*ratio) -
             im*gamma^2*lambda^2*(1-ratio)^2/omega) / det
        RH = 2e*gamma*lambda*ss*as*(1-ratio)/(omega*det)
        TT, TH = ratio*(1+R), ratio*RH
        a0, ag = sqrt(lambda^2+kap2[1]), sqrt(lambda^2+kap2[2])
        dm, ds = mu[2]*a0+mu[1]*ag, sh[2]*a0+sh[1]*ag
        pairs = (
            (solved, [R, RH, TT, TH]),
            (ks*(1+R)/as, -2k2[source]*mt/dm -
                2gamma^2*ss*(ms*as+mt*at)/(dm*ds)),
            (1+R-omega*ms*lambda*RH/(gamma*e*as), 2ss*at/ds),
            (-at*TT/as-omega*mt*lambda*TH/(gamma*e*as), -2ss*as/ds)
        )
        for (actual, expected) in pairs
            error = relative(actual, expected)
            worst = max(worst, error)
            @assert error < 1e-10
            checks += 1
        end
    end
    return (; checks, worst_relative_residual = worst)
end

function buried_kernels(f, method; radius, depth, separation, rho, epsilon, mu,
        Gamma = 0.0im, scalar_diagnostic = false, integration_options = (;))
    @assert iszero(Gamma) # This audit's independent reflection check uses Γ=0.
    s = complex(0.0, 2pi*f)
    sh0, shg = s*epsilon[1], inv(rho)+s*epsilon[2]
    k02, kg2 = s*mu[1]*sh0-Gamma^2, s*mu[2]*shg-Gamma^2
    kg = sqrt(kg2)
    # Manuscript eq:Ar-Fr, using the target circumference.
    ar = SF.besseli(0, kg*radius)
    fr = 2pi*shg*radius*SF.besseli(1, kg*radius)/(kg*ar)
    controls = computation_options(SpectralIntegral,
        (method = method, options = merge(
            (rtol = method == :quad ? 1e-10 : 1e-6,), integration_options)))
    matrices = ntuple(_ -> zeros(ComplexF64, 2, 2), 4)
    Z, Pphi, H, Ldirect = matrices
    for (row, col) in ((1,1), (1,2))
        y = row == col ? 0.0 : separation
        distance = row == col ? radius : separation
        image_distance = row == col ? 2depth : hypot(separation, 2depth)
        direct = special_besselk(0, kg*distance)
        # Manuscript eq:Lambda-self and eq:Lambda-mutual.
        lambda_trace = (row == col ? direct : ar*direct) -
                       ar*special_besselk(0, kg*image_distance)
        function spectral(kind)
            kernel = lambda -> begin
                a0, ag = sqrt(lambda^2+k02), sqrt(lambda^2+kg2)
                dm, ds = mu[2]*a0+mu[1]*ag, shg*a0+sh0*ag
                decay = exp(-2depth*kg2/(ag+lambda))
                if kind === :Q
                    return decay*mu[1]/dm
                elseif kind === :Mphi
                    return decay*shg*(mu[2]*ag+mu[1]*a0)/(dm*ds)
                elseif kind === :MU
                    return decay*shg*a0/(ag*ds)
                else
                    # Γ=0 independently reduces the electric-Hertz reflection.
                    @assert iszero(Gamma)
                    return decay*(mu[1]*ag-mu[2]*a0)/(dm*ag)
                end
            end
            integral = SpectralIntegral(Val(:cosine), kernel,
                (height = 2depth, separation = y), abs(kg);
                angle = min(pi/4, atan(depth / max(y, eps()))))
            try
                return integrate(controls.method, integral, controls.options, nothing)
            catch error
                throw(ErrorException("$kind ($row,$col): $(sprint(showerror,error))"))
            end
        end
        Q, MU = spectral(:Q), spectral(:MU)
        Mphi = scalar_diagnostic ? spectral(:Mphi) : zero(Q)
        # Manuscript eq:Z-pair-same, eq:Pphi-pair-same, and eq:H-gg.
        Z[row,col] = s*mu[2]/(2pi)*(lambda_trace + 2ar*Q)
        Pphi[row,col] = s/(2pi*shg)*(lambda_trace + 2ar*Mphi)
        H[row,col] = s/(2pi*shg)*(lambda_trace + 2ar*MU)
        reflected = spectral(:R)
        Ldirect[row,col] = row == col ?
            kg*radius*special_besselk(1, kg*radius) -
                kg*radius*SF.besseli(1, kg*radius)*reflected :
            -kg*radius*SF.besseli(1, kg*radius)*(direct+reflected)
    end
    # This reuse is specific to the equal-radius, equal-depth benchmark.
    for matrix in matrices
        matrix[2,2] = matrix[1,1]
        matrix[2,1] = matrix[1,2]
    end
    # Manuscript eq:K-definition at Γ=0, eq:D-implementation, and global solves.
    K = Z
    L = Matrix{ComplexF64}(I, 2, 2)/ar-fr*K
    P = H/L
    Zexternal = (K+Gamma^2/s*H)/L
    Y = s*L/H
    # Candidate scalar convention uses the same physical-current map.
    Yphi = scalar_diagnostic ? s*L/Pphi : nothing
    return (; Zexternal, Y, Yphi, P, H, K, L,
        closure_error = relative(L, Ldirect),
        inverse_error = relative(Y*P, s*Matrix{ComplexF64}(I,2,2)),
        condition_L = cond(L), condition_H = cond(H), radius_argument = abs(kg*radius))
end

function main()
    BLAS.set_num_threads(1)
    pec = "--pec" in ARGS
    output = joinpath(@__DIR__, "..", "..", ".linecablemodels", "qa", "unified-earth-audit")
    pec && (output = joinpath(output,"pec"))
    mkpath(output)
    loaded = GauntletSupport.load_case(:two_bare_wires;
        variation = GauntletSupport.ExactOverrides((frequencies = AUDIT_FREQUENCIES,
            core_material = pec ? :pec : :copper)))
    problem = loaded.problem
    # Optional comparison only: neither these results nor their formula state
    # supplies any coefficient of the manuscript calculation below.
    native = "--legacy-baseline" in ARGS ?
        compute(problem, Formulation(); options = (trace = true,)) : nothing
    # Read geometry/materials directly, without preparing an earth-formula workspace.
    engine = LineCableModels.Engine
    blueprints = engine.CableBlueprint{Float64}[
        engine.flatten(engine.LineCableModelsCoaxial(), design, Float64)
        for design in problem.system.designs]
    core_materials = [only(blueprint.conductors).material for blueprint in blueprints]
    air, earth = problem.earth_props.layers
    @assert isinf(air.rho) # Ideal air in the supplied two-wire benchmark.
    epsilon0, mu0 = 8.8541878128e-12, 4pi*1e-7
    atol_index = findfirst(a -> startswith(a, "--cim-atol="), ARGS)
    cim_atol = atol_index === nothing ? 0.0 : parse(Float64, split(ARGS[atol_index], '=')[2])
    rows = Any[]
    raw = Any[]
    geometry = map(p -> p.nominal, loaded.definition.parameters)
    for (fi, f) in enumerate(AUDIT_FREQUENCIES)
        params = (radius = geometry.core_radius, depth = abs(geometry.cable_y),
            separation = abs(geometry.second_x-geometry.first_x),
            rho = earth.rho, epsilon = (epsilon0*air.eps_r,epsilon0*earth.eps_r),
            mu = (mu0*air.mu_r,mu0*earth.mu_r))
        internal = pec ? zeros(ComplexF64,2,2) :
            Diagonal([manuscript_internal(f, geometry.core_radius,
                material, problem.temperature) for material in core_materials])
        quad = nothing
        for method in (:quad, :trapz, :cim)
            try
                timed = @timed buried_kernels(f, method; params...,
                    scalar_diagnostic = method == :quad,
                    integration_options = method == :cim ? (atol = cim_atol,) : (;))
                result = timed.value
                method == :quad && (quad = result)
                Z = internal + result.Zexternal
                row = (; frequency = f, method, status = "ok", seconds = timed.time,
                    bytes = timed.bytes, result.closure_error, result.inverse_error,
                    result.condition_L, result.condition_H, result.radius_argument,
                    Z_vs_quad = relative(result.Zexternal, quad.Zexternal),
                    Y_vs_quad = relative(result.Y, quad.Y),
                    Z_mutual_vs_quad = relative(result.Zexternal[1,2],quad.Zexternal[1,2]),
                    Y_mutual_vs_quad = relative(result.Y[1,2],quad.Y[1,2]),
                    Z_vs_native = native === nothing ? nothing : relative(Z, native.Z.values[:,:,fi]),
                    Y_vs_native = native === nothing ? nothing : relative(result.Y, native.Y.values[:,:,fi]))
                push!(rows,row)
                push!(raw,(frequency=f, method, Z=encoded(Z), Zearth=encoded(result.Zexternal),
                    Y=encoded(result.Y),
                    Yphi=result.Yphi === nothing ? nothing : encoded(result.Yphi)))
                println(row)
                flush(stdout)
            catch error
                method == :quad && rethrow()
                row = (; frequency = f, method, status = "failed", error = sprint(showerror,error))
                push!(rows,row)
                println(row)
                flush(stdout)
            end
        end
    end
    report = (git_commit = readchomp(`git rev-parse HEAD`),
        audit_sha256 = bytes2hex(sha256(read(@__FILE__))),
        case_sha256 = bytes2hex(sha256(read(loaded.source_file))),
        candidate = (earth_kernels = "supplied manuscript",
            internal_impedance = pec ? "PEC: Zc=0" : "eq:internal-impedance",
            core_material = pec ? "pec" : "copper",
            core_resistivity = getproperty.(core_materials,:rho),
            radius_m = geometry.core_radius, depth_m = abs(geometry.cable_y),
            separation_m = abs(geometry.second_x-geometry.first_x), earth_resistivity = earth.rho,
            Gamma = "0", placement = "two equal bare buried wires", cim_atol),
        frequencies = AUDIT_FREQUENCIES, interface = interface_audit(),
        rows, raw, native = native === nothing ? nothing :
            (Z=encoded(native.Z.values), Y=encoded(native.Y.values)))
    write(joinpath(output,"analytical.json"), JSON3.write(report))
    if pec
        dense_frequencies = collect(10.0 .^ range(-1,6; length=281))
        params = (radius=geometry.core_radius, depth=abs(geometry.cable_y),
            separation=abs(geometry.second_x-geometry.first_x), rho=earth.rho,
            epsilon=(epsilon0*air.eps_r,epsilon0*earth.eps_r),
            mu=(mu0*air.mu_r,mu0*earth.mu_r))
        dense = map(dense_frequencies) do f
            result = buried_kernels(f,:quad; params...)
            (frequency=f, Z=encoded(result.Zexternal), Y=encoded(result.Y))
        end
        write(joinpath(output,"analytical-dense.json"), JSON3.write((
            frequencies=dense_frequencies, candidate=report.candidate, raw=dense)))
    end
    if "--fem" in ARGS
        fem = compute(problem, Formulation(:LineCableModelsFEM;
            fem_options=(gmsh_verbosity=1, getdp_verbosity=1, frequency_workers=2,
                keep_run_directory=true));
            options=(trace=true,))
        write(joinpath(output,"fem.json"), JSON3.write((
            frequencies=AUDIT_FREQUENCIES, Z=encoded(fem.Z.values),
            Y=encoded(fem.Y.values), core_material=pec ? "pec" : "copper",
            core_resistivity=getproperty.(core_materials,:rho), details=fem.details)))
        println("FEM results written to ", output)
    end
    println("Report: ", joinpath(output,"analytical.json"))
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
