"""
$(TYPEDSIGNATURES)

Write the native model data of `model` to `path`, with the mesh controls in `options`.
Use 17 significant digits for numbers and braces for numerical lists. Region dimensions
come from [`section`](@ref), repeated regions from `model.strand`, passive annuli from
`model.layer`, and terminal names from the model's labels.

The material data evaluates

```math
\\sigma = \\operatorname{Re}\\kappa, \\qquad
\\varepsilon = \\frac{\\operatorname{Im}\\kappa}{2\\pi f}, \\qquad
\\mu = \\mu_0\\mu_r,
```

Here ``\\kappa`` is admittivity \\[S/m\\], ``f`` is frequency \\[Hz\\],
``\\sigma`` is conductivity \\[S/m\\], and ``\\varepsilon`` is permittivity \\[F/m\\].
Permeability is ``\\mu`` \\[H/m\\], and vacuum permeability is ``\\mu_0`` \\[H/m\\].
The factor ``\\mu_r`` is relative permeability \\[dimensionless\\].

Earth uses the formulation's properties from `model.earth` at each frequency. Air uses the
declared properties of the first earth layer, written once. Both use resistivity ``\\rho``
\\[Ω·m\\] with ``\\sigma = 1/\\rho``, and relative permittivity ``\\varepsilon_r``
\\[dimensionless\\] with ``\\varepsilon = \\varepsilon_0\\varepsilon_r``, where
``\\varepsilon_0`` is vacuum permittivity \\[F/m\\]. Coordinates and lengths are in \\[m\\],
and longitudinal wavenumbers are in \\[1/m\\]. Return `path`.
"""
function Base.write(path::AbstractString, model::FieldModel,
        options::ComputationOptions{GmshBackend})
    system = model.problem.system
    frequencies = model.problem.frequencies
    air = model.problem.earth.layers[1]
    controls = options.overrides
    open(path, "w") do io
        println(io, "// Resolved LineCableModels FEM inputs; solver logic lives in getdp/*.pro")
        for (name, number) in (("NumTerminals", length(model.label)),
            ("NumCables", length(system.designs)),
            ("NumMaterialRegions", length(model.materials)),
            ("CADGeometryTolerance", 1e-8), ("FrequencyCount", length(frequencies)))
            println(io, name, " = ", @sprintf("%.17g", Float64(number)), ";")
        end
        for (name, entries) in (("Frequencies", frequencies),
            ("GammaReValues", real.(model.gamma)), ("GammaImValues", imag.(model.gamma)))
            numbers = [@sprintf("%.17g", Float64(number)) for number in entries]
            println(io, name, "() = {", join(numbers, ", "), "};")
        end
        for (name, number) in (("PmlSideThicknessFactor", controls.PmlSideThicknessFactor),
            ("PmlSideLayers", controls.PmlSideLayers),
            ("PmlSideGrading", controls.PmlSideGrading),
            ("PmlTopThicknessFactor", controls.PmlTopThicknessFactor),
            ("PmlTopLayers", controls.PmlTopLayers),
            ("PmlTopGrading", controls.PmlTopGrading),
            ("PmlBottomThicknessFactor", controls.PmlBottomThicknessFactor),
            ("PmlBottomLayers", controls.PmlBottomLayers),
            ("PmlBottomGrading", controls.PmlBottomGrading),
            ("DomainSizeFactor", options.domain_size_factor),
            ("MeshSizeFactor", options.mesh_size_factor),
            ("ExteriorMeshSizeFactor", controls.ExteriorMeshSizeFactor),
            ("InterfaceRefinementFactor", controls.InterfaceRefinementFactor),
            ("PmlQuadrangles", controls.PmlQuadrangles),
            ("ThinLayerStrips", controls.ThinLayerStrips),
            ("VolumeQuadrature", controls.VolumeQuadrature),
            ("PhysicalVolumeQuadrature", controls.PhysicalVolumeQuadrature),
            ("FEMMumpsForwardErrorBudget", 0.01),
            ("PmlQuadrature", controls.PmlQuadrature),
            ("ConductorGeometryTolerance", controls.ConductorGeometryTolerance),
            ("ConductorSkinDepthElements", controls.ConductorSkinDepthElements),
            ("ConductorMeshGrowth", controls.ConductorMeshGrowth),
            ("ConductorSkinDepths", controls.ConductorSkinDepths),
            ("ConductorThicknessElements", controls.ConductorThicknessElements))
            println(io, name, " = ", @sprintf("%.17g", Float64(number)), ";")
        end
        radii = [support(design.outer) for design in system.designs]
        for (name, entries) in (("CableX", [at.x for at in system.placements]),
            ("CableY", [at.y for at in system.placements]), ("CableRadius", radii))
            numbers = [@sprintf("%.17g", Float64(number)) for number in entries]
            println(io, name, "() = {", join(numbers, ", "), "};")
        end
        println(io, "NumPhysicalRegions = ", @sprintf("%.17g", length(system.regions)), ";")
        for (i, region) in enumerate(system.regions)
            s = section(region)
            println(io, "RegionAnnularPassive~{", i - 1, "} = ", Int(model.layer[i]), ";")
            for (name, number) in (("RegionDimensionKinds", s.kind),
                ("RegionDimension1", s.a), ("RegionDimension2", s.b),
                ("RegionRepeated", model.strand[i]))
                println(io, name, "~{", i - 1, "}() = {",
                    @sprintf("%.17g", Float64(number)), "};")
            end
        end
        offset = 0
        for (k, design) in enumerate(system.designs)
            offset += length(design.regions)
            println(io, "CableSizingRegions~{", k - 1, "}() = {", offset - 1, "};")
        end
        println(io, "PmlReflection = ", @sprintf("%.17g", options.pml_reflection), ";")
        sampling_tag = ("VOLTAGE_SAMPLING_NODES", TAGS.sampling_nodes)
        for (name, tag) in (("AIR_EM", TAGS.air), ("EARTH_EM", TAGS.earth),
            ("AIR_PML", TAGS.air_pml), ("EARTH_PML", TAGS.earth_pml),
            ("DOMAIN_PML", TAGS.pml), ("OUTBND_EM", TAGS.outer_boundary),
            ("OUTBND_ELE_AIR", TAGS.outer_air_boundary),
            ("OUTBND_ELE_REF", TAGS.outer_earth_boundary),
            ("INTERFACE_AIR_SOIL", TAGS.interface),
            ("INNER_PML_BND", TAGS.pml_inner_boundary),
            ("TERMINAL", TAGS.terminal_base + 1),
            ("TERMINAL_CONTOUR", TAGS.terminal_contour_base + 1),
            ("GAUGE_TREE", TAGS.gauge_tree), sampling_tag,
            ("MEASUREMENT_LINE", TAGS.measurement_line_base + 1))
            println(io, name, " = ", @sprintf("%.17g", tag), ";")
        end
        receivers = [Int(model.air[terminal.cable]) for terminal in system.terminals]
        println(io, "ReceiverInAir() = {", join(receivers, ", "), "};")
        println(io, "CABLE_CONTOUR = ", TAGS.cable_contour_base + 1, ";")
        labels = ["\"" * replace(label, '\\' => "\\\\", '"' => "\\\"") * "\""
                  for label in model.label]
        println(io, "TerminalNames() = Str[", join(labels, ", "), "];")
        tags = [material.tag for material in model.materials]
        conductive = [material.conductor for material in model.materials]
        losses = [any(value -> !iszero(real(value)), material.admittivity)
                  for material in model.materials]
        permeabilities = [vacuum_permeability * material.mu_r
                          for material in model.materials]
        for (name, entries) in (("MaterialRegionTags", tags),
            ("MaterialIsConductor", conductive), ("MaterialHasLoss", losses),
            ("MaterialMu", permeabilities))
            numbers = [@sprintf("%.17g", Float64(number)) for number in entries]
            println(io, name, "() = {", join(numbers, ", "), "};")
        end
        for (i, material) in enumerate(model.materials)
            epsilon = imag.(material.admittivity) ./ (2π .* frequencies)
            for (name, entries) in
                (("Sigma", real.(material.admittivity)), ("Epsilon", epsilon))
                numbers = [@sprintf("%.17g", Float64(number)) for number in entries]
                println(io, "Material", name, "_", i, "() = {", join(numbers, ", "), "};")
            end
        end
        earth_sigma = [inv(state.rho) for state in model.earth]
        earth_epsilon = [vacuum_permittivity * state.eps_r for state in model.earth]
        earth_mu = [vacuum_permeability * state.mu_r for state in model.earth]
        for (name, entries) in (("EarthSigma", earth_sigma),
            ("EarthEpsilon", earth_epsilon), ("EarthMu", earth_mu))
            numbers = [@sprintf("%.17g", Float64(number)) for number in entries]
            println(io, name, "() = {", join(numbers, ", "), "};")
        end
        for (name, number) in (("AirSigma", inv(air.rho)),
            ("AirEpsilon", vacuum_permittivity * air.eps_r),
            ("AirMu", vacuum_permeability * air.mu_r),
            ("Xcenter", model.center), ("Ycenter", 0), ("Zcenter", 0))
            println(io, name, " = ", @sprintf("%.17g", Float64(number)), ";")
        end
        return
    end
    return path
end

"""
$(TYPEDSIGNATURES)

Return the real path of the executable named by `LINECABLEMODELS_GETDP` when set.
Otherwise, install the GetDP 3.5.0 artifact on first use and return its executable path.
Base throws when the selected path does not exist.
"""
function Base.realpath(::GetdpArtifact)
    if haskey(ENV, "LINECABLEMODELS_GETDP")
        return realpath(ENV["LINECABLEMODELS_GETDP"])
    end
    toml = joinpath(@__DIR__, "..", "..", "Artifacts.toml")
    root = ensure_artifact_installed("getdp", toml)
    folder = Sys.isapple() ? "getdp-3.5.0-MacOSX" : "getdp-3.5.0-Linux64"
    return realpath(joinpath(root, folder, "bin", "getdp"))
end

"""
$(TYPEDSIGNATURES)

Scan the frequencies with at most `options.workers` GetDP processes at once.
Each process uses one mesh and reuses its factorization for every terminal basis. Read its
outputs after it exits, keep its log under `logs/`, and remove its work directory. Write the
combined rows to `results/Z.tsv` and `results/P.tsv`.

# Arguments

- `directory`: run directory with `input/getdp/`, `input/model_data.pro`, and the meshes.
- `model`: field model with terminal labels and frequencies \\[Hz\\].
- `options`: executable, concurrency, and native solver controls.

# Returns

- `Z`: series impedance in \\[Ω/m\\], indexed by response, basis, and frequency.
- `P`: potential coefficient \\[Ω·m\\], with the same indices as `Z`.

# Errors

- Throws for an unavailable executable, a failed identity check, or a native path of
  256 bytes or more.
- Throws `ErrorException` for a failed process or output check.
- Rethrows `InterruptException` after all processes stop.

The error gives the frequency, run directory, log path, and last twenty log lines.
Process failures also give the exit code and signal. The run directory stays on failure.
"""
function scan(directory::AbstractString, model::FieldModel{T},
        options::ComputationOptions{GmshBackend}) where {T}
    executable = realpath(options.getdp_executable)
    isfile(executable) || error("getdp_executable: not a file: $executable")
    for (name, path) in (("getdp_executable", executable), ("model", "model.pro"),
        ("ModelDataPath", "model_data.pro"), ("RunDirectory", "."),
        ("RawDirectory", "results"), ("BasisListPath", "bases.pro"))
        ncodeunits(path) < 256 ||
            error("$name: native paths must be shorter than 256 bytes: $path")
    end
    for k in eachindex(model.problem.frequencies)
        path = @sprintf("../../mesh/frequency-%04d.msh", k)
        ncodeunits(path) < 256 ||
            error("mesh: native paths must be shorter than 256 bytes: $path")
    end
    output = try
        buffer = IOBuffer()
        run(pipeline(Cmd([executable, "-info"]); stdout = buffer, stderr = buffer))
        String(take!(buffer))
    catch exception
        exception isa InterruptException && rethrow()
        error("getdp_executable: failed to execute $executable -info: " *
              sprint(showerror, exception))
    end
    occursin("getdp", lowercase(output)) && occursin(r"\d+\.\d+", output) ||
        error("getdp_executable: identity check did not report GetDP: $executable")
    n = length(model.label)
    nf = length(model.problem.frequencies)
    Z = Array{Complex{T}, 3}(undef, n, n, nf)
    P = Array{Complex{T}, 3}(undef, n, n, nf)
    processes = Base.Process[]
    stopping = Ref(false)
    failures = Pair{Int, Exception}[]
    pool = @async asyncmap(eachindex(model.problem.frequencies); ntasks = options.workers) do k
        stopping[] && return
        try
            stem = @sprintf("frequency-%04d", k)
            folder = joinpath(directory, "work", stem)
            mkpath(joinpath(folder, "results", "jobs"))
            for source in readdir(joinpath(directory, "input", "getdp"); join = true)
                endswith(source, ".pro") && cp(source, joinpath(folder, basename(source)))
            end
            cp(joinpath(directory, "input", "model_data.pro"),
                joinpath(folder, "model_data.pro"))
            write(joinpath(folder, "bases.pro"),
                "RequestedBases() = {" * join(1:n, ", ") * "};\n")
            mesh_path = @sprintf("../../mesh/frequency-%04d.msh", k)
            controls = options.overrides
            tolerance = string(controls.MumpsBackwardErrorTolerance)
            arguments = [executable, "model.pro", "-solve", "LineCableModelsFEMScan",
                "-setnumber", "Physics", "1", "-setnumber", "PreparedGaugeTree", "1",
                "-msh", mesh_path, "-name", "solver", "-v", string(options.getdp_verbosity),
                "-setstring", "ModelDataPath", "model_data.pro",
                "-setstring", "RunDirectory", ".", "-setstring", "RawDirectory", "results",
                "-setstring", "BasisListPath", "bases.pro",
                "-setnumber", "FrequencyIndex", string(k),
                "-setnumber", "PlotFieldMaps", "0", "-setnumber", "ReuseFactorization", "1",
                "-setnumber", "LinearSolver", "0",
                "-setnumber", "MumpsOrdering", string(controls.MumpsOrdering),
                "-setnumber", "PetscPrealloc", string(controls.PetscPrealloc),
                "-setnumber", "MumpsErrorAnalysis", string(controls.MumpsErrorAnalysis),
                "-setnumber", "MumpsRefinementMax", string(controls.MumpsRefinementMax),
                "-setnumber", "MumpsBackwardErrorTolerance", tolerance]
            if options.getdp_verbosity >= 4
                append!(arguments, ["-cpu", "-ksp_view", "-log_view", ":petsc.log"])
            end
            threads = string(options.solver_threads)
            command = addenv(Cmd(Cmd(arguments); dir = folder),
                "OPENBLAS_NUM_THREADS" => threads, "OPENBLAS_DEFAULT_NUM_THREADS" => threads,
                "OMP_NUM_THREADS" => threads, "MKL_NUM_THREADS" => threads)
            log_path = joinpath(folder, "getdp.log")
            io = open(log_path, "w")
            try
                process = run(pipeline(command; stdout = io, stderr = io); wait = false)
                push!(processes, process)
                wait(process)
                if !success(process) && !stopping[]
                    error("exit code $(process.exitcode), signal $(process.termsignal)")
                end
            finally
                close(io)
            end
            stopping[] && return
            for path in readdir(joinpath(folder, "results", "jobs"); join = true)
                mv(path, joinpath(directory, "results", "jobs", basename(path)))
            end
            log_path = joinpath(directory, "logs", stem * ".log")
            mv(joinpath(folder, "getdp.log"), log_path)
            petsc_path = joinpath(folder, "petsc.log")
            if isfile(petsc_path)
                mv(petsc_path, joinpath(directory, "logs", stem * "-petsc.log"))
            end
            rm(folder; recursive = true)
            primitives = read(directory, model, k)
            Z[:, :, k] = primitives.Z
            P[:, :, k] = primitives.P
            if !iszero(controls.MumpsErrorAnalysis)
                text = read(log_path, String)
                blocks = eachmatch(r"(?ms)^KSP Object:.*?(?=^KSP Object:|\z)", text)
                for (basis, block) in enumerate(blocks)
                    pattern = r"RINFOG\(7\),\s*RINFOG\(8\)\s*\(backward error est\):\s*([^,\s]+),\s*(\S+)"
                    found = match(pattern, block.match)
                    if found === nothing
                        @warn "MUMPS diagnostics unavailable; inspect native log" frequency_hz=model.problem.frequencies[k] basis log=log_path
                        continue
                    end
                    omega1 = tryparse(Float64, found[1])
                    omega2 = tryparse(Float64, found[2])
                    if omega1 === nothing || omega2 === nothing
                        @warn "MUMPS diagnostics unavailable; inspect native log" frequency_hz=model.problem.frequencies[k] basis log=log_path
                        continue
                    end
                    backward_error = omega1 + omega2
                    target = controls.MumpsBackwardErrorTolerance
                    if !isfinite(backward_error) || backward_error > target
                        @warn "MUMPS backward-error target not met; result retained" frequency_hz=model.problem.frequencies[k] basis backward_error target log=log_path
                    end
                end
            end
        catch exception
            if !stopping[]
                if exception isa InterruptException
                    push!(failures, k => exception)
                else
                    stem = @sprintf("frequency-%04d", k)
                    log_path = joinpath(directory, "logs", stem * ".log")
                    if !isfile(log_path)
                        log_path = joinpath(directory, "work", stem, "getdp.log")
                    end
                    tail = isfile(log_path) ? join(last(readlines(log_path), 20), '\n') : ""
                    pattern = r"(?i)(unknown|syntax error).*?(GenerateRHS|SolveAgain)"
                    capability = occursin(pattern, tail)
                    requirement = capability ?
                                  "GetDP 3.5 with GenerateRHSGroup and SolveAgain is required.\n" :
                                  ""
                    message = "solve: frequency $(model.problem.frequencies[k]) Hz, " *
                              "run directory $directory, log $log_path\n" *
                              sprint(showerror, exception) * "\n" * requirement * tail
                    push!(failures, k => ErrorException(message))
                end
                stopping[] = true
            end
        end
        return
    end
    try
        while !istaskdone(pool) && !stopping[]
            sleep(0.02)
        end
    finally
        stopping[] = true
        deadline = time() + 2
        terminated = Set{Base.Process}()
        while !istaskdone(pool) || any(process_running, processes)
            for process in processes
                if process_running(process) &&
                   (!(process in terminated) || time() >= deadline)
                    kill(process, time() < deadline ? Base.SIGTERM : Base.SIGKILL)
                    push!(terminated, process)
                end
            end
            sleep(0.02)
        end
        for process in processes
            wait(process)
        end
        wait(pool)
    end
    isempty(failures) || throw(last(first(failures)))
    header = "frequency_index\tfrequency_hz\tresponse_terminal\tbasis_terminal\treal\timaginary"
    for name in ("Z", "P")
        open(joinpath(directory, "results", name * ".tsv"), "w") do io
            println(io, header)
            for k in eachindex(model.problem.frequencies), basis in 1:n

                stem = @sprintf("getdp-f%04d-b%04d", k, basis)
                path = joinpath(directory, "results", "jobs", stem * "-$name.tsv")
                for line in eachline(path)
                    println(io, line)
                end
            end
            return
        end
    end
    return (Z = Z, P = P)
end

"""
$(TYPEDSIGNATURES)

Read the primitive matrices of frequency `k` from the job files under `results/jobs/`.
Check each basis completion row and every response row for indices, finite values, and a
frequency within `16eps(T)` relative tolerance. Report the earth domain ceiling and exact
transverse cutoff from the native PML row as warnings.

# Arguments

- `directory`: run directory with the native output files.
- `model`: field model with terminal labels and frequencies \\[Hz\\].
- `k`: frequency index in the problem.

# Returns

- `Z`: series impedance in \\[Ω/m\\], indexed by response and basis terminal.
- `P`: potential coefficient \\[Ω·m\\], with the same indices as `Z`.

# Errors

- Throws `ErrorException` for an absent file or malformed row.
- Throws `ErrorException` for inconsistent indices or frequency, duplicate responses, or
  nonfinite values.

The error identifies the file and the failed check.
"""
function Base.read(directory::AbstractString, model::FieldModel{T}, k::Integer) where {T}
    n = length(model.label)
    frequency = model.problem.frequencies[k]
    Z = Matrix{Complex{T}}(undef, n, n)
    P = Matrix{Complex{T}}(undef, n, n)
    for basis in 1:n
        stem = @sprintf("getdp-f%04d-b%04d", k, basis)
        path = joinpath(directory, "results", "jobs", stem * ".done")
        isfile(path) || error("$path: completion file is missing")
        lines = readlines(path)
        length(lines) == 1 || error("$path: expected one completion row")
        parts = split(only(lines), '\t'; keepempty = true)
        length(parts) == 6 || error("$path: completion row must have six columns")
        entries = Float64[]
        for part in parts
            number = tryparse(Float64, part)
            number === nothing && error("$path: malformed completion value: $part")
            isfinite(number) || error("$path: nonfinite completion value: $part")
            push!(entries, number)
        end
        entries[1] == 2 || error("$path: completion version must be 2")
        entries[2] == k || error("$path: completion frequency index does not match $k")
        isapprox(T(entries[3]), frequency; rtol = 16eps(T), atol = zero(T)) ||
            error("$path: completion frequency does not match $frequency Hz")
        entries[4] == basis || error("$path: completion basis does not match $basis")
        entries[5] == n || error("$path: completion terminal count does not match $n")
        iszero(entries[6]) || error("$path: completion field map flag must be zero")
        for (name, target) in (("Z", Z), ("P", P))
            path = joinpath(directory, "results", "jobs", stem * "-$name.tsv")
            isfile(path) || error("$path: matrix file is missing")
            lines = readlines(path)
            length(lines) == n || error("$path: expected $n response rows")
            seen = falses(n)
            for (index, line) in enumerate(lines)
                parts = split(line, '\t'; keepempty = true)
                length(parts) == 6 || error("$path: row $index must have six columns")
                entries = Float64[]
                for part in parts
                    number = tryparse(Float64, part)
                    number === nothing && error("$path: malformed value in row $index")
                    isfinite(number) || error("$path: nonfinite value in row $index")
                    push!(entries, number)
                end
                entries[1] == k || error("$path: frequency index in row $index must be $k")
                isapprox(T(entries[2]), frequency; rtol = 16eps(T), atol = zero(T)) ||
                    error("$path: frequency in row $index does not match $frequency Hz")
                entries[3] in 1:n || error("$path: response in row $index is out of range")
                response = Int(entries[3])
                seen[response] && error("$path: duplicate response $response in row $index")
                entries[4] == basis || error("$path: basis in row $index must be $basis")
                real_part = T(entries[5])
                imaginary_part = T(entries[6])
                isfinite(real_part) && isfinite(imaginary_part) ||
                    error("$path: converted value in row $index is nonfinite")
                target[response, basis] = complex(real_part, imaginary_part)
                seen[response] = true
            end
            all(seen) || error("$path: one or more response terminals are missing")
        end
    end
    path = joinpath(directory, "results", "jobs", @sprintf("pml-f%04d.tsv", k))
    isfile(path) || error("$path: PML file is missing")
    lines = readlines(path)
    length(lines) == 1 || error("$path: expected one PML row")
    parts = split(only(lines), '\t'; keepempty = true)
    length(parts) == 25 || error("$path: PML row must have 25 columns")
    entries = Float64[]
    for part in parts
        number = tryparse(Float64, part)
        number === nothing && error("$path: malformed PML value: $part")
        isfinite(number) || error("$path: nonfinite PML value: $part")
        push!(entries, number)
    end
    entries[1] == k || error("$path: PML frequency index does not match $k")
    isapprox(T(entries[2]), frequency; rtol = 16eps(T), atol = zero(T)) ||
        error("$path: PML frequency does not match $frequency Hz")
    if !iszero(entries[19])
        @warn "earth too resistive for FEM domain sizing; results not qualified" frequency_hz=frequency
    end
    if !iszero(entries[6]) || !iszero(entries[7])
        @warn "exact transverse cutoff; the transverse problem is singular" frequency_hz=frequency
    end
    return (Z = Z, P = P)
end
