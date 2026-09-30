# Qualification only: exact metal-drive substitution w = u - Gamma^2 V.
# v is constant in each metal. The two metal v couplings therefore cancel
# algebraically after substituting u = w + Gamma^2 V. Physical output restores u.
# Frozen sources, meshes and production code remain untouched.
isdefined(@__MODULE__, :ROOT) || include("qualify.jl")

function terminal_shift_control(source; suffix="-terminal-shift", extra_args=String[])
    target = source * suffix
    if !isfile(joinpath(target, "sources.toml"))
        files = readlines(joinpath(source, ".onelab-export-files"))
        for file in [files; ".onelab-export-files"; "study.msh"; "mesh.toml"]
            mkpath(dirname(joinpath(target, file)))
            cp(joinpath(source, file), joinpath(target, file); force=true)
        end
        path = joinpath(target, "formulations/quasi-full.pro")
        text = read(path, String)
        old = join([
            "        Galerkin { [-gamma2[] * seZ[] * Dof{v} * Vector[0,0,1], {a}];",
            "          In DomainFields; Jacobian Vol; Integration I1; }",
            "        If(!PerfectConductors)",
            "          Galerkin { [-gamma2[] * seZ[] * Dof{v} * Vector[0,0,1], {ur}];",
            "            In ConductorMaterialRegions; Jacobian Vol; Integration I1; }",
            "        EndIf"], '\n')
        @assert occursin(old, text)
        text = replace(text, old => """        Galerkin { [-gamma2[] * seZ[] * Dof{v} * Vector[0,0,1], {a}];
          In DomainMedia_Ele; Jacobian Vol; Integration I1; }""")
        text = replace(text,
            "-Re[{U} / UnitSource]" => "-Re[({U}+gamma2[]*{V}) / UnitSource]",
            "-Im[{U} / UnitSource]" => "-Im[({U}+gamma2[]*{V}) / UnitSource]")
        # Restore the physical drive also in optional field-map expressions.
        text = replace(text, "  gamma2[] = gamma[] * gamma[];" =>
            "  gamma2[] = gamma[] * gamma[];\n  driveShift[DomainMedia_Ele] = 0.;\n  driveShift[ConductorMaterialRegions] = gamma2[];")
        parts = split(text, "PostProcessing {"; limit=2)
        text = parts[1] * "PostProcessing {" * replace(parts[2],
            "{ur}" => "({ur}+driveShift[]*{v}*Vector[0,0,1])")
        write(path, text)
        record(joinpath(target, "sources.toml"), Dict(file=>digest(joinpath(target,file)) for file in files))
    end
    @assert digest(joinpath(target,"study.msh")) == digest(joinpath(source,"study.msh"))
    solve!(target; extra_args)
    return target
end

if abspath(PROGRAM_FILE) == @__FILE__
    for candidate in (isempty(ARGS) ? ["baseline", "localized"] : ARGS)
        source = joinpath(ROOT, "mixed-f1.0e6-gamma0.99", candidate)
        say("FINITE GAMMA EXACT TERMINAL SHIFT ", candidate)
        terminal_shift_control(source)
    end
    say("COMPLETE exact terminal-drive controls")
end
