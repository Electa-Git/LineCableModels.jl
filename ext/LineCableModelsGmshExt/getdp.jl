const FEM_FIELD_QUANTITIES = (
    "az", "b", "bm", "e", "ez", "em", "jz", "jm", "rhoj2"
)

function _pro_string(value::AbstractString)
    escaped = replace(String(value), '\\' => "\\\\", '"' => "\\\"")
    return "\"$escaped\""
end

function _pro_number(value::Real)
    return @sprintf("%.17g", Float64(value))
end

_pro_array(values) = "{" * join(_pro_number.(values), ", ") * "}"

function _write_problem_snapshot(path::String, problem::LineParametersProblem)
    open(path, "w") do io
        JSON3.pretty(io, ImportExport.serialize_value(problem))
        write(io, '\n')
    end
    return path
end

function _write_model_data(path::String, model::FEMResolvedModel)
    earth = model.problem.earth_props.layers[2]
    materials = model.material_plans
    open(path, "w") do io
        println(io, "// Resolved LineCableModels FEM inputs; solver logic lives in getdp/*.pro")
        println(io, "NumTerminals = ", length(model.terminal_ids), ";")
        println(io, "NumCables = ", length(model.problem.system.designs), ";")
        println(io, "NumMaterialRegions = ", length(materials), ";")
        println(io, "FrequencyCount = ", length(model.problem.frequencies), ";")
        println(io, "AIR_EM = ", model.tags.air, ";")
        println(io, "EARTH_EM = ", model.tags.earth, ";")
        println(io, "AIR_INF = ", model.tags.air_infinite, ";")
        println(io, "EARTH_INF = ", model.tags.earth_infinite, ";")
        println(io, "DOMAIN_INF = ", model.tags.infinite_domain, ";")
        println(io, "OUTBND_EM = ", model.tags.outer_boundary, ";")
        println(io, "OUTBND_ELE_INS = ", model.tags.outer_air_boundary, ";")
        println(io, "OUTBND_ELE_REF = ", model.tags.outer_earth_boundary, ";")
        println(io, "INTERFACE_AIR_SOIL = ", model.tags.interface, ";")
        println(io, "INNER_INF_BND = ", model.tags.inner_shell_boundary, ";")
        println(io, "TERMINAL = ", model.tags.terminal_base + 1, ";")
        println(io, "TERMINAL_CONTOUR = ", model.tags.terminal_contour_base + 1, ";")
        println(io, "CABLE_CONTOUR = ", model.tags.cable_contour_base + 1, ";")
        println(io, "TerminalNames() = Str[",
            join(_pro_string.(model.terminal_ids), ", "), "];")
        println(io, "MaterialRegionTags() = ",
            _pro_array(getproperty.(materials, :physical_tag)), ";")
        println(io, "MaterialIsConductor() = ",
            _pro_array([material.kind === :conductor for material in materials]), ";")
        println(io, "MaterialHasLoss() = ", _pro_array([
            any(value -> !iszero(real(value)), material.admittivity)
            for material in materials]), ";")
        println(io, "MaterialMu() = ",
            _pro_array([material.mu_r * 4π * 1e-7 for material in materials]), ";")
        for (index, material) in pairs(materials)
            println(io, "MaterialSigma_", index, "() = ",
                _pro_array(real.(material.admittivity)), ";")
            println(io, "MaterialEpsilon_", index, "() = ",
                _pro_array(imag.(material.admittivity) ./ (2π .* model.problem.frequencies)), ";")
        end
        println(io, "sigma_earth = ", _pro_number(inv(earth.rho)), ";")
        println(io, "eps_earth = ", _pro_number(earth.eps_r * 8.8541878128e-12), ";")
        println(io, "mu_earth = ", _pro_number(earth.mu_r * 4π * 1e-7), ";")
        println(io, "DomainRadius = ", _pro_number(model.domain_radius), ";")
        println(io, "ShellOuterRadius = ", _pro_number(model.shell_outer_radius), ";")
        println(io, "Xcenter = ", _pro_number(model.centre[1]), ";")
        println(io, "Ycenter = ", _pro_number(model.centre[2]), ";")
        println(io, "Zcenter = 0.0;")
    end
    return path
end

function _getdp_assets(root::AbstractString = joinpath(@__DIR__, "getdp"))
    return (
        model = joinpath(root, "model.pro"),
        jacobian = joinpath(root, "jacobian_integration.pro"),
        materials = joinpath(root, "materials.pro"),
        quasi_tem = joinpath(root, "quasi_tem.pro")
    )
end

# Capture solver text with the loaded Julia implementation, not halfway through
# a long mesh/solve sequence. Include dependencies invalidate the Julia cache.
const FEM_GETDP_SOURCES = map(_getdp_assets()) do path
    Base.include_dependency(path)
    read(path, String)
end

const GETDP_ARTIFACT_VERSION = v"3.5.0"
const GETDP_ENVIRONMENT_VARIABLE = "LINECABLEMODELS_GETDP"
const GETDP_ARTIFACTS_TOML = normpath(joinpath(@__DIR__, "..", "..", "Artifacts.toml"))

function _artifact_getdp()
    hash = artifact_hash("getdp", GETDP_ARTIFACTS_TOML)
    hash === nothing && return nothing
    root = artifact"getdp"
    version = string(GETDP_ARTIFACT_VERSION)
    relative = if Sys.isapple()
        joinpath("getdp-$version-MacOSX", "bin", "getdp")
    else
        joinpath("getdp-$version-Linux64", "bin", "getdp")
    end
    return (path=joinpath(root, relative), source=:artifact, artifact_hash=string(hash))
end

function _getdp_selection(formulation::LineCableModelsFEM; run_directory=nothing)
    explicit = formulation.execution.getdp_executable
    environment = get(ENV, GETDP_ENVIRONMENT_VARIABLE, nothing)
    selection = if explicit !== nothing
        (path=abspath(explicit), source=:explicit, artifact_hash=nothing)
    elseif environment !== nothing
        (path=abspath(environment), source=:environment, artifact_hash=nothing)
    else
        artifact = _artifact_getdp()
        if artifact !== nothing
            artifact
        else
            executable = Sys.which("getdp")
            executable === nothing ? nothing :
            (path=executable, source=:path, artifact_hash=nothing)
        end
    end
    selection !== nothing && isfile(selection.path) || _fem_error(
        :getdp, "GetDP", :getdp_executable,
        "GetDP executable is unavailable for $(triplet(HostPlatform())); " *
        "pass getdp_executable, set " * GETDP_ENVIRONMENT_VARIABLE *
        ", or install getdp on PATH when no package artifact supports this platform";
        run_directory)
    return merge(selection, (; path=realpath(selection.path)))
end

function _getdp_identity(executable::String)
    output = try
        buffer = IOBuffer()
        Base.run(pipeline(
            Cmd([executable, "-info"]), stdout = buffer, stderr = buffer
        ))
        String(take!(buffer))
    catch exception
        _fem_error(
            :getdp,
            "GetDP",
            :getdp_executable,
            "failed to execute $executable -info: $(sprint(showerror, exception))"
        )
    end
    occursin("getdp", lowercase(output)) &&
    occursin(r"\d+\.\d+", output) || _fem_error(
        :getdp,
        "GetDP",
        :getdp_executable,
        "executable identity check did not report GetDP: $executable"
    )
    return (sha256=bytes2hex(open(sha256, executable)), info=output)
end

function _resolve_getdp(formulation::LineCableModelsFEM, run::FEMRun)
    selection = _getdp_selection(formulation; run_directory=run.path)
    identity = _getdp_identity(selection.path)
    recorded = JSON3.read(read(joinpath(run.path, "input", "computation.json"), String))
    recorded_identity = Dict(String(key)=>value for (key,value) in pairs(recorded.getdp_identity))
    pop!(recorded_identity, "path", nothing) # schema 4 compatibility
    expected_identity = Dict(String(key)=>value for (key,value) in
        pairs(JSON3.read(JSON3.write(identity))))
    _resume_value_matches(recorded_identity, expected_identity) ||
        _fem_error(:getdp, "GetDP", :getdp_executable,
            "GetDP executable identity changed after input preparation; start a new run";
            run_directory=run.path)
    return selection.path
end

function _job_raw_paths(root::AbstractString, job_name::String)
    return (Z=joinpath(root, "raw", "jobs", "$job_name-Z.tsv"),
        P=joinpath(root, "raw", "jobs", "$job_name-P.tsv"))
end

_job_raw_paths(run::FEMRun, job_name::String) = _job_raw_paths(run.path, job_name)

function _valid_job_raw(
        path::String,
        terminal_count::Int,
        frequency_index::Int,
        frequency::Real,
        basis_terminal::Int
)
    isfile(path) || return false
    rows = readlines(path)
    length(rows) == terminal_count || return false
    responses = Set{Int}()
    for row in rows
        columns = split(row, '\t'; keepempty = true)
        length(columns) == 6 || return false
        parsed_frequency_index = tryparse(Int, columns[1])
        parsed_frequency = tryparse(Float64, columns[2])
        response = tryparse(Int, columns[3])
        parsed_basis = tryparse(Int, columns[4])
        real_part = tryparse(Float64, columns[5])
        imaginary_part = tryparse(Float64, columns[6])
        parsed_frequency_index == frequency_index || return false
        response isa Int && response in 1:terminal_count || return false
        parsed_basis == basis_terminal || return false
        parsed_frequency === nothing && return false
        isapprox(
            parsed_frequency,
            Float64(frequency);
            rtol = 16eps(Float64),
            atol = 0.0
        ) || return false
        real_part === nothing && return false
        imaginary_part === nothing && return false
        isfinite(real_part) && isfinite(imaginary_part) || return false
        response in responses && return false
        push!(responses, response)
    end
    return length(responses) == terminal_count
end

function _write_scan_completion!(run::FEMRun, model::FEMResolvedModel)
    frequency_count = length(model.problem.frequencies)
    terminal_count = length(model.terminal_ids)
    expected_rows = frequency_count * terminal_count * terminal_count
    path = joinpath(run.path, "raw", "scan_complete.tsv")
    temporary = tempname(dirname(path))
    open(temporary, "w") do io
        println(io, join(FEM_COMPLETE_HEADER, '\t'))
        println(io, join((
            frequency_count,
            terminal_count,
            expected_rows,
            expected_rows,
            expected_rows,
            1
        ), '\t'))
    end
    mv(temporary, path; force = true)
    return path
end

function _log_tail(records; count::Int = 20)
    isempty(records) && return "(no GetDP log records)"
    first_index = max(firstindex(records), lastindex(records) - count + 1)
    return join(@view(records[first_index:lastindex(records)]), '\n')
end
