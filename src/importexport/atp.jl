"""$(TYPEDSIGNATURES)

Export a [`LineCableSystem`](@ref) to an **ATPDraw‑compatible** XML file (LCC component with input data).

This routine serializes the cable system geometry (positions and outer radii) and the
already‑computed, frequency‑specific equivalent parameters of each cable component to the
ATPDraw XML schema. The result is written to disk and the absolute file path is returned
on success.

# Arguments

- `::Val{:atp}`: Backend selector for the ATP/ATPDraw exporter.
- `cable_system::LineCableSystem`: The system to export. Each entry in `cable_system.cables` provides one phase position and its associated [`CableDesign`](@ref). The number of phases exported equals `length(cable_system.cables)`.
- `earth_props::EarthModel`: Ground model used to populate ATP soil parameters. The exporter
uses the **last** layer’s base resistivity as *Grnd resis*.
- `base_freq::Number = f₀` \\[Hz\\]: System frequency written to ATP (`SysFreq`) and stored in component metadata. *This exporter does not recompute R/L/C/G; it writes the values as
present in the groups/components at the time of export.*
- `file_name::String = "*_export.xml"`: Output file name or path. If a relative path is given, it is resolved against the exporter’s source directory. The absolute path of the saved file is returned.

# Behavior

1. Create the ATPDraw `<project>` root and header and insert a single **LCC** component with
   `NumPhases = length(cable_system.cables)`.
2. For each [`CablePosition`](@ref) in `cable_system.cables`:

   * Write a `<cable>` element with:

     * `NumCond` = number of [`CableComponent`](@ref)s in the design,
     * `Rout` = outermost radius of the design (m),
     * `PosX`, `PosY` = cable coordinates (m).
3. For each [`CableComponent`](@ref) inside a cable:

   * Write one `<conductor>` element with fields (all per unit length):

     * `Rin`, `Rout` — from the component’s conductor group,
     * `rho` — conductor equivalence via [`calc_equivalent_rho`](@ref),
     * `muC` — conductor relative permeability via [`calc_equivalent_mu`](@ref),
     * `muI` — insulator relative permeability (taken from the first insulating layer’s material),
     * `epsI` — insulation relative permittivity via [`calc_equivalent_eps`](@ref),
     * `Cext`, `Gext` — shunt capacitance and conductance from the component’s insulator group.
4. Soil resistivity is written as *Grnd resis* using `earth_props.layers[end].base_rho_g`.
5. The XML is pretty‑printed and written to `file_name`. On I/O error, the function logs an error and returns `nothing`.

# Units

Units are printed in the XML file according to the ATPDraw specifications:

- Radii (`Rin`, `Rout`, `Rout` of cable): \\[m\\]
- Coordinates (`PosX`, `PosY`): \\[m\\]
- Length (`Length` tag): \\[m\\]
- Frequency (`SysFreq`/`Freq`): \\[Hz\\]
- Resistivity (`rho`, *Grnd resis*): \\[Ω·m\\]
- Relative permittivity (`epsI`) / permeability (`muC`, `muI`): \\[dimensionless\\]
- Shunt capacitance (`Cext`): \\[F/m\\]
- Shunt conductance (`Gext`): \\[S/m\\]

# Notes

* The exporter assumes each component’s equivalent parameters (R/G/C and derived ρ/ε/μ) were
  already computed by the design/group constructors at the operating conditions of interest.
* Mixed numeric types are supported; values are stringified for XML output. When using
  uncertainty types (e.g., `Measurements.Measurement`), the uncertainty is removed.
* Overlap checks between cables are enforced when building the system, not during export.

# Examples

```julia
# Build or load a system `sys` and an earth model `earth`
file = $(FUNCTIONNAME)(Val(:atp), sys, earth; base_freq = 50.0,
                       file_name = "system_id_export.xml")
println("Exported to: ", file)
```

# See also

* [`LineCableSystem`](@ref), [`CablePosition`](@ref), [`CableComponent`](@ref)
* [`EarthModel`](@ref)
* [`calc_equivalent_rho`](@ref), [`calc_equivalent_mu`](@ref), [`calc_equivalent_eps`](@ref)
  """
function export_data(::Val{:atp},
    cable_system::LineCableSystem,
    earth_props::EarthModel;
    base_freq=f₀,
    file_name::Union{String,Nothing}=nothing,
)::Union{String,Nothing}

    function _set_attributes!(element::EzXML.Node, attrs::Dict)
        for (k, v) in attrs
            element[k] = string(v)
        end
    end
    # --- 1. Setup Constants and Variables ---
    if isnothing(file_name)
        # caller didn't supply a name -> derive from cable_system if present
        file_name = joinpath(@__DIR__, "$(cable_system.system_id)_export.xml")
    else
        # caller supplied a path/name -> respect directory, but prepend system_id to basename
        requested = isabspath(file_name) ? file_name : joinpath(@__DIR__, file_name)
        if isnothing(cable_system)
            file_name = requested
        else
            dir = dirname(requested)
            base = basename(requested)
            file_name = joinpath(dir, "$(cable_system.system_id)_$base")
        end
    end

    num_phases = length(cable_system.cables)

    # Create XML Structure and LCC Component
    doc = XMLDocument()
    project = ElementNode("project")
    setroot!(doc, project)
    _set_attributes!(project, Dict("Application" => "ATPDraw", "Version" => "7.3", "VersionXML" => "1"))
    header = addelement!(project, "header")
    _set_attributes!(header, Dict("Timestep" => 1e-6, "Tmax" => 0.1, "XOPT" => 0, "COPT" => 0, "SysFreq" => base_freq, "TopLeftX" => 200, "TopLeftY" => 0))
    objects = addelement!(project, "objects")
    variables = addelement!(project, "variables")
    comp = addelement!(objects, "comp")
    _set_attributes!(comp, Dict("Name" => "LCC", "Id" => "$(cable_system.system_id)_1", "Capangl" => 90, "CapPosX" => -10, "CapPosY" => -25, "Caption" => ""))
    comp_content = addelement!(comp, "comp_content")
    _set_attributes!(comp_content, Dict("PosX" => 280, "PosY" => 360, "NumPhases" => num_phases, "Icon" => "default", "SinglePhaseIcon" => "true"))
    for side in ["IN", "OUT"]
        y0 = -20
        for k in 1:num_phases
            y0 += 10
            node = addelement!(comp_content, "node")
            _set_attributes!(node, Dict("Name" => "$side$k", "Value" => "C$(k)$(side=="IN" ? "SND" : "RCV")", "UserNamed" => "true", "Kind" => k, "PosX" => side == "IN" ? -20 : 20, "PosY" => y0, "NamePosX" => 0, "NamePosY" => 0))
        end
    end

    line_length = to_nominal(cable_system.line_length)
    soil_rho = to_nominal(earth_props.layers[end].base_rho_g)
    for (name, value) in [("Length", line_length), ("Freq", base_freq), ("Grnd resis", soil_rho)]
        data_node = addelement!(comp_content, "data")
        _set_attributes!(data_node, Dict("Name" => name, "Value" => value))
    end

    # Populate the LCC Sub-structure with CORRECTLY Structured Cable Data
    lcc_node = addelement!(comp, "LCC")
    _set_attributes!(lcc_node, Dict("NumPhases" => num_phases, "IconLength" => "true", "LineCablePipe" => 2, "ModelType" => 1))
    cable_header = addelement!(lcc_node, "cable_header")
    _set_attributes!(cable_header, Dict("InAirGrnd" => 1, "MatrixOutput" => "true", "ExtraCG" => "$(num_phases)"))

    for (k, cable) in enumerate(cable_system.cables)
        cable_node = addelement!(cable_header, "cable")

        num_components = length(cable.design_data.components)
        outermost_radius = to_nominal(cable.design_data.components[end].insulator_group.radius_ext)

        _set_attributes!(cable_node, Dict(
            "NumCond" => num_components,
            "Rout" => outermost_radius,
            "PosX" => to_nominal(cable.horz),
            "PosY" => to_nominal(cable.vert)
        ))

        for component in cable.design_data.components
            conductor_node = addelement!(cable_node, "conductor")

            cond_group = component.conductor_group
            cond_props = component.conductor_props
            ins_group = component.insulator_group
            ins_props = component.insulator_props

            rho_eq = (cond_props.rho)
            mu_r_cond = (cond_props.mu_r)
            mu_r_ins = (ins_props.mu_r)
            eps_eq = (ins_props.eps_r)

            _set_attributes!(conductor_node, Dict(
                "Rin" => to_nominal(cond_group.radius_in),
                "Rout" => to_nominal(cond_group.radius_ext),
                "rho" => to_nominal(rho_eq),
                "muC" => to_nominal(mu_r_cond),
                "muI" => to_nominal(mu_r_ins),
                "epsI" => to_nominal(eps_eq),
                "Cext" => to_nominal(ins_group.shunt_capacitance),
                "Gext" => to_nominal(ins_group.shunt_conductance)
            ))
        end
    end

    # Finalize and Write to File
    _set_attributes!(variables, Dict("NumSim" => 1, "IOPCVP" => 0, "UseParser" => "false"))

    try
        open(file_name, "w") do fid
            prettyprint(fid, doc)
        end
        @info "XML file saved to: $(display_path(file_name))"
        return file_name
    catch e
        @error "Failed to write XML file '$(display_path(file_name))'" exception = (e, catch_backtrace())
        return nothing
    end
end

# TODO: Develop `.lis` import and tests
# Issue URL: https://github.com/Electa-Git/LineCableModels.jl/issues/12
function read_data end
# I TEST THEREFORE I EXIST
# I DON´T TEST THEREFORE GO TO THE GARBAGE 
# """
#     read_atp_data(file_name::String, cable_system::LineCableSystem)

# Reads an ATP `.lis` output file, extracts the Ze and Zi matrices, and dynamically
# reorders them to a grouped-by-phase format based on the provided `cable_system`
# structure. It correctly handles systems with a variable number of components per cable.

# # Arguments
# - `file_name`: The path to the `.lis` file.
# - `cable_system`: The `LineCableSystem` object corresponding to the data in the file.

# # Returns
# - `Array{T, 2}`: A 2D complex matrix representing the total reordered series
#   impedance `Z = Ze + Zi` for a single frequency.
# - `nothing`: If the file cannot be found, parsed, or if the matrix dimensions in the
#   file do not match the provided `cable_system` structure.
# """
# function read_data(::Val{:atp},
#     cable_system::LineCableSystem,
#     freq::AbstractFloat;
#     file_name::String="$(cable_system.system_id)_1.lis"
# )::Union{Array{COMPLEXSCALAR,2},Nothing}
#     # --- Inner helper function to parse a matrix block from text lines ---
#     function parse_block(block_lines::Vector{String})
#         data_lines = filter(line -> !isempty(strip(line)), block_lines)
#         if isempty(data_lines)
#             return Matrix{ComplexF64}(undef, 0, 0)
#         end
#         matrix_size = length(split(data_lines[1]))
#         real_parts = zeros(Float64, matrix_size, matrix_size)
#         imag_parts = zeros(Float64, matrix_size, matrix_size)
#         row_counter = 1
#         for i in 1:2:length(data_lines)
#             if i + 1 > length(data_lines)
#                 break
#             end
#             real_line, imag_line = data_lines[i], data_lines[i+1]
#             try
#                 real_parts[row_counter, :] = [parse(Float64, s) for s in split(real_line)[1:matrix_size]]
#                 imag_parts[row_counter, :] = [parse(Float64, s) for s in split(imag_line)[1:matrix_size]]
#             catch e
#                 @error "Parsing failed" exception = (e, catch_backtrace())
#                 return nothing
#             end
#             row_counter += 1
#             if row_counter > matrix_size
#                 break
#             end
#         end
#         return real_parts + im * imag_parts
#     end

#     # --- Main Function Logic ---
#     if !isfile(file_name)
#         @error "File not found: $file_name"
#         return nothing
#     end
#     lines = readlines(file_name)
#     ze_start_idx = findfirst(occursin.("Earth impedance [Ze]", lines))
#     zi_start_idx = findfirst(occursin.("Conductor internal impedance [Zi]", lines))
#     if isnothing(ze_start_idx) || isnothing(zi_start_idx)
#         @error "Could not find Ze/Zi headers."
#         return nothing
#     end

#     Ze = parse_block(lines[ze_start_idx+1:zi_start_idx-1])
#     Zi = parse_block(lines[zi_start_idx+1:end])
#     if isnothing(Ze) || isnothing(Zi)
#         return nothing
#     end

#     # --- DYNAMICALLY GENERATE PERMUTATION INDICES (Numerical Method) ---
#     component_counts = [length(c.design_data.components) for c in cable_system.cables]
#     total_conductors = sum(component_counts)
#     num_phases = length(component_counts)
#     max_components = isempty(component_counts) ? 0 : maximum(component_counts)

#     if size(Ze, 1) != total_conductors
#         @error "Matrix size from file ($(size(Ze,1))x$(size(Ze,1))) does not match total components in cable_system ($total_conductors)."
#         return nothing
#     end

#     num_conductors_per_type = [sum(c >= i for c in component_counts) for i in 1:max_components]
#     type_offsets = cumsum([0; num_conductors_per_type[1:end-1]])

#     permutation_indices = Int[]
#     sizehint!(permutation_indices, total_conductors)
#     instance_counters = ones(Int, max_components)
#     for phase_idx in 1:num_phases
#         for comp_type_idx in 1:component_counts[phase_idx]
#             instance = instance_counters[comp_type_idx]
#             original_idx = type_offsets[comp_type_idx] + instance
#             push!(permutation_indices, original_idx)
#             instance_counters[comp_type_idx] += 1
#         end
#     end

#     Ze_reordered = Ze[permutation_indices, permutation_indices]
#     Zi_reordered = Zi[permutation_indices, permutation_indices]

#     return Ze_reordered + Zi_reordered
# end


"""$(TYPEDSIGNATURES)

Export calculated [`LineParameters`](@ref) (series impedance **Z** and shunt admittance **Y**) to an **compliant** `ZY` XML file.

This routine writes the complex **Z** and **Y** matrices versus frequency into a compact XML
structure understood by external tools. Rows are emitted as comma‑separated complex entries
(`R+Xi` / `G+Bi`) with one `<Z>`/`<Y>` block per frequency sample.

# Arguments

- `::Val{:atp}`: Backend selector for the ATP/ATPDraw ZY exporter.
- `line_params::LineParameters`: Object holding the frequency‑dependent matrices `Z[:,:,k]` and `Y[:,:,k]`.
- `freq::AbstractVector{<:AbstractFloat}`: Frequency vector \\[Hz\\] aligned with the third dimension of `Z` and `Y`.
- `file_name::String = "ZY_export.xml"`: Output file name or path. If relative, it is resolved against the exporter’s source directory. The absolute path of the saved file is returned.
- `cable_system::Union{LineCableSystem,Nothing} = nothing`: Optional system used only to derive a default name. When provided and `file_name` is not overridden, the exporter uses `"\$(cable_system.system_id)_ZY_export.xml"`.

# Behavior

1. The root tag `<ZY>` includes `NumPhases`, `Length` (fixed to `1.0`), and format attributes `ZFmt="R+Xi"`, `YFmt="G+Bi"`.
2. For each frequency `fᵏ = freq[k]`:

   * Emit a `<Z Freq=...>` block with `num_phases` lines, each line the `k`‑th slice of row `i` formatted as `real(Z[i,j,k]) + imag(Z[i,j,k])i`.
   * Emit a `<Y Freq=...>` block in the same fashion (default `G+Bi`).
3. Close the `</ZY>` element and write to disk. On I/O error the function logs and returns `nothing`.

# Units

Units are printed in the XML file according to the ATPDraw specifications:

- `freq` (XML `Freq` attribute): \\[Hz\\]
- `Z` entries: \\[Ω/km\\] (per unit length)
- `Y` entries: \\[S/km\\] (per unit length) when `YFmt = "G+Bi"`
- XML `Length` attribute: \\[m\\]

# Notes

- The exporter assumes `size(line_params.Z, 1) == size(line_params.Z, 2) == size(line_params.Y, 1) == size(line_params.Y, 2)` and `length(freq) == size(Z,3) == size(Y,3)`.
- Numeric types are stringified; mixed numeric backends (e.g., with uncertainties) are acceptable as long as they can be printed via `@sprintf`.
- This exporter **does not** modify or recompute matrices; it serializes exactly what is in `line_params`.

# Examples

```julia
# Z, Y have already been computed into `lp::LineParameters`
file = $(FUNCTIONNAME)(:atp, lp, freq; file_name = "ZY_export.xml")
println("Exported ZY to: ", file)

# Naming based on a cable system
file2 = $(FUNCTIONNAME)(:atp, lp, freq; cable_system = sys)
println("Exported ZY to: ", file2)  # => "\$(sys.system_id)_ZY_export.xml"
```

# See also

* [`LineParameters`](@ref)
* [`LineCableSystem`](@ref)
* [`export_data(::Val{:atp}, cable_system, ...)`](@ref) — exporter that writes full LCC input data
  """
function export_data(::Val{:atp},
    line_params::LineParameters,
    freq::AbstractVector{<:AbstractFloat};
    file_name::Union{String,Nothing}=nothing,
    cable_system::Union{LineCableSystem,Nothing}=nothing
)::Union{String,Nothing}

    # Resolve final file_name while preserving any user-supplied path.
    if isnothing(file_name)
        # caller didn't supply a name -> derive from cable_system if present
        if isnothing(cable_system)
            file_name = joinpath(@__DIR__, "ZY_export.xml")
        else
            file_name = joinpath(@__DIR__, "$(cable_system.system_id)_ZY_export.xml")
        end
    else
        # caller supplied a path/name -> respect directory, but prepend system_id to basename if cable_system provided
        requested = isabspath(file_name) ? file_name : joinpath(@__DIR__, file_name)
        if isnothing(cable_system)
            file_name = requested
        else
            dir = dirname(requested)
            base = basename(requested)
            file_name = joinpath(dir, "$(cable_system.system_id)_$base")
        end
    end

    @info ("ZY export called",
        :method => "ZY",
        :cable_system_isnothing => isnothing(cable_system),
        :cable_system_type => (isnothing(cable_system) ? :nothing : typeof(cable_system)),
        :file_name_in => file_name)

    cable_length = isnothing(cable_system) ? 1.0 : to_nominal(cable_system.line_length)
    atp_format = "G+Bi"
    # file_name = isabspath(file_name) ? file_name : joinpath(@__DIR__, file_name)

    open(file_name, "w") do fid
        num_phases = size(line_params.Z, 1)
        y_fmt = (atp_format == "C") ? "C" : "G+Bi"

        @printf(fid, "<ZY NumPhases=\"%d\" Length=\"%.4f\" ZFmt=\"R+Xi\" YFmt=\"%s\">\n", num_phases, cable_length, y_fmt)

        # --- Z Matrix Printing ---
        for (k, freq_val) in enumerate(freq)
            @printf(fid, "  <Z Freq=\"%.16E\">\n", to_nominal(freq_val))
            for i in 1:num_phases
                row_str = join([@sprintf("%.16E%+.16Ei", to_nominal(real(line_params.Z[i, j, k])), to_nominal(imag(line_params.Z[i, j, k]))) for j in 1:num_phases], ",")
                println(fid, row_str)
            end
            @printf(fid, "  </Z>\n")
        end

        # --- Y Matrix Printing ---
        if atp_format == "C"
            freq1 = to_nominal(f[1])
            @printf(fid, "  <Y Freq=\"%.16E\">\n", freq1)
            for i in 1:num_phases
                row_str = join([@sprintf("%.16E", to_nominal(imag(line_params.Y[i, j, 1]) / (2 * pi * freq1))) for j in 1:num_phases], ",")
                println(fid, row_str)
            end
            @printf(fid, "  </Y>\n")
        else # Case for "G+Bi"
            for (k, freq_val) in enumerate(freq)
                @printf(fid, "  <Y Freq=\"%.16E\">\n", to_nominal(freq_val))
                for i in 1:num_phases
                    row_str = join([@sprintf("%.16E%+.16Ei", to_nominal(real(line_params.Y[i, j, k])), to_nominal(imag(line_params.Y[i, j, k]))) for j in 1:num_phases], ",")
                    println(fid, row_str)
                end
                @printf(fid, "  </Y>\n")
            end
        end

        # --- Footer ---
        println(fid, "</ZY>")
    end
    try
        # Use pretty print option for debugging comparisons if needed
        # open(filename, "w") do io; prettyprint(io, doc); end
        if isfile(file_name)
            @info "XML file saved to: $(display_path(file_name))"
        end
        return file_name
    catch e
        @error "Failed to write XML file '$(display_path(file_name))': $(e)"
        isa(e, SystemError) && println("SystemError details: ", e.extrainfo)
        return nothing
        rethrow(e) # Rethrow to indicate failure clearly
    end
end